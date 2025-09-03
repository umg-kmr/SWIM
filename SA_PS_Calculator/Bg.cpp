#include <iostream> 
#include <cmath> 
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/tools/roots.hpp>
#include <stdexcept>
#include <algorithm>    // std::reverse, replace_if
#include <cctype>    // for std::isspace

using namespace std;
using namespace boost::numeric::odeint;

vector <double> klist;
vector <double> Plist;
vector<double> Qfile;
vector<double> Gfile;
function<double(double)> logGasQ;

double Nend = 0.0;

double As_val = 0.0;
double ns_val = 0.0;
double Nreh = 0.0;

uintmax_t max_iter = 1000000; //upper limit for root finding algorithm
double maxiter_bg = 1e6; //Upper limit on the background integration iterations
double maxiter_dxdN = 1e5; //Upper limit on the dxdN integration iterations
double atol_bg = 1e-14; //Sets absolute tolerance for all the integrators
double rtol_bg = 1e-12; //Sets relative tolerance for all the integrators
double tend = 200.0; //Upper limit of Bg integration

//Function to terminate root finding algorithm with some epsilon.
struct root_stop  {
    bool operator() (double r1, double r2)  {
         return abs(r1 - r2) <= 1e-14 + 1e-14*max(abs(r1),abs(r2));
    }
};

//4 - point derivative approximation
double num_derv (function<double(double)> F, double h, double x) {
    return ((-11.0*F(x))/6.0 + 3.0*F(h + x) - (3.0*F(2.0*h + x))/2.0 + F(3.0*h + x)/3.0)/h;
}

//Passing model functions by reference as the background code doesn't modify the original model functions like potential and upsilon.
void bg_solver (const function<double(double)> &V,const function<double(double)> &Vd,const function<double(double,double)> &Ups,const function<double(double)> &GQ,double Cr,double Np,double phi_ini,double php_ini,double T_ini,int therm,double kp,double klow,double kup,int npts,int camb_mcmc, int want_Np_autocalc, int verbose) {
    double Qstar = 0.0;
    typedef boost::array< double , 3 > state_type; //For boost ode solver, number of differential equations (3)
    
    //Arrays to store background calculations
    vector <double> Nl;
    vector <double> phil;
    vector <double> phpl;
    vector <double> Tl;

    //Vectors to store reheating calculations
    vector <double> Nrh;
    vector <double> phirh;
    vector <double> phprh;
    vector <double> Trh;

    //dxdN array solutions
    vector <double> x_back;
    vector <double> N_back;
    vector <double> x_forw;
    vector <double> N_forw;
    vector <double> x_final;
    vector <double> N_final;

    int cntr_reject_initial = 0; //Reject storing initial condition


    using boost::math::interpolators::makima;
    //Initializer for interpolation
    function<double(double)> phiasN;
    function<double(double)> phpasN;
    function<double(double)> TasN;
    function<double(double)> NasX;

    double Npp = 0.0; //Npivot dummy variable
   
    auto gen_klist = [klow,kup,npts] (vector <double> &klist) -> void {
        vector <double> temp;
        double step = (kup-klow)/(npts-1.0);
        for (double i=klow; i<=kup; i+=step) {
            temp.push_back(i);
        }

        //Raise to power of 10
        for (unsigned long k = 0; k<temp.size(); k++) {
            klist.push_back(pow(10.0,temp[k]));
        }
    };
      
    //Generates the klist between, 10^kup and k^low of size npts
    gen_klist(klist);
 
    auto H = [V,Cr] (double phi,double phip,double T) -> double {
        return sqrt(2.0 * ( V(phi) + (Cr*pow(T,4.0))) *  pow( ( 6.0 - (phip*phip) ) ,-1.0)  );
    };

    auto Hp = [H,Cr] (double phi,double phip,double T) -> double {
        double H1;
        H1 = H(phi,phip,T);
        return - ( (H1*(phip*phip))/2.0 ) - ( (2.0/3.0) * ((Cr*pow(T,4.0))/H1) );
    };

    auto eH = [H,Hp] (double phi,double phip,double T) -> double {
        return -Hp(phi,phip,T)/H(phi,phip,T);
    };

    auto Q = [Ups,H] (double phi,double phip,double T) -> double {
        return Ups(phi,T)/(3.0*H(phi,phip,T));
    };
    
 
    auto func = [Cr,Ups,Vd,H,Hp] ( const state_type &x ,state_type &dxdt , double t ) -> void{
         
        double phi;
        double php;
        double T;
        double Hpi;
        double Hi;
        double Upsi;
        double Vdi;
        
        phi = x[0];
        php = x[1];
        T = x[2];

        Hi = H(phi,php,T);
        Hpi = Hp(phi,php,T);
        Upsi = Ups(phi,T);
        Vdi = Vd(phi);

        
        dxdt[0] = php; 
        dxdt[1] = -(php*( (Hpi/Hi) + 3.0 + (Upsi/Hi))) - ( Vdi/(Hi*Hi));
        dxdt[2] = -T + ((php*php) * ( (Upsi*Hi)/(4.0*Cr*(T*T*T)) ) );
    };
    
    //ODE Observer to view and store results 
    auto obsv = [eH,&Nl,&phil,&phpl,&Tl] ( const state_type &x ,const double t) -> void {
        if (eH(x[0],x[1],x[2])>1.0 || Nl.size()>maxiter_bg) { //caps upper limit on maximum iterations
          throw runtime_error("Inflation Ended"); //Signals the end of Inflation
        }
        //cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] <<" eH "<<eH(x[0],x[1],x[2])<< endl;
        Nl.push_back(t);
        phil.push_back(x[0]);
        phpl.push_back(x[1]);
        Tl.push_back(x[2]);
    };

    //ODE solver
    auto SolveODE = [phi_ini,php_ini,T_ini,func,obsv] () -> void {
        auto stepper = make_controlled( atol_bg , rtol_bg , runge_kutta_fehlberg78 <state_type>() );
        state_type x = { phi_ini,php_ini,T_ini }; // initial conditions
        try {
            integrate_adaptive( stepper ,func , x , 0.0 , tend, 1e-6 , obsv ); 
        }
        catch (const exception& e) {
           // cout<<"Threshold reached, stopping integration!"<<endl;
        }
    };

    //Interpolation Initializer
    auto get_Interp = [&phiasN,&phpasN,&TasN] (vector<double> Nlcop, vector<double> philcop, vector<double> phplcop, vector<double> Tlcop) -> void {
        vector <double> Nlcopy = Nlcop; //Need copy as the move function changes the original vector to null
        phiasN = makima(std::move(Nlcopy),std::move(philcop));
        Nlcopy = Nlcop;                //Reupdating the moved vector
        phpasN = makima(std::move(Nlcopy),std::move(phplcop));
        Nlcopy = Nlcop;
        TasN = makima(std::move(Nlcopy),std::move(Tlcop));
    };

    SolveODE(); //Call solver to initialize the bg vectors
    //cout<<Nl.size()<<endl;
    
    //Error handling
    if (Nl.size()<1)  {
        if (verbose == 1) {
            cout<<"Solver did not run"<<endl;
        }
        Plist.assign(npts,1); //Create a P vector of size npts filled with 1s.
        return;
    }

    Nend = Nl.back(); //Stores the value of N at the end of Inflation ; eH = 1.
		      //
    if ( (Nend>=tend) || (Nl.size()<=4) || ( (Nl.size()>=maxiter_bg ) && (eH(phil.back(),phpl.back(),Tl.back())<0.9999) ) )   { //Early exit if Upper limit of Bg integration is reached or number of datapoints are too less for interpolation
        if (verbose == 1) {
            cout<<"Exiting: Either not enough points for interpolation or Nend has exceeded the cap of 1000 or upper limit of integration steps reached."<<endl;
        }
        Plist.assign(npts,1);
        return;
    }

    //ODE section for reheating -> eH = 2 or w_eff =1/3
    auto w_eff = [V,H,Cr] (double phi,double phip,double T) -> double {
        double Hi = H(phi,phip,T);
        double Vi = V(phi);
        double p_tot = ((1.0/3.0)*Cr*T*T*T*T) + (Hi*phip*Hi*phip*0.5) - Vi;
        double rho_tot = (Cr*T*T*T*T) + (Hi*phip*Hi*phip*0.5) + Vi; 
        return p_tot/rho_tot;
    };
    //ODE Observer to view and store results 
    auto obsv_reh = [eH,&Nrh,&phirh,&phprh,&Trh] ( const state_type &x ,const double t) -> void {
        if (eH(x[0],x[1],x[2])>2.0 || Nrh.size()>maxiter_bg) { //caps upper limit on maximum iterations
          throw runtime_error("Reheating."); //Signals the end of Inflation
        }
        //cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] <<" eH "<<eH(x[0],x[1],x[2])<< endl;
        Nrh.push_back(t);
        phirh.push_back(x[0]);
        phprh.push_back(x[1]);
        Trh.push_back(x[2]);
    };

    //ODE solver for reheating
    auto SolveODE_reh = [&Nl,&phil,&phpl,&Tl,func,obsv_reh] () -> void {
        auto stepper = make_controlled( atol_bg , rtol_bg , runge_kutta_fehlberg78 <state_type>() );
        state_type x = { phil.back(),phpl.back(),Tl.back() }; // initial conditions
        try {
            integrate_adaptive( stepper ,func , x , Nend , tend, 1e-6 , obsv_reh ); 
        }
        catch (const exception& e) {
           // cout<<"Threshold reached, stopping integration!"<<endl;
        }
    };
    SolveODE_reh();
    Nreh = Nrh.back();
    //cout<<"Nreh: "<<Nreh<<endl;
     if ( (Nreh>=tend) || ( (Nrh.size()>=maxiter_bg ) && ( w_eff(phirh.back(),phprh.back(),Trh.back())<=0.3) ) )    { //Early exit if Upper limit of Bg integration is reached.
        if (verbose == 1) {
            cout<<"Exiting: Error in calculation of e-folds at reheating."<<endl;
        }
        return;
    }

    auto Calc_Nstar  = [H,kp] (auto phiasN,auto phpasN,auto TasN, auto Trh, double Nend, double Nreh) -> double {

        /*make a partial (lambda) function that only takes the root to be found and captures the other objects that are passed to the main function that can change during runtime
        like interpolation objects*/
        auto Nstar = [H,kp,Nend,phiasN,phpasN,TasN,Trh,Nreh] (double Ntry) -> double {
            double Htry = H(phiasN(Ntry),phpasN(Ntry),TasN(Ntry));
            return (kp*2.63e-57) - (exp(Ntry-Nreh)*pow( (43.0/(11.0*(106.75))),(1.0/3.0) )*((2.349e-13*4.11e-19)/(Trh.back()))*(Htry));
        };
    	auto res = boost::math::tools::bisect(Nstar, 0.0, Nend,root_stop(),max_iter);
    	return (res.second + res.first)/2;
    };

    
    try {
        get_Interp(Nl,phil,phpl,Tl); //Initialize interpolations after succesful computation of background, Np calc should follow this always.
    }

    catch (const exception& e) {
        if (verbose == 1) {
            cout<<"Error with interpolation in background"<<endl;
        }
        Plist.assign(npts,1);
        //Clear out bg integration data
        Nl = vector<double> (); 
        phil = vector<double> (); 
        phpl = vector<double> ();
        Tl = vector<double> (); 
        return;
    }
    
    //Clear out bg integration data
    Nl = vector<double> (); 
    phil = vector<double> (); 
    phpl = vector<double> ();
    Tl = vector<double> (); 

    //Check if user wants Np to be auto calculated or passed through model_calc.c
    if (want_Np_autocalc == 1) {
        try {
            Npp = Calc_Nstar(phiasN,phpasN,TasN,Trh,Nend,Nreh);
        }
	catch (const exception& e) {
            if (verbose == 1) {
                cout<<"Error while calculating N_pivot"<<endl;
            }
            return;
        }

    }
    else {
        Npp = Np;
    }
    //Constrain duration of Inflation
    try {
        if ( ((Nend-Npp) < 42.0) || (Npp<3.0)  )   {
            if (verbose == 1) {
                cout<<"Duration of Inflation less than 42.0"<<endl;
            }
            return;
        }
    }
    catch (const exception& e) {
        return;
    }
    //cout<<"Np: "<<Npp<<endl;
    /*Qstar = Q(phiasN(Npp),phpasN(Npp),TasN(Npp));
    if ( (Qstar>7.19) || (Qstar<5.77)  ) {
        return;
    }*/ 
    /*    Bg calculations done    */

    //Analytical Power Spectrum

    typedef boost::array< double , 1 > state_type2; //For boost ode solver, number of differential equations (1)
    //dxdN function ode where x = ln(k/kp)
    auto dxdN = [eH,phiasN,phpasN,TasN] (const state_type2 &x ,state_type2 &dxdt , double t) -> void{
        dxdt[0] = 1.0-eH(phiasN(t),phpasN(t),TasN(t));
    };

    auto obs_back = [&N_back,&x_back] ( const state_type2 &x ,const double t) -> void {
        if (N_back.size()>maxiter_dxdN) {
          throw runtime_error("Maxiter exceeded.");
        }
        //cout << t << '\t' << x[0]<< endl;
        N_back.push_back(t);
        x_back.push_back(x[0]);
    };

    auto obs_forw = [&N_forw,&x_forw,&cntr_reject_initial] ( const state_type2 &x ,const double t) -> void {
         if (N_forw.size()>maxiter_dxdN) {
            throw runtime_error("Maxiter exceeded.");
         }
        //cout << t << '\t' << x[0]<< endl;

        //reject initial condition : don't want duplicates in the final vector
        if (cntr_reject_initial !=0 ) {
            N_forw.push_back(t);
            x_forw.push_back(x[0]);
        }
        else {
            cntr_reject_initial+=1;
        }
    };

    //cout<<"Np: "<<Npp<<endl;
    auto solve_dxdN_back =[dxdN,Npp,obs_back] () -> void {
        auto stepper = make_controlled( atol_bg , rtol_bg, runge_kutta_fehlberg78 <state_type2> () );
        state_type2 x = {0.0}; // initial conditions
        integrate_adaptive( stepper ,dxdN , x , Npp , 0.0 , -1e-6 , obs_back );
    };

    auto solve_dxdN_forw =[dxdN,Npp,obs_forw] (double Nend) -> void {
        auto stepper = make_controlled( atol_bg , rtol_bg , runge_kutta_cash_karp54 <state_type2> () );
        state_type2 x = {0.0}; // initial conditions
        integrate_adaptive( stepper ,dxdN , x , Npp , Nend, 1e-6 , obs_forw );
    };
    
    try {
    	solve_dxdN_back();
        solve_dxdN_forw(Nend);
    }

    catch (const exception& e) {
        if (verbose == 1) {
            cout<<"Error with dxdN computation"<<endl;
        }
        Plist.assign(npts,1);
        N_back = vector<double> (); 
        x_back = vector<double> (); 
        N_forw = vector<double> ();
        x_forw = vector<double> (); 
        return;
    }


    auto final_x_N = [&N_final,&x_final,&N_forw,&x_forw,&N_back,&x_back] () -> void{
        //Reverse the backward integration
        std::reverse(N_back.begin(),N_back.end());
        std::reverse(x_back.begin(),x_back.end());

        //Store the backwards content
        N_final = std::move(N_back);
        x_final = std::move(x_back);

        //Store the forwards content
        N_final.insert(N_final.end(), make_move_iterator(N_forw.begin()),make_move_iterator( N_forw.end()));
        x_final.insert(x_final.end(), make_move_iterator(x_forw.begin()), make_move_iterator(x_forw.end()));

    };

    final_x_N();
    //Remove intermediate arrays
    N_back = vector<double> (); 
    x_back = vector<double> (); 
    N_forw = vector<double> ();
    x_forw = vector<double> (); 

    //return if not enough points for interpolation
    if (N_final.size()<=4) {
        if (verbose == 1) {
            cout<<"Not enough points for interpolation in NasX"<<endl;
        }
        Plist.assign(npts,1);

        return;
    }

    //Remove last elements from final vectors just for sanity for interpolation 
    x_final.pop_back();
    N_final.pop_back();
    
    //Interpolation Initializer for NasX
    try {
        NasX = makima(std::move(x_final),std::move(N_final));
    }

    catch (const exception& e) {
         if (verbose == 1) {
            cout<<"Error with interpolation in NasX"<<endl;
         }
         Plist.assign(npts,1);
         N_final = vector<double> ();
         x_final = vector<double> ();
         return;
    }

    //Remove final dxdN vectors
    N_final = vector<double> ();
    x_final = vector<double> (); 
    
    auto Nask = [NasX,kp] (double k) -> double {
        double X;
        X = log(k/kp);
        return NasX(X);
    };

    //Analytical Power Spectrum
    auto P = [Nask,phiasN,phpasN,TasN,H,Q,therm,GQ] (double k) -> double {
        double NN = Nask(k);
        double phin = phiasN(NN);
        double phpn = phpasN(NN);
        double Tn = TasN(NN);
        double Hn = H(phin,phpn,Tn);
        double Qn = Q(phin,phpn,Tn);
        double distrib = 0.0;
        if (therm == 1) {
            distrib = 1/tanh(Hn/(2*Tn));
        }
        else if (therm == 0) {
            distrib = 1.0;
        }
        return pow((Hn/(2.0*M_PI*phpn)),2.0) * (distrib + (Tn/Hn)*( 2.0*sqrt(3.0)*M_PI*Qn/sqrt(3.0+4.0*M_PI*Qn) )) * GQ(Qn);
    };


    auto PT = [Nask,phiasN,phpasN,TasN,H] (double k) -> double {
        double NN = Nask(k);
        double phin = phiasN(NN);
        double phpn = phpasN(NN);
        double Tn = TasN(NN);
        double Hn = H(phin,phpn,Tn);
        return (2.0*pow(Hn,2.0))/pow(M_PI,2.0);
    };
 
    auto ns = [P] (double k) -> double {
        double h = 1e-4;
        return 1.0 + ( (k/P(k)) * num_derv(P,h,k) );
    };

    auto r = [PT,P] (double k) -> double {
        return PT(k)/P(k);
    };
    
    //Constraint on r
    /*try {
        if (r(0.002)>0.056) {
            return;
        }   
    }
    catch (const exception& e) {
         if (verbose == 1) {
             cout<<"Tensor-to-Scalar ratio (r) could not be calculated, check for interpolation error."<<endl;
         }
         Plist.assign(npts,1);
         return;
    }*/

    

    auto alphs = [ns] (double k) -> double {
        double h = 1e-4;
        return k * (  num_derv(ns,h,k)  );
    };

    auto Calc_Plist = [P] (vector <double> &klist, vector <double> &Plist) -> void {
        //#pragma omp parallel for    //reduces performance if enabled, might be useful for very large value of npts; seldom required
        for (unsigned long i=0;i<klist.size();i++) {
            Plist.push_back(P(klist[i]));
        }

    };

    if (camb_mcmc == 1) {
        try {
            Calc_Plist(klist,Plist);
        }

        catch (const exception& e) {
            if (verbose == 1) {
                cout<<"Error in computing power spectrum for the full range of k values"<<endl;
            }
            Plist.assign(npts,1);
        }
    }

    else if (camb_mcmc == 0) {
        try {
            As_val = P(kp);
            ns_val = ns(kp);
        }
        catch (const exception& e) {
            if (verbose == 1) {
                cout<<"Error in computing power spectrum and/or spectral tilt at the specified k"<<endl;
            }
            Plist.assign(npts,1);

        }
    }
    else {
        cout<<"Please check the camb_mcmc setting it can be either 1 or 0"<<endl;
    }

}


void GQ_parser (string fname) {
    ifstream file(fname); //Read the file
    string line;
    string cell;
   
    //Clear out old entries (if any)
    Qfile = vector<double> ();
    Gfile = vector<double> ();

    while (getline(file, line)) {
        replace_if(line.begin(), line.end(), ::isspace, ','); //Convert spaces to comma
        stringstream ss(line);
        getline(ss, cell, ','); //Split the string on the bases of comma
        Qfile.push_back(stod(cell)); //Convert string to double and store
        getline(ss, cell, ',');
        Gfile.push_back(stod(cell));
    }
    file.close();
}

void GQ_interp () {
    using boost::math::interpolators::makima;
    vector<double> logQ(Qfile.size());
    vector<double> logG(Gfile.size());
    //log - transform data
    transform(Qfile.begin(), Qfile.end(), logQ.begin(), [](double x) {
        return std::log10(x);  
    });
    transform(Gfile.begin(), Gfile.end(), logG.begin(), [](double y) {
        return std::log10(y);  
    });
    logGasQ = makima(std::move(logQ),std::move(logG));
}

