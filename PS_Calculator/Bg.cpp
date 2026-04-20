#include <iostream> 
#include <fstream>
#include <cmath> 
#include <vector>
#include <string>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/tools/roots.hpp>
#include <stdexcept>
#include <random>
#include <chrono> // Include for time-based seeding
#include <omp.h>  //For parallel processing
#include <unistd.h> // For getpid()

using namespace std;
using namespace boost::numeric::odeint;

//generating random seeds using random device
chrono::high_resolution_clock clk;
random_device rd;
//RNG engine for thread safety defined with thread_local
static thread_local mt19937_64 rngen(clk.now().time_since_epoch().count() ^ rd() ^ omp_get_thread_num() ^ getpid());

//Global vectors to store power spectrum
vector <double> klist;
vector <double> Plist;
//Global vectors to store background
vector <double> N_glob_array;
vector <double> phi_glob_array;
vector <double> php_glob_array;
vector <double> T_glob_array;

double Nend = 0.0;
double tend = 200.0; //Upper limit of Bg integration
double maxiter_bg = 1e6; //Upper limit on the background integration iterations
uintmax_t max_iter = 1000000; //upper limit for root finding algorithm
double Nreh = 0.0;
double atol_bg = 1e-14; //Sets absolute tolerance for all the integrators
double rtol_bg = 1e-12; //Sets relative tolerance for all the integrators

double PT_kp = 0.0;

//Function to terminate root finding algorithm with some epsilon.
struct root_stop  {
    bool operator() (double r1, double r2)  {
         return abs(r1 - r2) <= 1e-16 + 1e-14*max(abs(r1),abs(r2));
    }
};

//Passing model functions by reference as the background code doesn't modify the original model functions like potential and upsilon.
void bg_solver (const function<double(double)> &V, const function<double(double)>& Vd,const function<double(double)>& Vdd,const function<double(double,double)>& Ups, const function<double(double,double)>& pT_Ups,const function<double(double,double)>& pph_Ups,double Cr,double Np,double phi_ini,double php_ini, double T_ini,int therm, double kp, double klow, double kup, double EM_step, int npts, int Nrealz, int want_Np_autocalc, int verbose,int rad_noise) {


    typedef boost::array< double , 3 > state_type; //For boost ode solver, number of differential equations (3)
    
    //Vectors to store background calculations
    vector <double> Nl;
    vector <double> phil;
    vector <double> phpl;
    vector <double> Tl;
    
    //Vectors to store reheating calculations
    vector <double> Nrh;
    vector <double> phirh;
    vector <double> phprh;
    vector <double> Trh;

    //Initializer for interpolation
    function<double(double)> phiasN;
    function<double(double)> phpasN;
    function<double(double)> TasN;

    using boost::math::interpolators::makima;

    double a0 = 0.0; //Scale factor normalization
    double Npp = 0.0; //Npivot dummy variable
    
    //Generates the klist when called
    auto gen_klist = [kup,klow,npts] (vector <double> &klist) -> void {
        vector <double> temp;
        double step = (kup-klow)/(npts-1.0);
        double kint = klow;

        //Loop from klow to kup in step size calculated above
        for (int i=0; i<npts; i++) {
            temp.push_back(kint);
            kint+=step;
        }

        //Raise to power of 10
        for (int k = 0; k<npts; k++) {
            klist.push_back(exp10(temp[k]));
        }
    };

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
        dxdt[2] = -T + ( (php*php) * ( (Upsi*Hi)/(4.0*Cr*(T*T*T)) ) );
    };
    
    //ODE Observer to view and store results 
    auto obsv = [eH,&Nl,&phil,&phpl,&Tl] ( const state_type &x ,const double t) -> void {
        if (eH(x[0],x[1],x[2])>1.0 || Nl.size()>maxiter_bg) { //prevents out of bounds
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
        auto stepper = make_controlled( atol_bg , rtol_bg , runge_kutta_fehlberg78 < state_type >() );
        state_type x = { phi_ini,php_ini,T_ini }; // initial conditions
        try {
            integrate_adaptive( stepper ,func , x , 0.0 , tend, 1e-6 , obsv ); //make_dense_output<runge_kutta_dopri5< state_type >>( 0.0 , 1.0e-14 )
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
     //Error handling
    if (Nl.size()<1) {
        if (verbose == 1) {
            cout<<"Solver did not run!"<<endl;
        }
        Plist.assign(npts,1); //Create a P vector of size npts filled with 1s.
        return;
    }

    Nend = Nl.back();
    //Copy background integration to global vectors
    N_glob_array = Nl;
    phi_glob_array = phil;
    php_glob_array = phpl;
    T_glob_array = Tl;
    
    if ( (Nend>=tend) || (Nl.size()<=4) || ( (Nl.size()>=maxiter_bg ) && (eH(phil.back(),phpl.back(),Tl.back())<0.9999) ) ) { //Early exit if Upper limit of Bg integration is reached or number of datapoints are too less for interpolation
        if (verbose == 1) {
            cout<<"Exiting: Either not enough points for interpolation or Nend has exceeded the cap of 1000 or upper limit of integration steps reached without eH = 1."<<endl;
        }
        Plist.assign(npts,1);
        return;
    }
    

    //cout<<"Nend: "<<Nend<<endl;

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
          throw runtime_error("Reheating."); //Signals the end of reheating
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
    //Copy reheating integration to global vectors, +1 avoids duplicates
    N_glob_array.insert(N_glob_array.end(), Nrh.begin(), Nrh.end());
    phi_glob_array.insert(phi_glob_array.end(), phirh.begin(), phirh.end());
    php_glob_array.insert(php_glob_array.end(), phprh.begin(), phprh.end());
    T_glob_array.insert(T_glob_array.end(), Trh.begin(), Trh.end());
    
    //cout<<"Nreh: "<<Nreh<<endl;
     if ( (Nreh>=tend) || ( (Nrh.size()>=maxiter_bg ) && ( w_eff(phirh.back(),phprh.back(),Trh.back())<=0.3)  ) ) { //Early exit if Upper limit of Bg integration is reached.
        if (verbose == 1) {
            cout<<"Exiting: Error in calculation of e-folds at reheating."<<endl;
        }
        Plist.assign(npts,1);
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
        return;
    }

    try {
        //Check if user wants Np to be auto calculated or passed through model_calc.c
        if (want_Np_autocalc == 1) {
            Npp = Calc_Nstar(phiasN,phpasN,TasN,Trh,Nend,Nreh);
        }
        else {
            Npp = Np;
        }
    }
    catch (const exception& e) {
        if (verbose == 1) {
            cout<<"Error in Npivot Calculation."<<endl;
        }
        Plist.assign(npts,1);
        return;
    }

    if (verbose == 1) {
            cout<<"N_end: "<<Nend<<endl;
            cout<<"N_reh: "<<Nreh<<endl;
	    cout<<"N_pivot: "<<Npp<<endl;
    }
    //Set the scale factor normalization
    try {
        a0 = kp/(exp(Npp)*H(phiasN(Npp),phpasN(Npp),TasN(Npp))); //Initialize interpolations after succesful computation of background, Np calc should follow this always.
    }

    catch (const exception& e) {
        if (verbose == 1) {
            cout<<"Error with Npivot parameter"<<endl;
        }
        Plist.assign(npts,1);
        return;
    }

    /* Constrain duration of Inflation. We want atleast 15 e-folds before pivot-scale exit to initialize the larger(smaller) scales(modes) as well. */
    try {
        if ( ((Nend-Npp) < 30.0) || (Npp<15.0)  )   {
            if (verbose == 1) {
                cout<<"Duration of Inflation less than 30.0"<<endl;
            }
            return;
        }
    }
    catch (const exception& e) {
        return;
    }

    //Perturbations
    auto a = [a0] (double N) -> double {
        return a0*exp(N);
    };

    struct Ni_Ne {
        double Ni;
        double Ne;
    };

    auto Calc_Ni_Ne = [Nl,phil,phpl,Tl,a,H] (double k) -> Ni_Ne {
        Ni_Ne nn;
        nn.Ni=0.0;
        nn.Ne = 0.0;
        for (unsigned long i = 0; i < Nl.size(); i++) {
            if ( ( k/(a(Nl[i])*H(phil[i],phpl[i],Tl[i])) ) >= 1000.0  ) {
                nn.Ni = Nl[i];
            }
            else if ( ( k/(a(Nl[i])*H(phil[i],phpl[i],Tl[i])) ) >= 0.1 ) {
                nn.Ne = Nl[i];
            }
            else {
                break;
            }
        }
        return nn;
    };

    //For calculating N and horizon crossing
    auto Calc_Nh = [Nl,phil,phpl,Tl,a,H] (double k) -> double {
        double Ntry=0.0;
        for (unsigned long i = 0; i < Nl.size(); i++) {
            if ( ( k/(a(Nl[i])*H(phil[i],phpl[i],Tl[i])) ) >= 1.0 ) {
                Ntry = Nl[i];
            }
            else {
                break;
            }
        }
        return Ntry;
    };
    
    //Clear out bg integration data
    Nl = vector<double> (); 
    phil = vector<double> (); 
    phpl = vector<double> ();
    Tl = vector<double> (); 

    //SDE Functions

    //Perturbation equations
    auto dpsi_dN = [] (double psi,double dph,double dqr,double Hn,double php) -> double {
        return -psi -( (1.0/2.0) * ( -(php*dph) + (dqr/Hn) ) );
    };

    auto ddqr_dN =[Cr] (double psi,double dph,double dqr,double drr,double Upsn,double Hn,double T,double php) -> double {
        return -(3.0*dqr) -(Upsn*php*dph) -(drr/(3.0*Hn)) - ((4.0*Cr*pow(T,4.0)*psi)/(3.0*Hn));
    };

    auto ddrr_dN = [Cr] (double k,double psi,double dph,double dqr,double drr,double dphp,double Upsn,double Hn,double T,double php,double UpsT,double Upsph,double ai,double psip) -> double {
        double CrT4 = Cr*pow(T,4.0);
        double php2 = php*php;

        double t1 = - (4.0 - ( (UpsT*Hn*(php2)*T)/(4.0*CrT4) ))*drr;
        double t2 = ((k*k)*dqr)/((ai*ai)*Hn);
        double t3 = 2.0*Upsn*Hn*php*dphp;
        double t4 = Upsph*Hn*(php2)*dph;
        double t5 = 4.0 * CrT4*psip;
        double t6 = -(Upsn*Hn*(php2))*psi;

        return t1 + t2 + t3 + t4 + t5 + t6;
    };

    auto ddph_dN = [] (double dphp0) -> double {
        return dphp0;
    };

    auto ddphp_dN = [Cr,Vd,Vdd] (double k,double psi,double dph,double dqr,double drr,double dphp,double Upsn,double Hn,double T,double php,double UpsT,double Upsph,double ai,double psip,double Hpi,double phi) -> double {
        double H2 = Hn*Hn;

        double t1 = -( 3.0 + (Upsn/Hn) + (Hpi/Hn) )*dphp;
        double t2 = -(  ( (k/(ai*Hn))*(k/(ai*Hn)) )  + (Vdd(phi)/(H2)) + (Upsph*php/Hn) )*dph;
        double t3 = -( (UpsT*T*php*drr)/(4.0*Hn*Cr*pow(T,4.0)) );
        double t4 = 4.0*psip*php;
        double t5 = -( (Upsn*php/Hn) + (2.0*Vd(phi)/(H2)) )*psi;

        return t1 + t2 + t3 + t4 + t5;
    };

    //Noise terms
    auto nT = [] (double Upsn,double T,double ai,double Hn) -> double {
        double aiHn = (ai*Hn);
        return sqrt( (2.0*Upsn*T)/(aiHn*aiHn*aiHn) );
    };

    auto nq = [therm] (double Upsn,double T,double ai,double Hn) -> double {

        if (therm == 1) {
            return ( pow(( (9.0*Hn) + (4.0*M_PI*Upsn) ),(1.0/4.0))* sqrt( 1.0/tanh(Hn/(2.0*T)) ) )/ sqrt(M_PI*pow(ai,3.0)*pow(Hn,(3.0/2.0)));
        }
        else if (therm == 0) {
            return ( pow(( (9.0*Hn) + (4.0*M_PI*Upsn) ),(1.0/4.0))* sqrt( 1.0 ) )/ sqrt(M_PI*pow(ai,3)*pow(Hn,(3.0/2.0)));
        }

        else {
            cout<<"Please set therm to either 1 (bose-einstein) or 0 (no thermalization)"<<endl;
            return 0.0;
        }
    };

    auto dW = [] (double h) -> double {
        //Create distribution object from the rng engine
        normal_distribution dist(0.0,sqrt(h));
        return dist(rngen);
    };

    //Euler-Maruyama

    //Struct to store computation products
    struct EM_vars {
    	double psiE;
    	double dqrE;
    	double dphE;
    	double drrE;
    	double dphpE;
    	double Ne;
    };/*
    auto EM = [Calc_Ni_Ne,phiasN,phpasN,TasN,Ups,pT_Ups,pph_Ups,H,Hp,a,dpsi_dN,ddqr_dN,ddph_dN,ddrr_dN,ddphp_dN,dW,nT,nq,rad_noise] (double k,double h) -> EM_vars {
        Ni_Ne nn = Calc_Ni_Ne(k);
        EM_vars em;
        double Ni = nn.Ni;
        double Ne = nn.Ne;

        //Initial Conditions
        double psi0=0.0;
        double dqr0=0.0;
        double drr0=0.0;
        double dph0 = -1.0/(sqrt(2.0*k)*a(Ni));
        double dphp0 = 1.0/(sqrt(2.0*k)*a(Ni));
        double dnT = dW(h);
        double dnq = dW(h);

        //Reused Variables
        double Hi;
        double Upsni;
        double Ti;
        double phpi;
        double UpsTi;
        double Upsphi;
        double ai;
        double psip;
        double Hpi;
        double phi;

        while (Ni<=Ne) {
            //Reused Variables
            phi = phiasN(Ni);
            phpi = phpasN(Ni);
            Ti = TasN(Ni);
            Upsni = Ups(phi,Ti);
            Hi = H(phi,phpi,Ti);
            Hpi = Hp(phi,phpi,Ti);
            UpsTi = pT_Ups(phi,Ti);
            Upsphi = pph_Ups(phi,Ti);
            ai = a(Ni);
            psip = dpsi_dN(psi0,dph0,dqr0,Hi,phpi);

            psi0 = psi0 + (h*dpsi_dN(psi0,dph0,dqr0,Hi,phpi));
            dqr0 = dqr0 + (h*ddqr_dN(psi0,dph0,dqr0,drr0,Upsni,Hi,Ti,phpi));
            dph0 = dph0 + (h*ddph_dN(dphp0));
            
            if (rad_noise == 1) {
                drr0 = drr0 + (h*ddrr_dN(k,psi0,dph0,dqr0,drr0,dphp0,Upsni,Hi,Ti,phpi,UpsTi,Upsphi,ai,psip)) + (dnT*nT(Upsni,Ti,ai,Hi)*(-Hi*Hi*phpi));
            }
            else {
                drr0 = drr0 + (h*ddrr_dN(k,psi0,dph0,dqr0,drr0,dphp0,Upsni,Hi,Ti,phpi,UpsTi,Upsphi,ai,psip));
            }
            
            dphp0 = dphp0 +  (dnT*nT(Upsni,Ti,ai,Hi)) + (h*ddphp_dN(k,psi0,dph0,dqr0,drr0,dphp0,Upsni,Hi,Ti,phpi,UpsTi,Upsphi,ai,psip,Hpi,phi)) + (dnq*nq(Upsni,Ti,ai,Hi));
            Ni+=h;

            dnT = dW(h);
            dnq = dW(h);

        }
        em.psiE = psi0;
        em.dqrE = dqr0;
        em.dphE = dph0;
        em.drrE = drr0;
        em.dphpE = dphp0;
        em.Ne = Ne;
        return em;
    };*/

    auto RI1W1 = [Calc_Ni_Ne,phiasN,phpasN,TasN,Ups,pT_Ups,pph_Ups,H,Hp,a,dpsi_dN,ddqr_dN,ddph_dN,ddrr_dN,ddphp_dN,dW,nT,nq, rad_noise] (double k,double h) -> EM_vars {
        Ni_Ne nn = Calc_Ni_Ne(k);
        EM_vars em;
        double Ni = nn.Ni;
        double Ne = nn.Ne;
        //Reused Variables
        double phi;
        double phpi;
        double Ti;
        double Hi;
        double Upsni;
        double UpsTi;
        double Upsphi;
        double ai;
        double psip;
        double Hpi;
        
        
        //Initial Conditions
        double psi0=0.0;
        double dqr0=0.0;
        double drr0=0.0;
        double dph0 = (-1.0/(sqrt(2.0*k)*a(Ni)));
        double dphp0 = (1.0/(sqrt(2.0*k)*a(Ni)));
        double dnT = dW(h);
        double dnq = dW(h);


        while (Ni<=Ne) {
            //Reused Variables
            phi = phiasN(Ni);
            phpi = phpasN(Ni);
            Ti = TasN(Ni);
            Upsni = Ups(phi,Ti);
            Hi = H(phi,phpi,Ti);
            Hpi = Hp(phi,phpi,Ti);
            UpsTi = pT_Ups(phi,Ti);
            Upsphi = pph_Ups(phi,Ti);
            ai = a(Ni);
            
            double kpsi1 = dpsi_dN(psi0,dph0,dqr0,Hi,phpi);
            double kdqr1 = ddqr_dN(psi0,dph0,dqr0,drr0,Upsni,Hi,Ti,phpi);
            double kdph1 = ddph_dN(dphp0);
            double kdrr1 = ddrr_dN(k,psi0,dph0,dqr0,drr0,dphp0,Upsni,Hi,Ti,phpi,UpsTi,Upsphi,ai,kpsi1);
            double kdphp1 = ddphp_dN(k,psi0,dph0,dqr0,drr0,dphp0,Upsni,Hi,Ti,phpi,UpsTi,Upsphi,ai,kpsi1,Hpi,phi);
            double kdphp1nT = nT(Upsni,Ti,ai,Hi);
            double kdphp1nq = nq(Upsni,Ti,ai,Hi);
            double kdrr1nT = kdphp1nT*(-Hi*Hi*phpi);

            double N2 = Ni+((2.0/3.0)*h);
            double phi2 = phiasN(N2);
            double phpi2 = phpasN(N2);
            double Ti2 = TasN(N2);
            double Upsni2 = Ups(phi2,Ti2);
            double Hi2 = H(phi2,phpi2,Ti2);
            double Hpi2 = Hp(phi2,phpi2,Ti2);
            double UpsTi2 = pT_Ups(phi2,Ti2);
            double Upsphi2 = pph_Ups(phi2,Ti2);
            double ai2 = a(N2);
            double psi2 = psi0+((2.0/3.0)*h*kpsi1);
            double dph2 = dph0+((2.0/3.0)*h*kdph1);
            double dqr2 = dqr0+((2.0/3.0)*h*kdqr1);
            double dphp2 = dphp0+((2.0/3.0)*h*kdphp1) + (dnT*kdphp1nT)  + (dnq*kdphp1nq)  ;
            double drr2 = 0.0;
            if (rad_noise == 1) {
                drr2 = drr0+((2.0/3.0)*h*kdrr1) + ((dnT*kdphp1nT)*(-Hi*Hi*phpi));
            }
            else {
                drr2 = drr0+((2.0/3.0)*h*kdrr1);
            }
            double kpsi2 = dpsi_dN(psi2 , dph2 , dqr2 , Hi2 , phpi2);
            double kdqr2 = ddqr_dN(psi2,dph2,dqr2,drr2,Upsni2,Hi2,Ti2,phpi2);
            double kdph2 = ddph_dN(dphp2);
            double kdrr2 = ddrr_dN(k,psi2,dph2,dqr2,drr2,dphp2,Upsni2,Hi2,Ti2,phpi2,UpsTi2,Upsphi2,ai2,kpsi2);
            double kdphp2 = ddphp_dN(k,psi2,dph2,dqr2,drr2,dphp2,Upsni2,Hi2,Ti2,phpi2,UpsTi2,Upsphi2,ai2,kpsi2,Hpi2,phi2);
            double N2n = Ni+h;
            double phi2n = phiasN(N2n);
            double phpi2n = phpasN(N2n);
            double Ti2n = TasN(N2n);
            double Upsni2n = Ups(phi2n,Ti2n);
            double Hi2n = H(phi2n,phpi2n,Ti2n);
            double ai2n = a(N2n);
            double kdphp2nT = nT(Upsni2n,Ti2n,ai2n,Hi2n);
            double kdphp2nq = nq(Upsni2n,Ti2n,ai2n,Hi2n);
            double kdrr2nT = kdphp2nT*(-Hi2n*Hi2n*phpi2n);

            double N3 = Ni+((2.0/3.0)*h);
            double phi3 = phiasN(N3);
            double phpi3 = phpasN(N3);
            double Ti3 = TasN(N3);
            double Upsni3 = Ups(phi3,Ti3);
            double Hi3 = H(phi3,phpi3,Ti3);
            double Hpi3 = Hp(phi3,phpi3,Ti3);
            double UpsTi3 = pT_Ups(phi3,Ti3);
            double Upsphi3 = pph_Ups(phi3,Ti3);
            double ai3 = a(N3);
            double psi3 = psi0-((1.0/3.0)*h*kpsi1)+(h*kpsi2);
            double dph3 = dph0-((1.0/3.0)*h*kdph1)+(h*kdph2);
            double dqr3 = dqr0-((1.0/3.0)*h*kdqr1)+(h*kdqr2);
            double dphp3 = dphp0-((1.0/3.0)*h*kdphp1)+(h*kdphp2);
            double drr3 = drr0-((1.0/3.0)*h*kdrr1)+(h*kdrr2);
            
            double kpsi3 = dpsi_dN(psi3 , dph3 , dqr3 , Hi3 , phpi3);
            double kdqr3 = ddqr_dN(psi3,dph3,dqr3,drr3,Upsni3,Hi3,Ti3,phpi3);
            double kdph3 = ddph_dN(dphp3);
            double kdrr3 = ddrr_dN(k,psi3,dph3,dqr3,drr3,dphp3,Upsni3,Hi3,Ti3,phpi3,UpsTi3,Upsphi3,ai3,kpsi3);
            double kdphp3 = ddphp_dN(k,psi3,dph3,dqr3,drr3,dphp3,Upsni3,Hi3,Ti3,phpi3,UpsTi3,Upsphi3,ai3,kpsi3,Hpi3,phi3);
                            
            psi0 = psi0 + ((h/4.0)*(kpsi1+(2.0*kpsi2)+kpsi3));
            dqr0 = dqr0 + ((h/4.0)*(kdqr1+(2.0*kdqr2)+kdqr3));
            dph0 = dph0 + ((h/4.0)*(kdph1+(2.0*kdph2)+kdph3));

            if (rad_noise == 1) {
                drr0 = drr0 + ((h/4.0)*(kdrr1+(2.0*kdrr2)+kdrr3)) + ((dnT/2.0)*(kdrr1nT+kdrr2nT));
            }
            else {
                drr0 = drr0 + ((h/4.0)*(kdrr1+(2.0*kdrr2)+kdrr3));
            }

            dphp0 = dphp0 + ((h/4.0)*(kdphp1+(2.0*kdphp2)+kdphp3)) + ((dnT/2.0)*(kdphp1nT+kdphp2nT)) + ((dnq/2.0)*(kdphp1nq+kdphp2nq));
            
            Ni+=h;

            dnT = dW(h);
            dnq = dW(h);

        }
        em.psiE = psi0;
        em.dqrE = dqr0;
        em.dphE = dph0;
        em.drrE = drr0;
        em.dphpE = dphp0;
        em.Ne = Ne;
        return em;
    };

    auto R = [RI1W1,EM_step,Cr,phiasN,phpasN,TasN,H,V] (double k, double h) -> double {
        
        EM_vars sol_r = RI1W1(k,h);
        double Ne = sol_r.Ne;
        double dphr = sol_r.dphE;
        double dqrr = sol_r.dqrE;
        double psir = sol_r.psiE;

        double phiNe = phiasN(Ne);
        double TNe = TasN(Ne);
        double phpNe = phpasN(Ne);
        double HNe = H(phiNe,phpNe,TNe);
        double VNe = V(phiNe);

        double CrT4 = Cr*pow(TNe,4.0);
        double HNephpNe2 = (HNe*phpNe)*(HNe*phpNe);

        double rhotot = CrT4 + ( (HNephpNe2/2.0) + VNe );
        double ptot = ((1.0/3.0)*CrT4) + ( (HNephpNe2/2.0) - VNe );

        double R_sol = ( (HNe/(rhotot+ptot)) * ( -(HNe*phpNe*dphr) + dqrr) ) - psir;

        return R_sol;
    };


    auto P_num = [Nrealz,R] (double k,double h) -> double {
        double sum_R2 = 0.0;
        double R1 = 0.0;
        // Parallelize the loop using OpenMP
        #pragma omp parallel for reduction(+:sum_R2)
        for (int i = 0; i < Nrealz; i++) {
            R1 = R(k,h);
            sum_R2 = sum_R2 + (R1*R1);
        }
        double Rmean = sum_R2/Nrealz;
        return ( (k*k*k)/(2.0*(M_PI*M_PI)) ) * Rmean;
    };

    //cout<<"As: "<<P_num(kp)<<endl;
    

    //Generate power spectrum between 10^klow and 10^kup
    auto Calc_P = [P_num,Q,Calc_Nh,EM_step,phiasN,phpasN,TasN,Calc_Ni_Ne,verbose] (vector <double> &klist, vector <double> &Plist) -> void {
        double h = 0.0;
        double Nh = 0.0;
        vector <double> kl_int;
        //#pragma omp parallel for //Doesn't improve execution time for small npts<=25 atleast
        for (unsigned long i=0;i<klist.size();i++) {
            double ki = klist[i];
            Nh = Calc_Nh(ki);
            Ni_Ne nn = Calc_Ni_Ne(ki);
            //Exit if there aren't enough efolds for the mode to freeze
            if (nn.Ne == Nend) {
                if (verbose == 1) {
                    cout<<"Exiting. Not enough efolds for mode: "<<ki<<" to freeze."<<endl;
                }
                break;
            }

            else if (nn.Ni == 0.0) {
                if (verbose == 1) {
                    cout<<"Exiting. Not enough efolds for mode: "<<ki<<" to initialize."<<endl;
                }
                continue;
            }
            //Adapts the step-size according to the value of Q.
            double Q_val = Q(phiasN(Nh),phpasN(Nh),TasN(Nh));
            
            //Empirical prescription to change step-size with Q to keep accuracy
            if (Q_val<=13600) {
        	h = 4.71550024e-05 - 7.85635682e-08*pow(Q_val,6.29070211e-01) + 9.72292031e-04*pow(9.95511012e-01,Q_val);
            }
            else {
                h = EM_step;
            }
            
            /*if (Q_val<=0.01) {
            	h = EM_step;
            }
            else if (Q_val>0.01 && Q_val<=310) {
                h = EM_step*0.1;
            }
            else {
                h = EM_step*0.01;
            }*/
            //cout<<"N: "<<Nh<<" Q: "<<Q_val<<" h: "<<h<<" k: "<<ki<<endl;
            Plist.push_back(P_num(ki,h));
            kl_int.push_back(ki);
        }
        klist = kl_int;
    };

    try {
        Calc_P(klist,Plist); //Calculate the power spectrum.
    }

    catch (const exception& e) {
        if (verbose == 1) {
            cout<<"Error with power spectrum calculation"<<endl;
        }
        Plist.assign(npts,1);
        return;
    }
    
    auto PT = [phiasN,phpasN,TasN,H] (double NN) -> double {
        double phin = phiasN(NN);
        double phpn = phpasN(NN);
        double Tn = TasN(NN);
        double Hn = H(phin,phpn,Tn);
        return (2.0*pow(Hn,2.0))/pow(M_PI,2.0);
    };

    try {
        PT_kp = PT(Npp); //Calculate the amplitude of tensor power spectrum at kp
    }
    catch (const exception& e) {
        if (verbose == 1) {
            cout<<"Error with tensor power spectrum calculation"<<endl;
        }
        PT_kp = 0.0;
        return;
    }    
}


