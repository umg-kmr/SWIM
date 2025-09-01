#include <iostream> 
#include <cmath> 
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/tools/roots.hpp>
#include <stdexcept>
#include <random>
#include <chrono> // Include for time-based seeding
#include <omp.h>  //For parallel processing
#include <unistd.h> // For getpid()
#include <fstream> 

using namespace std;
using namespace boost::numeric::odeint;


double Nend = 0.0;
double tend = 100.0; //Upper limit of Bg integration
double maxiter_bg = 1e8; //Upper limit on the background integration iterations
uintmax_t max_iter = 1000000; //upper limit for root finding algorithm

//Function to terminate root finding algorithm with some epsilon.
struct root_stop  {
    bool operator() (double r1, double r2)  {
         return abs(r1 - r2) <= 1e-16 + 1e-14*max(abs(r1),abs(r2));
    }
};

//Passing model functions by reference as the background code doesn't modify the original model functions like potential and upsilon.
void bg_solver (const function<double(double)> &V, const function<double(double)>& Vd,const function<double(double,double)>& Ups,double Cr,double phi_ini,double php_ini, double T_ini, double kp, int verbose) {


    typedef boost::array< double , 3 > state_type; //For boost ode solver, number of differential equations (3)
    
    //Arrays to store background calculations
    vector <double> Nl;
    vector <double> phil;
    vector <double> phpl;
    vector <double> Tl;


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
        auto stepper = make_controlled( 1e-10 , 1e-8 , runge_kutta_fehlberg78 < state_type >() );
        state_type x = { phi_ini,php_ini,T_ini }; // initial conditions
        try {
            integrate_adaptive( stepper ,func , x , 0.0 , tend, 1e-6 , obsv ); //make_dense_output<runge_kutta_dopri5< state_type >>( 0.0 , 1.0e-14 )
        }
        catch (const exception& e) {
           // cout<<"Threshold reached, stopping integration!"<<endl;
        }
    };


    SolveODE(); //Call solver to initialize the bg vectors
     //Error handling
    if (Nl.size()<1) {
        if (verbose == 1) {
            cout<<"Solver did not run!"<<endl;
        }
        Nend = 0.0; //Create a P vector of size npts filled with 1s.
        return;
    }

    Nend = Nl.back();

    if ( (Nend>=tend) || (Nl.size()<=4) || ( (Nl.size()>=maxiter_bg ) && (eH(phil.back(),phpl.back(),Tl.back())<0.9999) ) ) { //Early exit if Upper limit of Bg integration is reached or number of datapoints are too less for interpolation
        if (verbose == 1) {
            cout<<"Exiting: Either not enough points for interpolation or Nend has exceeded the cap of 1000 or upper limit of integration steps reached without eH = 1."<<endl;
        }
        return;
    }

}


