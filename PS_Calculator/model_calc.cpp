#include "Bg.cpp"

/* Set from the python scripts */
double kp = 0.05;
//int therm : Thermalization of inflaton, 1 = yes = Bose-Einstein; 0 = no

double EM_step = 1e-5; //SDE solver step-size after Q>1e4 
int Nrealz = 2500; //Number of realizations for SDE

//Klist lower and upper bounds in log10
double klow = -6.0;
double kup = 2.0;
int npts = 40; 

int want_Np_autocalc = 1; //Set to 1 if you want the solver to calculate Np itself otherwise set 0. In both cases supply Np
int verbose = 0; //Set to one if you want to see the error messages
int want_FP = 1; //enable the Fokker-Planck (deterministic) approach
/* ######################### */

extern "C" {

    void model (double phi_ini,double gst,double Q_ini,double V0,double Np,int p,int c,int therm,int rad_noise) {

        double Cr =  ((M_PI*M_PI) / 30.0) * gst;
        double php_ini=0.0;
        double T_ini=0.0;
        double Cy=0.0;
        double M = 5.6e-5;
        double alph = sqrt(1.0/8.0);
        double g = 0.47;

        //#### Model Definition here ####//

        auto V = [V0] (double phi) -> double {
            return 0.5*V0*phi*phi;
        };

        auto Vd = [V0] (double phi) -> double {
            return V0*phi;
        };

        auto Vdd = [V0] (double phi) -> double {
            return V0;
        };

        auto set_php_ini = [Vd,V,phi_ini,&php_ini] (double Q_ini) -> void {
            php_ini = -Vd(phi_ini)/(V(phi_ini)*(1.0+Q_ini));
        };

        set_php_ini(Q_ini);

        auto set_T_ini = [V,phi_ini,php_ini,Cr,&T_ini] (double Q_ini) -> void {
            T_ini = pow( ((Q_ini*V(phi_ini)*php_ini*php_ini)/(4.0*Cr))  , (1.0/4.0) );
        };

        set_T_ini(Q_ini);

        /*  Specify your WI dissipation here if not of the form T^p \phi^c  */
        auto Ups_wo_Cy = [p,c,M,g,alph] (double phi,double T) -> double {  //Form of Upsilon without the constant
            double mx = sqrt(((g*g)*(M*M))/2.0 + ((alph*alph)*(T*T)));
            return exp(-mx/T)*(pow(g,4.0)*(M*M)*(T*T)/( pow(mx,3.0) ) ) * (1.0 + (1.0/(sqrt(2*M_PI)))*pow((mx/T),(3.0/2.0)) );
        };

        Cy = 3.0*sqrt(V(phi_ini)/3.0)*Q_ini/Ups_wo_Cy(phi_ini,T_ini);
        
        auto Ups = [Cy,Ups_wo_Cy] (double phi,double T) -> double {
            return Cy * Ups_wo_Cy(phi,T);
        };
        
        /*  Define the partial derivatives of Upsilon if not of the form T^p \phi^c  */
        auto pT_Ups = [p,c,Cy,M,g,alph] (double phi, double T) -> double {
            return (Cy*pow(g,4.0)*M*M*(8.0*pow(T,4.0)*alph*alph*(-2.0*sqrt(M_PI) - (pow(2.0,0.75)*alph*alph)/sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T)) + (pow(2.0,0.25)*pow(g,4.0)*pow(M,4.0)*(sqrt(2.0)*T + 2.0*sqrt(g*g*M*M + 2.0*T*
            T*alph*alph)))/(T*sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T)) + 2.0*g*g*M*M*T*(8.0*sqrt(M_PI)*T - (pow(2.0,0.75)*T*alph*alph)/sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T) + 2.0*pow(2.0,0.25)*T*alph*alph*sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T) + 2.0*sqrt(2.0*M_PI)*sqrt(g*g*M*M + 2.0*T*T*alph*alph))))/(2.0*exp(sqrt((g*g*M*M)/2.0 + T*T*alph*alph)/T)*sqrt(2.0*M_PI)*T*pow((g*g*M*M + 2.0*T*T*alph*alph),2.5));
        };

        auto pph_Ups = [Cy,p,c] (double phi, double T) -> double {
            return 0.0;
        };

        /*
        //Find a value of Q_ini from phi_ini
        //DO NOT USE THIS FUNCTION IF YOUR UPSILON CAN NOT BE WRITTEN IN THE FORM: T^p * \phi^c . SUPPLY A Q_INITIAL INSTEAD
        auto Qi_find = [V,Vd,Cr,p,c,Cy,phi_ini] (double Qi) -> double {
            Qi = pow(10.0,Qi);
            return ( pow((1.0+Qi),(2.0*p))*pow(Qi,(4.0-p)) ) - ( (pow(Cy,4.0)/(9.0*pow(4.0,p)*pow(Cr,p))) * pow(phi_ini,(4.0*c)) * ( pow(Vd(phi_ini),(2.0*p))/pow(V(phi_ini),(p+2.0))  ) ) ;
        };

        try {
            auto res = boost::math::tools::bisect(Qi_find, -20.0, 5.0,root_stop(),max_iter);
            Q_ini = (res.second + res.first)/2;
            Q_ini = pow(10.0,Q_ini);
        }
        catch (const exception& e) {
            if (verbose==1) {
                cout<<"Q_initial couldn't be found"<<endl;
            }
            return;
        }
        */
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        bg_solver (V,Vd,Vdd,Ups,pT_Ups,pph_Ups,Cr,Np,phi_ini,php_ini,T_ini,therm,kp,klow,kup,EM_step,npts,Nrealz,want_Np_autocalc,want_FP,verbose,rad_noise); //Calculates the power-spectrum
    }  //Model Specification ends

    //Function to set the global variables
    void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity, int FP_approach) {
        kp = kpivot;
        EM_step = Em_h;
        Nrealz = N_realizations;
        klow = kmin;
        kup = kmax;
        npts = points_bw_k;
        want_Np_autocalc = Np_calc;
        verbose = verbosity;
        want_FP = FP_approach;
    }

    // Define a function that returns a pointer to the global array
    int get_npts () {
        int npts = 0;
        try {
            npts = Plist.size();
        }
        catch (const exception& e) {
            npts = 0;
        }
        return npts;
    }

    double* get_Plist() {
        double* pt = Plist.data();
        return pt;
    }

    double* get_klist() {
        double* pt = klist.data();
        return pt;
    }
    
    void write_Bg(const char* fname) {
       ofstream fl(fname); //Stores in the format: N,phi,phi',T
       //Output header
       fl<<"N (e-folds), phi (in Mpl), phi' (in Mpl), T (in Mpl)"<<endl;
       for (int i=0;i<N_glob_array.size();i++) {
           fl<<N_glob_array[i]<<","<<phi_glob_array[i]<<","<<php_glob_array[i]<<","<<T_glob_array[i]<<endl;
       }
       fl.close();
    }
    
    void clear_P() {
        Plist = vector<double> ();
        Nend = 0.0;
        rngen.seed(std::chrono::system_clock::now().time_since_epoch().count() ^std::random_device{}() ^omp_get_thread_num() ^getpid()); //Reset global RNG after exit
    }
    void clear_k() { //Also clears bg vectors
        klist = vector<double> ();
        N_glob_array = vector<double> ();
        phi_glob_array = vector<double> ();
        php_glob_array = vector<double> ();
        T_glob_array = vector<double> ();
    }
}

