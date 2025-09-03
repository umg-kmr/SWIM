#include "Bg.cpp"

double kp = 0.05;
//int therm : Thermalization of inflaton, 1 = yes = Bose-Einstein; 0 = no

double EM_step = 1e-3; //SDE solver step-size (fixed)
int Nrealz = 2500; //Number of realizations for SDE

//Klist lower and upper bounds in log10
double klow = -6.0;
double kup = 2.0;
int npts = 28; 

int want_Np_autocalc = 1; //Set to 1 if you want the solver to calculate Np itself (only when smooth transition to RD after Inflation) otherwise set 0. In both cases supply Np
int verbose = 0; //Set to one if you want to see the error messages

extern "C" {

    void model (double phi_ini,double gst,double V0, double alph, double n,double Cy,double Np,int p,int c,int therm,int rad_noise) {

        double Cr =  ((M_PI*M_PI) / 30.0) * gst;
        double php_ini=0.0;
        double T_ini=0.0;
        double Q_ini=0.0;

        //#### Model Definition here ####//

        auto V = [V0,alph,n] (double phi) -> double {
            return V0*exp(-alph*(pow(phi,n)));
            //return lmbd*pow(phi,4.0);
        };

        auto Vd = [V0,alph,n] (double phi) -> double {
            return -((n*V0*alph*pow(phi,-1.0 + n))/exp(alph*pow(phi,n)));
            //return 4.0*lmbd*pow(phi,3.0);
        };

        auto Vdd = [V0,alph,n] (double phi) -> double {
            return -(((-1.0 + n)*n*V0*alph*pow(phi,(-2.0 + n)))/exp(alph*pow(phi,n))) + (n*n*V0*alph*alph*pow(phi,(-2.0 + 2.0*n)))/exp(alph*pow(phi,n));
            //return 12.0*lmbd*pow(phi,2.0);
        };

        auto Ups = [Cy,p,c] (double phi,double T) -> double {
            return Cy * pow(T,p) * pow(phi,c);
        };

        auto pT_Ups = [Cy,p,c] (double phi, double T) -> double {
            return p * Cy * pow(T,p-1.0) * pow(phi,c);
        };

        auto pph_Ups = [Cy,p,c] (double phi, double T) -> double {
            return c * Cy * pow(T,p) * pow(phi,c-1.0);
        };

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
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        auto set_php_ini = [Vd,V,phi_ini,&php_ini,Q_ini] () -> void {
            php_ini = -Vd(phi_ini)/(V(phi_ini)*(1.0+Q_ini));
        };
        //#################################//

        set_php_ini();

        auto set_T_ini = [V,phi_ini,Q_ini,php_ini,Cr,&T_ini] () -> void {
            T_ini = pow( ((Q_ini*V(phi_ini)*php_ini*php_ini)/(4.0*Cr))  , (1.0/4.0) );
        };

        set_T_ini();

        bg_solver (V,Vd,Vdd,Ups,pT_Ups,pph_Ups,Cr,Np,phi_ini,php_ini,T_ini,therm,kp,klow,kup,EM_step,npts,Nrealz,want_Np_autocalc,verbose,rad_noise); //Calculates the power-spectrum
    }  //Model Specification ends

    //Function to set the global variables
    void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity) {
        kp = kpivot;
        EM_step = Em_h;
        Nrealz = N_realizations;
        klow = kmin;
        kup = kmax;
        npts = points_bw_k;
        want_Np_autocalc = Np_calc;
        verbose = verbosity;
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

