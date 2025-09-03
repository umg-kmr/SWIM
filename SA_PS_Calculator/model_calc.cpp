#include "Bg.cpp"

double kp = 0.05;
//int therm : Thermalization of inflaton, 1 = yes = Bose-Einstein; 0 = no

//Klist lower and upper bounds in log10
double klow = -6.0;
double kup = 2.0;
int camb_mcmc = 1; //If you want full power spectrum for CAMB_MCMC then 1 otherwise 0 will give As and ns
int want_Np_autocalc = 1; //Set to 1 if you want the solver to calculate Np itself (only when smooth transition to RD after Inflation) otherwise set 0. In both cases supply Np
int verbose = 1; //Set to one if you want to see the error messages
int GQ_from_file = 1; //If you want to interpolate GQ from a file then 1 otherwise 0 for analytical forms
string GQfname = "../GQ_Calculator/GQ_smooth.dat";

extern "C" {
    int npts = 2000; //No of points to evaluate between klow and kup
    void model (double phi_ini,double gst,double alph, double n,double Cy,double V0,double Np,int c,int p,int therm) {
	double Cr = (M_PI*M_PI / 30.0) * gst;
        double php_ini=0.0;
        double T_ini=0.0;
        double Q_ini=0.0;

        //#### Model Def here ####//
        auto V = [V0,alph,n] (double phi) -> double {
            return V0*exp(-alph*(pow(phi,n)));//0.5*m*m*phi*phi;
        };

        auto Vd = [V0,alph,n] (double phi) -> double {
            return -((n*V0*alph*pow(phi,-1.0 + n))/exp(alph*pow(phi,n)));//m*m*phi;
        };

        auto Ups = [Cy,p,c] (double phi,double T) -> double {
            return Cy * pow(T,p) * pow(phi,c);
        };

        auto set_php_ini = [Vd,V,phi_ini,&php_ini] (double Q_ini) -> void {
            php_ini = -Vd(phi_ini)/(V(phi_ini)*(1.0+Q_ini));
        };

        auto Qi_find = [V,Vd,Cr,p,c,Cy,phi_ini] (double Qi) -> double {
            Qi = pow(10.0,Qi);
            return ( pow((1.0+Qi),(2.0*p))*pow(Qi,(4.0-p)) )  - ( (pow(Cy,4.0)/(9.0*pow(4.0,p)*pow(Cr,p))) * pow(phi_ini,(4.0*c)) * ( pow(Vd(phi_ini),(2.0*p))/pow(V(phi_ini),(p+2.0))  ) ) ;
        };

        try {
            auto logQ = boost::math::tools::bisect(Qi_find, -18.0, 6.0,root_stop(),max_iter);
            Q_ini = (logQ.second + logQ.first)/2;
            Q_ini = pow(10.0,Q_ini);
        }
        catch (const exception& e) {
            if (verbose==1) {
                cout<<"Q_initial couldn't be found"<<endl;
            }
            return;
        }

        set_php_ini(Q_ini);

        auto set_T_ini = [V,phi_ini,php_ini,Cr,&T_ini] (double Q_ini) -> void {
            T_ini = pow( ((Q_ini*V(phi_ini)*php_ini*php_ini)/(4.0*Cr))  , (1.0/4.0) );
        };
        //#################################//

        set_T_ini(Q_ini);
       //Check for T/H>1 where H is slow-roll approximated
        if ( (T_ini*sqrt(3.0)/sqrt(V(phi_ini)) )<=1.0 ) {
            if (verbose==1) {
                cout<<"T<=H, not warm inflation. Exiting."<<endl;
            }
            return;
        }
        //Initialiaze interpolation object for GQ (if requested)
        if (GQ_from_file == 1) {
            try {
                GQ_parser(GQfname); //Readfile
                GQ_interp();//Interpolation object
           } 
           catch (const exception& e) {
                if (verbose==1){
                    cout<<"Unable to initialize GQ file, check if the file is present, fname string is correctly initialized."<<endl;
                }
            } 
        }
        auto GQ = [p] (double Q) -> double {
            if (GQ_from_file == 1) {
                return exp10(logGasQ(log10(Q)));
            }
            else {
                if (p == 3) {
            		return ( ( 1 + (6.12 * pow(Q,2.73) ) ) / pow( (1 + (6.96*pow(Q,0.78))),(0.72) ) ) + ( (0.01*pow(Q,4.61)*(1 + (4.82e-6*pow(Q,3.12)) ) ) / pow( (1 + (6.83e-13*pow(Q,4.12))),2.0 ) );//1.0 + 4.981*pow(Q,1.946) + 0.127*pow(Q,4.330);
                }
                else if (p == 1) {
                    return 1.0 + 0.335*pow(Q,1.364)+ 0.0185*pow(Q,2.315);
                }
                else if (p == -1) {
                    return (1.0 + 0.4*pow(Q,0.77))/pow((1.0 + 0.15*pow(Q,1.09)),2.0);
                }
                else {
                    printf("GQ not found in database. Please modify model_calc.cpp");
                    return 0.0;
                }
            }
        };
        bg_solver (V,Vd,Ups,GQ,Cr,Np,phi_ini,php_ini,T_ini,therm,kp,klow,kup,npts,camb_mcmc,want_Np_autocalc,verbose); //Calculates the power-spectrum

    }

    // Define a function that returns a pointer to the global array
    int get_npts () {
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
    void clear_P (){
        Plist = vector<double> ();
        As_val = 0.0;
        ns_val = 0.0;
        Nend = 0.0;
    }
    void clear_k () {
        klist = vector<double> ();
    }
    //Function to set the global variables
    void set_globals (double kpivot, double kmax, double kmin, int Np_calc, int verbosity,int full_spectrum, int GQ_dat_file) {
        kp = kpivot;
        klow = log10(kmin);
        kup = log10(kmax);
        want_Np_autocalc = Np_calc;
        verbose = verbosity;
        camb_mcmc = full_spectrum;
        GQ_from_file = GQ_dat_file;
    }
}

/*
int main () {
    model(0.295933,5291.1,1.00461,2.41651e+16,2.74485e-37,1.0,0,3,0); //phi_ini, gst, lmbd,Cy, V0, Np, c, p, therm // phi_ini,gst,m,Cy,Np,c,p,therm
    cout<<"Nend: "<<Nend<<endl;
    cout<<"ns: "<<ns_val<<endl;
    cout<<"As: "<<As_val<<endl;
    //for (int i=0;i<=10;i++) {
    //    printf("%e\n",Plist[i]);
    //}
   return 0;
}*/
