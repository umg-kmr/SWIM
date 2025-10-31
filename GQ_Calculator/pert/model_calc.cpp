#include "Bg.cpp"

double kp = 0.05;
//int therm : Thermalization of inflaton, 1 = yes = Bose-Einstein; 0 = no

double EM_step = 1e-5; //SDE solver step-size after Q>1e4 (fixed)
int Nrealz = 2048; //Number of realizations for SDE

int want_Np_autocalc = 0; //Set to 1 if you want the solver to calculate Np itself (only when smooth transition to RD after Inflation) otherwise set 0. In both cases supply Np
int verbose = 0; //Set to one if you want to see the error messages

extern "C" {

    void model (double phi_ini,double Q_ini, double gst, double V0, double alph, double n, int p, int c, int therm, int rad_noise, int hybrid_inf) {

        double Cr = (M_PI*M_PI / 30.0) * gst;
        double php_ini=0.0;
        double T_ini=0.0;
        double Cy = 0.0;
        double Np = 0.0;


        //#### Model Def here ####//

        auto V = [V0,alph,n] (double phi) -> double {
            //return V0*log(lmbd*phi)*log(lmbd*phi);
            //return 0.5*V0*phi*phi;
            return V0*exp(-alph*(pow(phi,n)));
            //return (pow(M,4.0)/(4.0*lmbd)) + ((V0/2.0)*phi*phi);
            //return V0 * (1.0 - (gma*phi*phi) );
            //return V0 * (1.0 + cos(phi/f));
            //return V0 * ( F - (4.0*exp(-phi/sqrt(3.0))) + exp(-4.0*phi/sqrt(3.0)) + ( R * ( exp(2.0*phi/sqrt(3.0)) - 1.0 ) ) ) ;
            //return (V0/4.0)*pow(phi,4.0);
        };

        auto Vd = [V0,alph,n] (double phi) -> double {
            //return (2.0*V0*log(lmbd*phi))/phi;
            //return V0*phi;
            return -((n*V0*alph*pow(phi,-1.0 + n))/exp(alph*pow(phi,n)));
            //return V0*phi;
            //return -2.0 * V0 * gma * phi;
            //return (-V0/f) * sin(phi/f);
            //return V0 * (2.0/sqrt(3.0)) * exp(-4.0*phi/sqrt(3.0)) * (-2.0 + (2.0*exp(sqrt(3.0)*phi)) + (R*exp(2.0*sqrt(3.0)*phi)) ) ; 
            //return V0*phi*phi*phi;
        };

        auto Vdd = [V0,alph,n] (double phi) -> double {
            //return (-2.0*V0*(-1.0 + log(lmbd*phi)))/(phi*phi);
            //return V0;
            return -(((-1.0 + n)*n*V0*alph*pow(phi,(-2.0 + n)))/exp(alph*pow(phi,n))) + (n*n*V0*alph*alph*pow(phi,(-2.0 + 2.0*n)))/exp(alph*pow(phi,n));
            //return V0;
            //return -2.0 * V0 * gma;
            //return (-V0/(f*f)) * cos(phi/f);
            //return V0 * (4.0/3.0) * exp(-4.0*phi/sqrt(3.0)) * (4.0 - exp(sqrt(3.0)*phi) + (R*exp(2.0*sqrt(3.0)*phi)) ) ;
            //return 3.0*V0*phi*phi;

        };

        auto set_php_ini = [Vd,V,phi_ini,&php_ini] (double Q_ini) -> void {
            php_ini = -Vd(phi_ini)/(V(phi_ini)*(1.0+Q_ini));
        };

        set_php_ini(Q_ini);

        auto set_T_ini = [V,phi_ini,php_ini,Cr,&T_ini] (double Q_ini) -> void {
            T_ini = pow( ((Q_ini*V(phi_ini)*php_ini*php_ini)/(4.0*Cr))  , (1.0/4.0) );
        };

        set_T_ini(Q_ini);

        auto Ups_wo_Cy = [p,c] (double phi,double T) -> double {  //Form of Upsilon without the constant
            return pow(T,p) * pow(phi,c);
            //double mx = sqrt(((g*g)*(M*M))/2.0 + ((alph*alph)*(T*T)));
            //return exp(-mx/T)*(pow(g,4.0)*(M*M)*(T*T)/( pow(mx,3.0) ) ) * (1.0 + (1.0/(sqrt(2*M_PI)))*pow((mx/T),(3.0/2.0)) );
        };
        
        Cy = 3.0*sqrt(V(phi_ini)/3.0)*Q_ini/Ups_wo_Cy(phi_ini,T_ini);

        auto Ups = [Cy,Ups_wo_Cy] (double phi,double T) -> double {
            return Cy * Ups_wo_Cy(phi,T);
        };

        auto pT_Ups = [Cy,p,c] (double phi, double T) -> double {
            /*return (Cy*pow(g,4.0)*M*M*(8.0*pow(T,4.0)*alph*alph*(-2.0*sqrt(M_PI) - (pow(2.0,0.75)*alph*alph)/sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T)) + (pow(2.0,0.25)*pow(g,4.0)*pow(M,4.0)*(sqrt(2.0)*T + 2.0*sqrt(g*g*M*M + 2.0*T*
            T*alph*alph)))/(T*sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T)) + 2.0*g*g*M*M*T*(8.0*sqrt(M_PI)*T - (pow(2.0,0.75)*T*alph*alph)/sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T) + 2.0*pow(2.0,0.25)*T*alph*alph*sqrt(sqrt(g*g*M*M + 2.0*T*T*alph*alph)/T) + 2.0*sqrt(2.0*M_PI)*sqrt(g*g*M*M + 2.0*T*T*alph*alph))))/(2.0*exp(sqrt((g*g*M*M)/2.0 + T*T*alph*alph)/T)*sqrt(2.0*M_PI)*T*pow((g*g*M*M + 2.0*T*T*alph*alph),2.5));*/
            return p * Cy * pow(T,p-1.0) * pow(phi,c);
        };

        auto pph_Ups = [Cy,p,c] (double phi, double T) -> double {
            //return 0.0;
            return c * Cy * pow(T,p) * pow(phi,c-1.0);
        };


        //#################################//

        bg_solver (V,Vd,Vdd,Ups,pT_Ups,pph_Ups,Cr,Np,phi_ini,php_ini,T_ini,therm,kp,EM_step,Nrealz,want_Np_autocalc,verbose,rad_noise,hybrid_inf); //Calculates the power-spectrum
    }

    void clear_vars () {
        Nend = 0.0;
        GQ_val = 0.0;
        Q_val = 0.0;
    }
    void set_phi_crit (double x) {
        ph_crit = x;
    }
    
    //Function to set the global variables
    void set_globals (int N_realizations, double Nstar, int verbosity) {
        Nrealz = N_realizations;
        Nevol = Nstar;
        verbose = verbosity;
    }

}

