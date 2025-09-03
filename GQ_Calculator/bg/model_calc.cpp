#include "Bg.cpp"

double kp = 0.05;

int verbose = 0; //Set to one if you want to see the error messages

extern "C" {

    void model (double phi_ini,double Q_ini, double gst, double V0, double alph, double n,int p, int c, int hybrid_inf) {

        double Cr = (M_PI*M_PI / 30.0) * gst;
        double php_ini=0.0;
        double T_ini=0.0;
        double Cy = 0.0;


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
            //return 0.5*(m*m)*(phi*phi);
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
            //return m*m*phi;
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


        //#################################//


        bg_solver (V,Vd,Ups,Cr,phi_ini,php_ini,T_ini,kp,verbose,hybrid_inf); //Calculates the power-spectrum
    }

    void clear_Nend () {
        Nend = 0.0;
    }
    
    void set_phi_crit (double x) {
        ph_crit = x;
    }


}

