import numpy as np
from ps_script import V0, gst, p, c, Q0, phi0

Cr = ( (np.pi**2) / 30.0 )*gst

#Define model here again
def V(phi):
    return 0.5*V0*phi*phi
def Vd(phi):
    return V0*phi

def Ups_wo_Cy(phi,T):
    M = 5.6e-5;
    alph = np.sqrt(1.0/8.0);
    g = 0.47;
    mx = np.sqrt(((g*g)*(M*M))/2.0 + ((alph*alph)*(T*T)));
    return np.exp(-mx/T)*((g**4.0)*(M*M)*(T*T)/( (mx**3.0) ) ) * (1.0 + (1.0/(np.sqrt(2.0*np.pi)))*((mx/T)**(3.0/2.0)) )
    
php0 = -Vd(phi0)/(V(phi0)*(1.0+Q0))
T0 = ( ((Q0*V(phi0)*php0*php0)/(4.0*Cr))**(1.0/4.0) );
Cy = 3.0*np.sqrt(V(phi0)/3.0)*Q0/Ups_wo_Cy(phi0,T0)

def Ups(phi,T):
    return Cy * Ups_wo_Cy(phi,T)


def H(phi,phip,T):
    return np.sqrt(2 * ( V(phi) + (Cr*(T**4))) * ( ( 6 - (phip**2) )**(-1) ) )
def Hp(phi,phip,T):
    return - ( (H(phi,phip,T)*(phip**2))/2 ) - ( (2/3) * ((Cr*(T**4))/H(phi,phip,T)) )
def eH(phi,phip,T):
    return -Hp(phi,phip,T)/H(phi,phip,T)
def Q(phi,phip,T):
    return Ups(phi,T)/(3*H(phi,phip,T))
def weff(phi,phip,T):
    Hi = H(phi,phip,T);
    Vi = V(phi);
    p_tot = ((1.0/3.0)*Cr*(T**4)) + (Hi*phip*Hi*phip*0.5) - Vi;
    rho_tot = (Cr*(T**4)) + (Hi*phip*Hi*phip*0.5) + Vi; 
    return p_tot/rho_tot;

