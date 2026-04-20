import numpy as np
from ps_script import V0, gst, p, c, Q0, phi0

Cr = ( (np.pi**2) / 30.0 )*gst

#Define model here again
def V(phi):
    return (V0/4.0)*(phi**4)
def Vd(phi):
    return V0*(phi**3)

def Ups_wo_Cy(phi,T):
    return (T**p) * (phi**c)
    
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
