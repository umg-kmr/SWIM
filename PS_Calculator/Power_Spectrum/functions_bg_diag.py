import numpy as np
from ps_script import V0, alph, n, gst, p, c, Cy

#Define model here again
def V(phi):
    return V0*np.exp(-alph*(phi**n))
def Vd(phi):
    return -((n*V0*alph*(phi**(-1.0 + n)))/np.exp(alph*(phi**n)))
def Ups(phi,T):
    return Cy * (T**p) * (phi**c)
Cr = ( (np.pi**2) / 30.0 )*gst

def H(phi,phip,T):
    return np.sqrt(2 * ( V(phi) + (Cr*(T**4))) * ( ( 6 - (phip**2) )**(-1) ) )
def Hp(phi,phip,T):
    return - ( (H(phi,phip,T)*(phip**2))/2 ) - ( (2/3) * ((Cr*(T**4))/H(phi,phip,T)) )
def eH(phi,phip,T):
    return -Hp(phi,phip,T)/H(phi,phip,T)
def Q(phi,phip,T):
    return Ups(phi,T)/(3*H(phi,phip,T))
