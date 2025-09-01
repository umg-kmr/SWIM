import numpy as np
from ps_script import lmbd, gst, p, c, Cy

#Define model here again
def V(phi):
    return lmbd*(phi**4.0)
def Vd(phi):
    return 4.0*lmbd*(phi**3.0)
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
