import os
#os.environ["OMP_NUM_THREADS"] = "2"

from cffi import FFI
ffi = FFI()
import numpy as np
import scipy as sp

Np_autocalc = int(1) #Automatically calculate the pivot scale exit.
verbosity = int(0)
want_full_spectrum = int(0) #1 if you want the full power-spectrum and 0 if you want just As and ns to be computed. Keep this to 0 in this case.
read_GQ_from_file = int(1)
kp = 0.05

#these parameters have no effect in this case
kmax = 100.0
kmin = 1e-6

ffi.cdef("void model (double phi_ini,double gst,double alph, double n,double Cy,double V0,double Np,int c,int p,int therm);extern int npts; double* get_Plist(); double* get_klist(); void clear_P(); void clear_k();void set_globals (double kpivot, double kmax, double kmin, int Np_calc, int verbosity,int full_spectrum, int GQ_dat_file);",override=True)

lib = ffi.dlopen("./libmodel.so")
lib.set_globals(kp, kmax, kmin, Np_autocalc, verbosity,want_full_spectrum, read_GQ_from_file) #sets global variables within C++ code

#P-ACT results https://arxiv.org/pdf/2503.14452 Table 5
y = np.array([3.062,0.9752])
yerr = np.array([0.011,0.0030])


def logp(phi0,gst,alph,n,Cy,V0,Np,c,p,therm):
    p = int(p)
    c = int(c)
    therm = int(therm)

    try:
        lib.model(phi0,gst,alph,n,Cy,V0,Np,c,p,therm)
        As = lib.As_val
        ns = lib.ns_val
        if (As == 0.0 or ns == 0.0):
            raise ValueError("Argument must be greater than zero.")

        model_fin = np.array([np.log(1e10*As),ns])
        sigma2 = np.power(yerr,2)
        log_lik = -0.5* np.sum(np.log(2*np.pi*sigma2) + (((y-model_fin)**2)/(sigma2)))
        return log_lik
    except:
        return -np.inf
    finally:
        lib.clear_k()
        lib.clear_P()


