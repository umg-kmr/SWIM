from cffi import FFI
ffi = FFI()
import numpy as np
import scipy as sp

Np_autocalc = int(1) #Enable automatic calculation of the pivot scale exit with 1 and disable with 0. If set to zero, set the correct value of pivot scale using Np parameter. Note that the integration starts from N=0.
verbosity = int(0)
want_full_spectrum = int(0) #1 if you want the full power-spectrum and 0 if you want just As and ns to be computed. Keep this to 0 in this case.
read_GQ_from_file = int(1) #Read pre-computed G(Q) from data file when 1 otherwise internal analytical approximations used with 0. Recommended to use GQ_Calculator to pre-compute G(Q)
kp = 0.05 #pivot scale in Mpc^-1

#these parameters have no effect in this case but need to be supplied
kmax = 100.0
kmin = 1e-6

#Modify this according to the model_calc.cpp "model" function signature (copy-paste the arguments of "void model" in model_calc.cpp). Nothing else needs to be modified here.
ffi.cdef("void model (double phi_ini,double gst,double Q_ini,double V0,double Np,int c,int p,int therm) ; extern int npts; extern double As_val; extern double ns_val; double* get_Plist(); double* get_klist(); void clear_P(); void clear_k();void set_globals (double kpivot, double kmax, double kmin, int Np_calc, int verbosity,int full_spectrum, int GQ_dat_file);",override=True)

lib = ffi.dlopen("./libmodel.so")
lib.set_globals(kp, kmax, kmin, Np_autocalc, verbosity,want_full_spectrum, read_GQ_from_file) #sets global variables within C++ code

#Update: results from SPT+Planck+ACT+DESI https://arxiv.org/pdf/2506.20707 Table 6
#P-ACT results https://arxiv.org/pdf/2503.14452 Table 5
y = np.array([3.0586,0.9726])
yerr = np.array([0.0094,0.0028])


#The function signature of logp should match with "void model" function of C++ library
def logp(phi0,gst,Q0,V0,Np,c,p,therm):
    p = int(p)
    c = int(c)
    therm = int(therm)

    try:
        lib.model(phi0,gst,Q0,V0,Np,c,p,therm) #match the function signature
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


