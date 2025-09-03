import os
import gc

from cffi import FFI
ffi = FFI()
import numpy as np
from cobaya.theory import Theory
import scipy as sp

## Don't forget to change camb_mcmc to 1 in model_calc.c and recompile otherwise it will not work.
Np_autocalc = int(1)
verbosity = int(1)
want_full_spectrum = int(1) #1 if you want the full power-spectrum and 0 if you want just As and ns to be computed. Keep this to 1 in this case.
read_GQ_from_file = int(1)

ffi.cdef("void model (double phi_ini,double gst,double alph, double n,double Cy,double V0,double Np,int c,int p,int therm);extern int npts; double* get_Plist(); double* get_klist(); void clear_P(); void clear_k();void set_globals (double kpivot, double kmax, double kmin, int Np_calc, int verbosity,int full_spectrum, int GQ_dat_file);",override=True)

def feature_power_spectrum(phi0,gst,alph,n,Cy,V0,Np,c,p,therm,  #model params 
                           kmin=1e-6, kmax=100, # generous, for transfer integrals 
                            kp=0.05):

    lib = ffi.dlopen("./libmodel.so")
    lib.set_globals(kp, kmax, kmin, Np_autocalc, verbosity,want_full_spectrum, read_GQ_from_file) #sets global variables within C++ code
    npts = lib.npts
    c = int(c)
    p = int(p)
    therm = int(therm)
    try:
        lib.model(phi0,gst,alph,n,Cy,V0,Np,c,p,therm)
        ks = np.array(ffi.unpack(lib.get_klist(),npts)) #np.array((ctypes.c_double * 2000).in_dll(lib,"klist"))
        #print(ks)
        Pks = np.array(ffi.unpack(lib.get_Plist(),npts)) #np.array((ctypes.c_double * 2000).in_dll(lib,"Plist"))
        return ks, Pks
    except:
        return [0,0,0,0],[0,0,0,0]
    finally:
        lib.clear_k()
        lib.clear_P()
        ffi.dlclose(lib)
        del lib
        gc.collect()
    


class FeaturePrimordialPk(Theory):
    params = {"phi0":None,"gst":None,"alph":None,"n":None,"Cy":None,"V0":None,"Np":None,"c":None,"p":None,"therm":None}

    def calculate(self, state, want_derived=True, **params_values_dict):
        phi0,gst,alph,n,Cy,V0,Np,c,p,therm = \
            [params_values_dict[p] for p in
             ["phi0","gst","alph","n","Cy","V0","Np","c","p","therm"]]
        ks, Pks = feature_power_spectrum(phi0,gst,alph,n,Cy,V0,Np,c,p,therm)
        state['primordial_scalar_pk'] = {'k': ks, 'Pk': Pks, 'log_regular': False}
    def get_primordial_scalar_pk(self):
        return self.current_state['primordial_scalar_pk']


