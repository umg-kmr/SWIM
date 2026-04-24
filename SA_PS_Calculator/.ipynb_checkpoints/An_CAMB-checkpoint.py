from cffi import FFI
ffi = FFI()
import numpy as np
from cobaya.theory import Theory
import scipy as sp

Np_autocalc = int(1) #Enable automatic calculation of the pivot scale exit with 1 and disable with 0. If set to zero, set the correct value of pivot scale using Np parameter. Note that the integration starts from N=0.
verbosity = int(0)
want_full_spectrum = int(1) #1 if you want the full power-spectrum and 0 if you want just As and ns to be computed. Keep this to 1 in this case.
read_GQ_from_file = int(1) #Read pre-computed G(Q) from data file when 1 otherwise internal analytical approximations used with 0. Recommended to use GQ_Calculator to pre-compute G(Q)

#Modify this according to the model_calc.cpp "model" function signature (copy-paste the arguments of "void model" in model_calc.cpp). Nothing else needs to be modified here.
ffi.cdef("void model (double phi_ini,double gst,double Q_ini,double V0,double Np,int c,int p,int therm) ; extern int npts; double* get_Plist(); double* get_klist(); void clear_P(); void clear_k();void set_globals (double kpivot, double kmax, double kmin, int Np_calc, int verbosity,int full_spectrum, int GQ_dat_file);",override=True)

lib = ffi.dlopen("./libmodel.so")

#The function signature of feature_power_spectrum should match with "void model" function of C++ library
def feature_power_spectrum(phi0,gst,Q0,V0,Np,c,p,therm,  #model params 
                           kmin=1e-6, kmax=100, # generous, for transfer integrals 
                            kp=0.05):
    
    lib.set_globals(kp, kmax, kmin, Np_autocalc, verbosity,want_full_spectrum, read_GQ_from_file) #sets global variables within C++ code
    npts = lib.npts
    c = int(c)
    p = int(p)
    therm = int(therm)
    try:
        lib.model(phi0,gst,Q0,V0,Np,c,p,therm) #match the function signature
        ks = np.array(ffi.unpack(lib.get_klist(),npts)) 
        Pks = np.array(ffi.unpack(lib.get_Plist(),npts)) 
        return ks, Pks
    except:
        return [0,0,0,0],[0,0,0,0]
    finally:
        lib.clear_k()
        lib.clear_P()

#Following Cobaya guide, change the following section to include your WI model parameters
class FeaturePrimordialPk(Theory):
    #define parameter names here
    params = {"phi0":None,"gst":None,"Q0":None,"V0":None,"Np":None,"c":None,"p":None,"therm":None}
    kp = 0.05
    def calculate(self, state, want_derived=True, **params_values_dict):
        #list the parameters here
        phi0,gst,Q0,V0,Np,c,p,therm = \
        [params_values_dict[itr] for itr in ["phi0","gst","Q0","V0","Np","c","p","therm"]] 
        ks, Pks = feature_power_spectrum(phi0,gst,Q0,V0,Np,c,p,therm,kp=self.kp) #match the function signature
        state['primordial_scalar_pk'] = {'k': ks, 'Pk': Pks, 'log_regular': False}
    def get_primordial_scalar_pk(self):
        return self.current_state['primordial_scalar_pk']


