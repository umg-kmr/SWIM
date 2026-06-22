from cffi import FFI
ffi = FFI()
import numpy as np
from cobaya.theory import Theory
from scipy.optimize import curve_fit

kp = 0.05
#Power spectrum fitting function with running
def fitting_fn(lnk,lnAs,ns,alphs,betas): 
    return lnAs + (lnk-np.log(kp))*(ns-1 + 0.5*alphs*(lnk-np.log(kp)) + (1/6)*betas*((lnk-np.log(kp))**2))

#unused parameters. Have no effect on the solver but need to be passed to C++ library
em_step = 1e-5 
Nrealz = int(2048) 

#edit these if required
points_k = int(1000) #number of points to be calculated between the k values specified
Np_autocalc = int(1) #Enable automatic calculation of the pivot scale exit with 1 and disable with 0. If set to zero, set the correct value of pivot scale using Np parameter. Note that the integration starts from N=0.
verbosity = int(0)
want_FP = int(1) #Do not change for this case

#Modify this according to the model_calc.cpp "model" function signature (copy-paste the arguments of "void model" in model_calc.cpp). Nothing else needs to be modified here.
ffi.cdef("void model (double phi_ini,double gst,double Q_ini,double V0,double Np,int p,int c,int therm,int rad_noise) ; void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity,int FP_approach);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();void write_Bg(const char* fname); extern double PT_kp;",override=True)

lib_pert = ffi.dlopen("../libmodel.so")

#The function signature of feature_power_spectrum should match with "void model" function of C++ library
def feature_power_spectrum(phi0,gst,Q0,V0,Np,p,c,therm,rad_noise,  #model params 
                           kmin=1e-6, kmax=100, # generous, for transfer integrals 
                            kp=0.05):
    
    lib_pert.set_globals(kp,em_step,Nrealz,np.log10(kmax),np.log10(kmin),points_k,Np_autocalc,verbosity,want_FP) #sets global variables within C++ code
    c = int(c)
    p = int(p)
    therm = int(therm)
    rad_noise = int(rad_noise)
    try:
        lib_pert.model(phi0,gst,Q0,V0,Np,p,c,therm,rad_noise) #match the function signature
        
        npts = lib_pert.get_npts()
        Pptr = lib_pert.get_Plist()
        kptr = lib_pert.get_klist()
        Pks = ffi.unpack(Pptr,npts)
        ks = ffi.unpack(kptr,npts)
        if (np.any(np.isclose(Pks, 1.0))) or (len(Pks) < 5):
            raise ValueError("Power Spectrum = 1")
        popt,pcov = curve_fit(fitting_fn,np.log(ks),np.log(Pks))
        logAs, ns, alphs, betas = popt
            
        PT = lib_pert.PT_kp
        r = PT/np.exp(logAs)
        if (r>0.038):
            raise ValueError("r greater than upper bound!")
        return ks, Pks
    except:
        return [0,0,0,0],[0,0,0,0]
    finally:
        lib_pert.clear_k()
        lib_pert.clear_P()

#Following Cobaya guide, change the following section to include your WI model parameters
class FeaturePrimordialPk(Theory):
    #define parameter names here
    params = {"phi0":None,"gst":None,"Q0":None,"V0":None,"Np":None,"p":None,"c":None,"therm":None,"rad_noise":None}
    kp = 0.05
    def calculate(self, state, want_derived=True, **params_values_dict):
        #list the parameters here
        phi0,gst,Q0,V0,Np,p,c,therm,rad_noise = \
        [params_values_dict[itr] for itr in ["phi0","gst","Q0","V0","Np","p","c","therm","rad_noise"]] 
        ks, Pks = feature_power_spectrum(phi0,gst,Q0,V0,Np,p,c,therm,rad_noise,kp=self.kp) #match the function signature
        state['primordial_scalar_pk'] = {'k': ks, 'Pk': Pks, 'log_regular': False}
    def get_primordial_scalar_pk(self):
        return self.current_state['primordial_scalar_pk']


