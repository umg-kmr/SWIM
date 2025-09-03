import numpy as np
from cffi import FFI
from scipy.optimize import basinhopping, minimize
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

ffi = FFI()

#Power spectrum fitting function with running
def fitting_fn(lnk,lnAs,ns,alphs,betas):
    return lnAs + (lnk-np.log(kp))*(ns-1 + 0.5*alphs*(lnk-np.log(kp)) + (1/6)*betas*((lnk-np.log(kp))**2))

#Observables to aim for
#P-ACT results https://arxiv.org/pdf/2503.14452 Table 5
y = np.array([3.062,0.9752]) #(As,ns)
yerr = np.array([0.011,0.0030])

#set globals
kp = 0.05 #pivot scale in Mpc^-1
em_step = 1e-3 #step-size for SDE solver
Nrealz = int(2400) #number of realizations over which to average, higher number leads to more compute time
kmax = np.log10(kp+0.01) #in log10 -> actual kmax used internally is 10^kmax
kmin = np.log10(kp-0.01) #in log10
points_k = int(10) #number of points to be calculated between the k values specified
Np_autocalc = int(1) # can be set to either 1 (for internal automatic calculation of N_pivot) or 0 (specify an N_pivot value) | in both cases a value for the Np parameter needs to be passed.
verbosity = int(1) #can be set to either 1 or 0, when set to 1 the error messages will be printed if encountered any.

#python script parameters
write_bg = False
fname_bg = "bg.dat"
fname_ps = "ps.dat"

#set model parameters (remember to specify the model in model_calc.cpp file):
p = int(3)
c = int(0)
therm = int(1)
Np = 1.0
rad_noise = int(0)
n = 2.0

ffi.cdef("void model (double phi_ini,double gst,double V0, double alph, double n,double Cy,double Np,int p,int c,int therm,int rad_noise);void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();void write_Bg(const char* fname);",override=True)
lib_pert.set_globals(kp,em_step,Nrealz,kmax,kmin,points_k,Np_autocalc,verbosity)

def logp(x):
    gst = 10**x[0]
    alph = 10**[1]
    V0 = 10**[2]
    Cy = 10**x[3]
    phi0 = (( 1.05*(n-1)/(n*alph) )**(1/n))
    
    lib_pert = ffi.dlopen("../libmodel.so")
    try:
        lib_pert.model(phi0,gst,V0,alph,n,Cy,Np,p,c,therm,rad_noise)
        npts = lib_pert.get_npts()
        Pptr = lib_pert.get_Plist()
        kptr = lib_pert.get_klist()
        Plist = ffi.unpack(Pptr,npts)
        klist = ffi.unpack(kptr,npts)
        
        if np.any(np.isclose(Plist, 1.0)):
            raise ValueError("Power Spectrum = 1")
            
        #Filter raw data, fit to power-spectrum function and obtain observables
        Plist = savgol_filter(Plist,5,3)
        popt,pcov = curve_fit(fitting_fn,np.log(klist),np.log(Plist))
        logAs, ns, alphs, betas = popt
        
        ln10_As = logAs + ( 10.0*np.log(10.0) )
        
        model_fin = np.array([ln10_As,ns])
        sigma2 = np.power(yerr,2)
        log_lik = -0.5* np.sum(np.log(2*np.pi*sigma2) + (((y-model_fin)**2)/(sigma2)))
        nve_log_lik = -log_lik #Maximize log-likelihood
        
        return negtve_log_lik
        
    except:
        return np.inf
    finally:
        lib_pert.clear_k()
        lib_pert.clear_P()
        ffi.dlclose(lib_pert)
        del lib_pert
  
local_minimizer = {
    "method": "L-BFGS-B"}
    
        
def callback(intermediate_result):
    print("Current x:", intermediate_result.x)
    print("Current f(x):", intermediate_result.fun)

x0 = [2.1,0.8,-34.0,14.0]#initial guess in log10 #gst,alph,V0,Cy
soln = basinhopping(func=logp,x0=x0,minimizer_kwargs=local_minimizer,niter=100)
result_params = soln.x
result_likelihood = soln.fun

print("Best-fit parameters: ",result_params)
print("Likelihood_value: ", result_likelihood)
        


