from scipy.stats import qmc
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from cffi import FFI
ffi = FFI()

ffi.cdef("void model (double phi_ini,double gst,double V0, double alph, double n,double Cy,double Np,int p,int c,int therm,int rad_noise);void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();void write_Bg(const char* fname);",override=True)

# Define bounds for parameter space
#gst, V0, alph, Cy
bounds = [[0.0,3.0],[-50.0,-10.0],[-6.0, 2.0],[0.0, 30.0]]
n_samples = 14#power of 2

d = len(bounds) #number of parameters
# Sobol Sampling
sampler = qmc.Sobol(d, scramble=True)
sample = sampler.random_base2(m=n_samples)

# Scale to parameter bounds
l_bounds = [b[0] for b in bounds]
u_bounds = [b[1] for b in bounds]
scaled_sample = qmc.scale(sample, l_bounds, u_bounds)

print("Number of points generated: ",len(scaled_sample))

def fitting_fn(lnk,lnAs,ns,alphs,betas):
    return lnAs + (lnk-np.log(kp))*(ns-1 + 0.5*alphs*(lnk-np.log(kp)) + (1/6)*betas*((lnk-np.log(kp))**2))
    
#set globals
kp = 0.05 #pivot scale in Mpc^-1
em_step = 1e-3 #step-size for SDE solver
Nrealz = int(1024) #number of realizations over which to average, higher number leads to more compute time
kmax = 2.0 #in log10 -> actual kmax used internally is 10^kmax
kmin = -6.0 #in log10
points_k = int(28) #number of points to be calculated between the k values specified
Np_autocalc = int(1) # can be set to either 1 (for internal automatic calculation of N_pivot) or 0 (specify an N_pivot value) | in both cases a value for the Np parameter needs to be passed.
verbosity = int(0) #can be set to either 1 or 0, when set to 1 the error messages will be printed if encountered any.

#set model parameters (remember to specify the model in model_calc.cpp file):
p = int(3)
c = int(0)
therm = int(1)
Np = 1.0
rad_noise = int(0)
n = 2

cntr = 0

print("Generating observables for the given model...")
for i in scaled_sample:
    #Open C++ library
    lib_pert = ffi.dlopen("../libmodel.so")
    
    #Set global parameters:
    lib_pert.set_globals(kp,em_step,Nrealz,kmax,kmin,points_k,Np_autocalc,verbosity)

    gst, V0, alph, Cy = (np.power(10.0,i)) #Convert from log10
    phi0 = (( 1.05*(n-1)/(n*alph) )**(1/n))
    
    #pass model paramters
    lib_pert.model(phi0,gst,V0,alph,n,Cy,Np,p,c,therm,rad_noise)

    try:

        npts = lib_pert.get_npts()
        Pptr = lib_pert.get_Plist()
        kptr = lib_pert.get_klist()
        Plist = ffi.unpack(Pptr,npts)
        klist = ffi.unpack(kptr,npts)

        if np.any(np.isclose(Plist, 1.0)):
            raise ValueError("Power Spectrum = 1")

        #Filter raw data, fit to power-spectrum function and obtain observables
        Plist = savgol_filter(Plist,11,3)
        popt,pcov = curve_fit(fitting_fn,np.log(klist),np.log(Plist))
        logAs, ns, alphs, betas = popt
    
        if not(logAs>= 0.0):
            cntr=+1
            #writing data to file
            with open('gen_data.txt', 'a') as f:
                f.write(",".join(str(p) for p in i)+","+str(logAs)+","+str(ns)+","+str(alphs)+","+str(betas)+"\n")
    except:
        continue
    finally:
        lib_pert.clear_k()
        lib_pert.clear_P()
        ffi.dlclose(lib_pert)
        del lib_pert
print("Final valid points generated: ",cntr)
