import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from cffi import FFI
import gc
ffi = FFI()

kp = 0.05

def fitting_fn(lnk,lnAs,ns,alphs,betas):
    return lnAs + (lnk-np.log(kp))*(ns-1 + 0.5*alphs*(lnk-np.log(kp)) + (1/6)*betas*((lnk-np.log(kp))**2))
#fixed model parameters:
p = int(3)
c = int(-2)
therm = int(1)
Np = 1.0
rad_noise = int(0)

ffi.cdef("void model (double phi_ini,double gst,double lmbd,double Cy,double Np,int p,int c,int therm,int rad_noise);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();",override=True)

y = [3.044,0.9649,0.002,0.01]
yerr = [0.014,0.0042,0.010,0.013]
sigma2 = np.power(yerr,2)

def logp(phi0,lmbd,Cy,gst):
    params = [phi0,lmbd,Cy,gst]
    lib_pert = ffi.dlopen("./libpert.so")
    lib_pert.model(phi0,gst,lmbd,Cy,Np,p,c,therm,rad_noise)
    try:
        npts = lib_pert.get_npts()
        if npts==0:
            raise ValueError("Power Spectrum error in calculation.")
        Pptr = lib_pert.get_Plist()
        kptr = lib_pert.get_klist()
        Plist = ffi.unpack(Pptr,npts)
        klist = ffi.unpack(kptr,npts)
        lib_pert.clear_P()
        lib_pert.clear_k()
        Plist = savgol_filter(Plist,11,3)
        popt,pcov = curve_fit(fitting_fn,np.log(klist),np.log(Plist))
        logAs, ns, alphs, betas = popt
        model_cp = [logAs + 10*np.log(10),ns,alphs,betas]
        if (logAs>= 0.0):
            raise ValueError("Power Spectrum greater than or equal to 1")
        log_lik = -0.5* np.sum(np.log(np.multiply(2*np.pi,sigma2)) + ((np.power((np.subtract(y,model_cp)),2))/(sigma2))) 
        #writing data to file
        with open('gen_data.txt', 'a') as f:
            f.write(",".join(str(p) for p in params)+","+str(logAs)+","+str(ns)+","+str(alphs)+","+str(betas)+"\n")
        return log_lik
    except:
        return -np.inf
    # finally:
    #     ffi.dlclose(lib_pert)

    #     del lib_pert
    #     gc.collect()
    #     print("Run correctly without errors")
        
        