import os
import gc
#os.environ["OMP_NUM_THREADS"] = "2"

from cffi import FFI
ffi = FFI()
import numpy as np
import scipy as sp

ffi.cdef("void model (double phi_ini,double gst,double lmbd,double Cy,double V0,double Np,int c,int p,int therm);extern double As_val; extern double ns_val; void clear_k();void clear_P();")
#lib = ffi.dlopen("./libmodel.so")

#P-ACT results https://arxiv.org/pdf/2503.14452 Table 5
y = np.array([3.062,0.9752])
yerr = np.array([0.011,0.0030])


def logp(phi0,gst,lmbd,Cy,V0,Np,c,p,therm):
    p = int(p)
    c = int(c)
    therm = int(therm)
    lib = ffi.dlopen("./libmodel.so")
    try:
        lib.model(phi0,gst,lmbd,Cy,V0,Np,c,p,therm)
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
        ffi.dlclose(lib)
        del lib
        gc.collect()


