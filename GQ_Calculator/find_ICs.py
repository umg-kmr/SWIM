import numpy as np
from scipy import optimize
from multiprocessing import Pool
from cffi import FFI
ffi = FFI()

#fixed model parameters (add more depending on the WI model):
gst = 100.0 
V0 = 1e-14

#Upsilon parameters
p = int(3)
c = int(0)

hybrid_inf = int(0)
#ph_crit = M/g  

#Q range:
Qlow = 1e-9
Qup = 1e4
npts = 150 #number of points in Q range

dur_N = 60.0 #Duration of inflation parameter
Nprocs = 4  #Set to the number of CPUs you wish to use for parallelization

#provide a range of values where you expect the phi_initial to be around (for the whole range of Q). For small-field models one can reduce the upper_bound. The bounds are in log10.
lower_bound = np.log10(0.01)
upper_bound = np.log10(30.0)

ranges = ((lower_bound, upper_bound),)

Qs = np.logspace(np.log10(Qlow),np.log10(Qup),npts)

# def pll_comp(i):

#Modify this according to the model_calc.cpp "model" function signature (copy-paste the arguments of "void model" in bg/model_calc.cpp). Nothing else needs to be modified here.
ffi.cdef("void model (double phi_ini,double Q_ini, double gst, double V0,int p, int c, int hybrid_inf) ; void clear_Nend();extern double Nend;void set_phi_crit (double x);",override=True)

lib = ffi.dlopen("./bg/libbg.so") #opens the compiled C++ binary

#For hybrid inflation
if hybrid_inf==1:
    lib.set_phi_crit(ph_crit)

def objfn(x,Q0):
    phi0 = 10**x[0]
    lib.model(phi0, Q0, gst, V0, p, c,hybrid_inf) #Modify to match the signature of the C++ library
    Nend = lib.Nend
    lib.clear_Nend()
    return np.abs(Nend-dur_N)
    
def solve_one_Q(Q0):

    # global search with DIRECT
    soln = optimize.direct(objfn,bounds=[(lower_bound, upper_bound)],args=(Q0,),locally_biased=True,maxfun=350,len_tol=1e-5)
    nfev_total = soln.nfev

    if soln.fun <= 1e-4:
        #if solution is already good enough - skip polishing
        x_best = soln.x[0]
        residual = soln.fun
    else:
        # local polishing around the DIRECT solution
        local = optimize.minimize(objfn,soln.x,args=(Q0,),method="Nelder-Mead",options={"xatol": 1e-6,"fatol": 1e-4,"maxiter": 100})
        nfev_total += local.nfev

        # keep whichever result is better
        if local.fun < soln.fun:
            x_best = local.x[0]
            residual = local.fun
        else:
            x_best = soln.x[0]
            residual = soln.fun

    phi_in = 10**x_best
    return Q0, phi_in, residual, nfev_total

def main():
    with Pool(processes=Nprocs) as pool:
        with open("ics.dat", "a") as fl:
            for Q0, phi_in, residual, nfev in pool.imap(solve_one_Q,Qs,chunksize=1):
                print(f"Q = {Q0:.4e} | "f"phi_initial = {phi_in:.12e} | "f"Function value = {residual:.4e} | "f"Nfev = {nfev}",flush=True)

                fl.write(f"{phi_in},{Q0}\n")
                fl.flush()
        
if __name__ == "__main__":
    main()
