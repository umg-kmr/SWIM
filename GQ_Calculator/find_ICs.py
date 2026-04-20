import numpy as np
from scipy import optimize
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
Nprocs = 24  #Set to the number of CPUs you wish to use for parallelization

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
for i in Qs:
    print(f"Finding initial condition for Q = {i:.4e}")  
    soln = optimize.brute(objfn,ranges=ranges, args=(i,), Ns=1000, full_output=1, disp=False, workers=Nprocs,finish=optimize.fmin)
    phi_in = 10**(soln[0][0])
    print(f"Function value at {i:.4e}: {soln[1]:.4e}")
    # logging.info(f"Function value at {i:.4e}: {soln[1]:.4e}")
    with open("ics.dat", "a") as fl:
        fl.write(str(phi_in)+","+str(i)+"\n")
        
ffi.dlclose(lib)
