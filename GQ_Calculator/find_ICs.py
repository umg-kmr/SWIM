import numpy as np
from scipy import optimize
from cffi import FFI
ffi = FFI()

# logging.basicConfig(level=logging.INFO, format='%(message)s')

#fixed model parameters:
gst = 100.0
V0 = 1.0e-14
lmbd = 1.0
p = int(3)
c = int(0)
hybrid_inf = int(0)
#ph_crit = M/g
#Q range:
Qlow = 1e4
Qup = 1e5
npts = 20 #number of points

dur_N = 60.0 #Duration of inflation parameter
Nprocs = 24

#provide a range of values where you expect the phi_initial to be around. For small-field models one can reduce the upper_bound. The bounds are in log10.
lower_bound = -5.0
upper_bound = 2.0
ranges = ((lower_bound, upper_bound),)

Qs = np.logspace(np.log10(Qlow),np.log10(Qup),npts)

# def pll_comp(i):

ffi.cdef("void model (double phi_ini,double Q_ini, double gst, double V0, double lmbd,int p, int c, int hybrid_inf);void clear_Nend();extern double Nend;void set_phi_crit (double x);",override=True)
lib = ffi.dlopen("./bg/libbg.so")

#For hybrid inflation
if hybrid_inf==1:
    lib.set_phi_crit(ph_crit)

def objfn(x,Q0):
    phi0 = 10**x[0]
    lib.model(phi0, Q0, gst, V0, lmbd, p, c,hybrid_inf)
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

# return (phi_in,i)
        
# with Pool() as pool:  
#     results = pool.map(pll_comp, Qs)
# with open("ics.dat", "w") as f:
#     for phi_in, Qi in results:
#         f.write(str(phi_in)+","+str(Qi)+"\n")

