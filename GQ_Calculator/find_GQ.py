import numpy as np
from cffi import FFI
ffi = FFI()


#fixed model parameters:
gst = 100.0
V0 = 1.0e-14
lmbd = 1.0
p = int(3)
c = int(0)
therm = int(0)
rad_noise = int(0)
hybrid_inf = int(0)
#ph_crit = M/g

phi_ini,Qini = np.loadtxt("ics.dat",unpack=True,delimiter=",")

ffi.cdef("void model (double phi_ini,double Q_ini, double gst, double V0, double lmbd, int p, int c, int therm, int rad_noise,int hybrid_inf);void clear_vars();void set_phi_crit (double x);extern double Nend;extern double GQ_val;extern double Q_val;",override=True)
lib_pert = ffi.dlopen("./pert/libpert.so")

#For hybrid inflation
if hybrid_inf==1:
    lib_pert.set_phi_crit(ph_crit)

for i,j in zip(phi_ini,Qini):
    lib_pert.model(i,j,gst,V0,lmbd,p,c,therm,rad_noise,hybrid_inf)
    GQ = lib_pert.GQ_val
    Q = lib_pert.Q_val
    lib_pert.clear_vars()
    print(f"Calculating GQ at Q= {Q}")
    with open('GQ.dat', 'a') as fl:
        fl.write(str(Q)+","+str(GQ)+"\n")
