import numpy as np
from cffi import FFI
ffi = FFI()

#set globals
kp = 0.05 #pivot scale in Mpc^-1
em_step = 1e-3 #step-size for SDE solver
Nrealz = int(2400) #number of realizations over which to average, higher number leads to more compute time
kmax = 2.0 #in log10 -> actual kmax used internally is 10^kmax
kmin = -6.0 #in log10
points_k = int(48) #number of points to be calculated between the k values specified
Np_autocalc = int(1) # can be set to either 1 (for internal automatic calculation of N_pivot) or 0 (specify an N_pivot value) | in both cases a value for the Np parameter needs to be passed.
verbosity = int(1) #can be set to either 1 or 0, when set to 1 the error messages will be printed if encountered any.

#python script parameters
write_bg = True
fname_bg = "bg.dat"
fname_ps = "ps.dat"

#set model parameters (remember to specify the model in model_calc.cpp file):
p = int(3)
c = int(0)
therm = int(1)
Np = 1.0
rad_noise = int(0)
n = 2.0
alph = 9.22
phi0 = (( 1.05*(n-1)/(n*alph) )**(1/n))
gst = 129
V0 = 10**-14.167655244919862
Cy = 10**5.542267426870692

if __name__ == "__main__":
    #Open C++ library
    lib_pert = ffi.dlopen("../libmodel.so")
    ffi.cdef("void model (double phi_ini,double gst,double V0, double alph, double n,double Cy,double Np,int p,int c,int therm,int rad_noise);void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();void write_Bg(const char* fname);",override=True)
    
    #Set global parameters:
    lib_pert.set_globals(kp,em_step,Nrealz,kmax,kmin,points_k,Np_autocalc,verbosity)
    
    #pass model paramters
    lib_pert.model(phi0,gst,V0,alph,n,Cy,Np,p,c,therm,rad_noise)
    
    if (write_bg == True):
        lib_pert.write_Bg(fname_bg.encode('utf-8'))
    
    npts = lib_pert.get_npts()
    Pptr = lib_pert.get_Plist()
    kptr = lib_pert.get_klist()
    Plist = ffi.unpack(Pptr,npts)
    klist = ffi.unpack(kptr,npts)
    
    #write raw power spectrum data to file
    with open(fname_ps, "w") as f_ps:
        f_ps.write("k (in Mpc^-1), P(k)\n") #header
        for i,j in zip(klist,Plist):
            f_ps.write(str(i)+","+str(j)+"\n")
