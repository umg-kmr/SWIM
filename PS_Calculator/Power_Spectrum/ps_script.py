import numpy as np
from cffi import FFI
ffi = FFI()

#set globals
kp = 0.05 #pivot scale in Mpc^-1
em_step = 1e-5 #step-size for SDE solver after Q > 10^4
Nrealz = int(2048) #number of realizations over which to average, higher number leads to more compute time
kmax = 2.0 #in log10 -> actual kmax used internally is 10^kmax
kmin = -6.0 #in log10
points_k = int(50) #number of points to be calculated between the k values specified
Np_autocalc = int(1) # can be set to either 1 (for internal automatic calculation of N_pivot) or 0 (specify an N_pivot value) | in both cases a value for the Np parameter needs to be passed.
verbosity = int(1) #can be set to either 1 or 0, when set to 1 the error messages will be printed if encountered any.

#python script parameters
write_bg = True
fname_bg = "bg.dat"
fname_ps = "ps.dat"

#set model parameters (remember to specify the model in model_calc.cpp file):

#fixed parameters:
V0 = 2.05237e-15
gst = 538.664
Q0 = 0.311553
phi0 = 23.6895
Np = 1.0
therm = int(0)
rad_noise = int(0)

#Upsilon parameters:
p = int(3)
c = int(0)


if __name__ == "__main__":
    #Open C++ library
    lib_pert = ffi.dlopen("../libmodel.so")
    #Modify this according to the model_calc.cpp "model" function signature (copy-paste the arguments of "void model" in model_calc.cpp). Nothing else needs to be modified here.
    ffi.cdef("void model (double phi_ini,double gst,double Q_ini,double V0,double Np,int p,int c,int therm,int rad_noise) ; void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();void write_Bg(const char* fname); extern double PT_kp;",override=True)
    
    #Set global parameters:
    lib_pert.set_globals(kp,em_step,Nrealz,kmax,kmin,points_k,Np_autocalc,verbosity)
    
    #pass model paramters | match the function signature
    lib_pert.model(phi0,gst,Q0,V0,Np,p,c,therm,rad_noise)
    
    if (write_bg == True):
        lib_pert.write_Bg(fname_bg.encode('utf-8'))
    
    npts = lib_pert.get_npts()
    Pptr = lib_pert.get_Plist()
    kptr = lib_pert.get_klist()
    Plist = ffi.unpack(Pptr,npts)
    klist = ffi.unpack(kptr,npts)

    PT = lib_pert.PT_kp
    #Write the tensor power spectrum amplitude to calculate r
    with open("PT_kp.dat","w") as f_pt:
        f_pt.write(str(PT))
    
    #write raw power spectrum data to file
    with open(fname_ps, "w") as f_ps:
        f_ps.write("k (in Mpc^-1), P(k)\n") #header
        for i,j in zip(klist,Plist):
            f_ps.write(str(i)+","+str(j)+"\n")
