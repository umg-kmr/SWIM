import numpy as np
from sklearn.ensemble import RandomForestRegressor
import random
from cffi import FFI
from scipy.optimize import curve_fit
import joblib

# Global parameters
X_data = []
y_data = []

rf = None
rf_trained = False

load_previous_rf = False #change to true if a previously trained RF model is present

if load_previous_rf:
    try:
        data = joblib.load("rf_model.pkl")
        rf = data["rf"]
        rf_trained = True
        print("Loaded previous RF model")
    except:
        print("No previous RF model found")

update_frequency = 100           # retrain RF every 50 true calls
min_points_before_rf = 100      # wait before first RF fit
max_training_points = 3000      # sliding window limit

rf_uncertainty_tol = 3.0        # trust RF if ratio of tree variance/observational variance
forced_true_fraction = 0.05     # occasionally force true evaluation

rf_calls = 0
solver_calls = 0
new_valid_points = 0
rf_std_values = []
best_logL_seen = -np.inf

#Update: results from SPT+Planck+ACT+DESI https://arxiv.org/pdf/2506.20707 Table 6
#P-ACT results https://arxiv.org/pdf/2503.14452 Table 5
yobs = np.array([3.0586,0.9726]) #,-0.0045])
yerr = np.array([0.0094,0.0028]) #,0.0067])

#Determined by running the solver 50times over the same input parameters for 2048 realizations and 50 points b/w k
yerr_solver = np.array([0.0071,0.0021]) #,0.0009])

scale = yerr + yerr_solver

# =========================================================
# YOUR TRUE SOLVER (PUT YOUR EXPENSIVE CODE HERE)
# =========================================================
ffi = FFI() 

#Power spectrum fitting function with running
def fitting_fn(lnk,lnAs,ns,alphs,betas): 
    return lnAs + (lnk-np.log(kp))*(ns-1 + 0.5*alphs*(lnk-np.log(kp)) + (1/6)*betas*((lnk-np.log(kp))**2))

    
#set globals 
kp = 0.05 #pivot scale in Mpc^-1 
em_step = 1e-5 #step-size for SDE solver 
Nrealz = int(2048) #number of realizations over which to average, higher number leads to more compute time 
kmax = np.log10(1e2) #in log10 -> actual kmax used internally is 10^kmax 
kmin = np.log10(1e-6) #in log10 
points_k = int(50) #number of points to be calculated between the k values specified 
Np_autocalc = int(1) # can be set to either 1 (for internal automatic calculation of N_pivot) or 0 (specify an N_pivot value) | in both cases a value for the Np parameter needs to be passed. 
verbosity = int(0) #can be set to either 1 or 0, when set to 1 the error messages will be printed if encountered any.

#Modify this according to the model_calc.cpp "model" function signature (copy-paste the arguments of "void model" in model_calc.cpp). Nothing else needs to be modified here.
ffi.cdef("void model (double phi_ini,double gst,double Q_ini,double V0,double Np,int p,int c,int therm,int rad_noise);void set_globals (double kpivot, double Em_h, int N_realizations, double kmax, double kmin, int points_bw_k, int Np_calc, int verbosity);int get_npts ();double* get_klist();double* get_Plist();void clear_P();void clear_k();void write_Bg(const char* fname);extern double PT_kp;",override=True)


lib_pert = ffi.dlopen("../../libmodel.so")
lib_pert.set_globals(kp,em_step,Nrealz,kmax,kmin,points_k,Np_autocalc,verbosity) 

lib_pert.clear_k() 
lib_pert.clear_P()

#True Numerical Solver
def true_solver_obs(theta):

    phi0,gst,Q0,V0,Np,p,c,therm,rad_noise = theta
    c = int(c)
    p = int(p)
    therm = int(therm)
    rad_noise = int(rad_noise)

    try:
        lib_pert.model(phi0,gst,Q0,V0,Np,p,c,therm,rad_noise) #match function signature

        npts = lib_pert.get_npts()
        Pptr = lib_pert.get_Plist()
        kptr = lib_pert.get_klist()

        Plist = ffi.unpack(Pptr,npts)
        klist = ffi.unpack(kptr,npts)

        if (np.any(np.isclose(Plist, 1.0))) or (len(Plist) < 5):
            raise ValueError("Power Spectrum = 1")

        popt,pcov = curve_fit(fitting_fn,np.log(klist),np.log(Plist))
        logAs, ns, alphs, betas = popt

        PT = lib_pert.PT_kp
        r = PT/np.exp(logAs)
        if (r>0.038): #95%, P-ACT-LB-BK18
            raise ValueError("Tensor-to-scalar ratio out of bounds.")

        ln10_As = logAs + (10.0*np.log(10.0))

        return np.array([ln10_As, ns, alphs, betas])

    except:
        return None
    finally:
        lib_pert.clear_k()
        lib_pert.clear_P()


#RF Trainer
def train_rf():
    global rf, rf_trained

    X_all = np.array(X_data)
    y_all = np.array(y_data)

    if len(X_all) < min_points_before_rf:
        return

    # Sliding window
    if len(X_all) > max_training_points:
        X_train = X_all[-max_training_points:]
        y_train = y_all[-max_training_points:]
    else:
        X_train = X_all
        y_train = y_all

    rf = RandomForestRegressor(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=1,
        n_jobs=-1
    )

    rf.fit(X_train, y_train)

    rf_trained = True
    print("RF trained on", len(X_train), "points")



#Likelihood function

#The function signature of logp should match with "void model" function of C++ library
def logp(phi0,gst,Q0,V0,Np,p,c,therm,rad_noise):

    global rf_trained, rf_std_values
    global rf_calls, solver_calls, new_valid_points
    global best_logL_seen

    #Modify according to your model parameters
    theta_true = np.array([phi0,gst,Q0,V0,Np,p,c,therm,rad_noise],dtype=float)
    theta_rf = np.array([phi0, gst, Q0, V0], dtype=float)
    theta_rf = np.log10(theta_rf)

    # Try RF
    if rf_trained:
        

        preds = np.array([
        tree.predict(theta_rf.reshape(1, -1))[0]
        for tree in rf.estimators_
        ])

        mean_pred = np.mean(preds, axis=0)
        std_pred  = np.std(preds, axis=0)

        # normalize uncertainty
        rel_std = std_pred[:len(yobs)] / scale
        score = np.mean(rel_std)

        rf_std_values.append(score)
    
        model_fin = mean_pred[:len(yobs)]
    
        # RF trust gate
        if score < rf_uncertainty_tol:
            if random.random() > forced_true_fraction:
    
                sigma2 = (yerr**2) + (yerr_solver**2)
    
                log_lik = -0.5 * np.sum(np.log(2*np.pi*sigma2) +((yobs - model_fin)**2)/sigma2)
    
                if log_lik > best_logL_seen:
                    pass
                else:
                    rf_calls += 1
                    return log_lik

    # True solver fallback
    solver_calls += 1
    obs = true_solver_obs(theta_true)

     # Print usage statistics every 100 calls
    total_calls = rf_calls + solver_calls
    if total_calls > 0 and total_calls % 100 == 0:

        print("\n------ Adaptive Stats ------")
        print("Total calls:", total_calls)
        print("Valid Points Stored for RF:", len(X_data))
        print("RF usage:", round(100*rf_calls/total_calls,2), "%")
        print("Solver usage:", round(100*solver_calls/total_calls,2), "%")
    
        if len(rf_std_values) > 0:
            print("Mean RF std:", round(np.mean(rf_std_values),4))
            print("Median RF std:", round(np.median(rf_std_values),4))
            print("Max RF std:", round(np.max(rf_std_values),4))
    
        print("----------------------------\n")

    if obs is None:
        return -np.inf

    # Store training data
    X_data.append(theta_rf)
    y_data.append(obs)
    new_valid_points += 1

    if len(X_data) > 2*max_training_points:
        X_data.pop(0)
        y_data.pop(0)

    # Retrain if needed
    if len(X_data) >= min_points_before_rf and new_valid_points >= update_frequency:
        train_rf()
        new_valid_points = 0
        
        # Save RF model occasionally
        joblib.dump({"rf": rf,}, "rf_model.pkl")
        
    # Use only first two observables for likelihood
    model_fin = obs[:len(yobs)]

    sigma2 = (yerr**2) + (yerr_solver**2)
    log_lik = -0.5 * np.sum(np.log(2*np.pi*sigma2)+ ((yobs - model_fin)**2)/sigma2)

    if log_lik > best_logL_seen:
        best_logL_seen = log_lik

    return log_lik
