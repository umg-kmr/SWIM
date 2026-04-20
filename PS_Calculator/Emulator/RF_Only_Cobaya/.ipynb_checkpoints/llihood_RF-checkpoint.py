import numpy as np
import joblib
from sklearn.ensemble import RandomForestRegressor


data = joblib.load("../RF_Acc_Cobaya/rf_model.pkl")
rf = data["rf"]

#Update: results from SPT+Planck+ACT+DESI https://arxiv.org/pdf/2506.20707 Table 6
#P-ACT results https://arxiv.org/pdf/2503.14452 Table 5
yobs = np.array([3.0586,0.9726]) #,-0.0045])
yerr = np.array([0.0094,0.0028]) #,0.0067])

def logp(phi0,gst,Q0,V0):

    theta = np.array([phi0, gst, Q0, V0], dtype=float)

    theta = np.log10(theta)

    preds = np.array([
        tree.predict(theta.reshape(1, -1))[0]
        for tree in rf.estimators_
    ])

    mean_pred = np.mean(preds, axis=0)
    std_pred  = np.std(preds, axis=0)

    # Use only required observables
    model_fin = mean_pred[:len(yobs)]

    sigma2 = (yerr**2) 

    log_lik = -0.5 * np.sum(
        np.log(2*np.pi*sigma2) +
        ((yobs - model_fin)**2) / sigma2
    )

    return log_lik