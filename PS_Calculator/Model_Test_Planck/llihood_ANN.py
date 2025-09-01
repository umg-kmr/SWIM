import torch.nn as nn
import torch
from torch.utils.data import TensorDataset, DataLoader
import numpy as np
from neuralnet_arch import Net #Import the neural net architecture
from cobaya.theory import Theory

NUM_INPUTS = 4
NUM_OUTPUTS = 4
MODEL_PATH = "test.pth"

# -----------------------------
# Load the model
# -----------------------------
def load_model(checkpoint_path, device='cpu'):
    # Initialize the model
    model = Net(input_dim=NUM_INPUTS, output_dim=NUM_OUTPUTS)
    model.load_state_dict(torch.load(checkpoint_path, map_location=device))
    model.to(device)  # Ensure model is on the right device (CPU or GPU)
    model.eval()  # Set model to evaluation mode
    return model

# -----------------------------
# Make Predictions
# -----------------------------
def make_predictions(model, X_input, device='cpu'):
    # Convert the input to a tensor and move to the correct device
    X_input_tensor = torch.tensor(X_input, dtype=torch.float32).to(device)
    
    # Make the prediction (no gradients needed for inference)
    with torch.no_grad():
        predictions = model(X_input_tensor)

    return predictions.cpu().numpy()  # Convert back to NumPy array


#define power-spectrum function
def Pr(lnAs, ns, alphas, betas, kp, k):
    P = np.exp(lnAs)*(k/kp)**(ns-1 + (0.5*alphas*np.log(k/kp)) + (betas/6)*(np.log(k/kp)**2) )
    return P

# Load model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = load_model(MODEL_PATH, device=device)
        
def feature_power_spectrum(phi0,lmbd,Cy,gst,  #model params 
                           kmin=1e-6, kmax=100, # generous, for transfer integrals ignored in this case
                            kp=0.05):
    X_input = np.array([phi0,lmbd,Cy,gst])
    X_input = np.log10(X_input)
    predictions = make_predictions(model, X_input, device=device)
    logAs, ns, alphs, betas = predictions
    ks = np.array(ffi.unpack(lib.get_klist(),npts)) #np.array((ctypes.c_double * 2000).in_dll(lib,"klist"))
    #print(ks)
    Pks = np.array([Pr(logAs,ns,alphs,betas,kp,i) for i in ks])
    #print(Pks[1000])
    return ks, Pks


class FeaturePrimordialPk(Theory):
    params = {"phi0":None,"lmbd":None,"Cy":None,"gst":None}
    kp = 0.05

    def calculate(self, state, want_derived=True, **params_values_dict):
        phi0,lmbd,Cy,gst = \
            [params_values_dict[p] for p in
             ["phi0","lmbd","Cy","gst"]]
        ks, Pks = feature_power_spectrum(
            phi0,lmbd,Cy,gst,kmin=1e-6, kmax=100,
            kp=self.kp)
        state['primordial_scalar_pk'] = {'k': ks, 'Pk': Pks, 'log_regular': False}
    def get_primordial_scalar_pk(self):
        return self.current_state['primordial_scalar_pk']


