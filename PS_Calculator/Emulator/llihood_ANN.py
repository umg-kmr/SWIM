import torch.nn as nn
import torch
from torch.utils.data import TensorDataset, DataLoader
import numpy as np
from neuralnet_arch import Net #Import the neural net architecture

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


# Load model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = load_model(MODEL_PATH, device=device)

y = [3.044,0.9649,0.002,0.01]
yerr = [0.014,0.0042,0.010,0.013]
sigma2 = np.power(yerr,2)

def logp(phi0,lmbd,Cy,gst):
    X_input = np.array([phi0,lmbd,Cy,gst])
    X_input = np.log10(X_input)
    try:
        predictions = make_predictions(model, X_input, device=device)
        logAs, ns, alphs, betas = predictions
        model_cp = [logAs + 10*np.log(10),ns,alphs,betas]
        if (logAs>= 0.0):
            raise ValueError("Power Spectrum greater than or equal to 1")
        log_lik = -0.5* np.sum(np.log(np.multiply(2*np.pi,sigma2)) + ((np.power((np.subtract(y,model_cp)),2))/(sigma2))) 
        return log_lik
    except:
        return -np.inf
