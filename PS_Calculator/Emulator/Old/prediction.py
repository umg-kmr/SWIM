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

# Example input (new data to predict on)
# Make sure this input has the same shape as the training data
X_input = np.array([25.0, 1e-14, 1e6, 129])  # Replace with actual data

# Log-transform the input if necessary (since you applied log10 during training)
X_input = np.log10(X_input)

# Make prediction
predictions = make_predictions(model, X_input, device=device)
print("Predictions: ", predictions)