import torch.nn as nn
import torch
from torch.utils.data import TensorDataset, DataLoader
import numpy as np

from neuralnet_arch import Net #Import the neural net architecture

# -----------------------------
# Global Parameters
# -----------------------------
NUM_INPUTS = 4 #Change according to the number of input (model) parameters
NUM_OUTPUTS = 4 #No need to change usually (predicts: As, ns, alphas, betas)
BATCH_SIZE = 64
NUM_EPOCHS = 1000
LEARNING_RATE = 1e-3
VAL_SPLIT = 0.2 #Default 20% split
RANDOM_SEED = 42
DATA_PATH = "gen_data.txt"  # <-- change this to your data file path
MODEL_PATH = "test.pth"

# -----------------------------
# Load and Split Data
# -----------------------------
def load_and_split_data(file_path, val_split=0.2, random_seed=42):
    data = np.loadtxt(file_path, delimiter=",")
    np.random.seed(random_seed)
    indices = np.random.permutation(len(data))
    data = data[indices]

    X = data[:, :NUM_INPUTS]
    y = data[:, NUM_INPUTS:NUM_INPUTS + NUM_OUTPUTS]

    #Log transform x data
    X = np.log10(X)

    n_total = len(X)
    n_val = int(n_total * val_split)
    n_train = n_total - n_val

    X_train, y_train = X[:n_train], y[:n_train]
    X_val, y_val = X[n_train:], y[n_train:]

    return X_train, X_val, y_train, y_val

# -----------------------------
# Training Function
# -----------------------------
def train_model(model, train_loader, val_loader, num_epochs, lr, device,checkpoint_path,patience=20):
    model.to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    best_val_loss = float('inf')  # Initialize with a large value
    epochs_no_improve = 0

    for epoch in range(num_epochs):
        model.train()
        train_loss = 0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * X_batch.size(0)
        train_loss /= len(train_loader.dataset)

        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                outputs = model(X_batch)
                loss = criterion(outputs, y_batch)
                val_loss += loss.item() * X_batch.size(0)
        val_loss /= len(val_loader.dataset)

        print(f"Epoch [{epoch+1}/{num_epochs}] - Train Loss: {train_loss:.4f} - Val Loss: {val_loss:.4f}")
        
        # Save the model checkpoint if validation loss improved
        if val_loss < best_val_loss:
            print(f"Validation loss decreased ({best_val_loss:.4f} -> {val_loss:.4f}). Saving model...")
            best_val_loss = val_loss
            torch.save(model.state_dict(), checkpoint_path)
            epochs_no_improve = 0  # Reset counter
        else:
            epochs_no_improve += 1
            print(f"No improvement for {epochs_no_improve} epoch(s).")

        # Early stopping
        if epochs_no_improve >= patience:
            print(f"Early stopping triggered after {epoch+1} epochs. Best val loss: {best_val_loss:.4f}")
            break


# -----------------------------
# Main Execution
# -----------------------------
def main():
    # Load and prepare data
    X_train, X_val, y_train, y_val = load_and_split_data(DATA_PATH, val_split=VAL_SPLIT, random_seed=RANDOM_SEED)

    # Convert to PyTorch tensors
    X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train, dtype=torch.float32)
    X_val_tensor   = torch.tensor(X_val, dtype=torch.float32)
    y_val_tensor   = torch.tensor(y_val, dtype=torch.float32)

    # Create datasets and loaders
    train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
    val_dataset   = TensorDataset(X_val_tensor, y_val_tensor)

    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader   = DataLoader(val_dataset, batch_size=BATCH_SIZE)

    # Model
    model = Net(input_dim=NUM_INPUTS, output_dim=NUM_OUTPUTS)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Train
    train_model(model, train_loader, val_loader, num_epochs=NUM_EPOCHS, lr=LEARNING_RATE, device=device,checkpoint_path=MODEL_PATH)
    
main()
