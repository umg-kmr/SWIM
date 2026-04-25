# `PS` Calculator

This module computes the full Warm Inflation (WI) primordial power spectrum numerically using a stochastic approach. It also provides a framework for parameter inference using a Random Forest Regression (RFR) emulator.

This is the most complete and physically accurate module in SWIM, as it does not rely on semi-analytical approximations and other modules of SWIM.

---

## Directory Structure

```
PS_Calculator/
├── model_calc.cpp
├── Bg.cpp
├── Power_Spectrum/
│   ├── ps_script.py
│   ├── functions_bg_diag.py
│   └── Plotting_NB.ipynb
└── Emulator/
    ├── RF_Acc_Cobaya/
    └── RF_Only_Cobaya/
```

---

## Overview

This module has two main components:

1. **Direct Numerical Power Spectrum**
   - Computes $P_{\mathcal{R}}(k)$ for a given WI model and parameters  

2. **Emulator-Based Parameter Inference**
   - Uses the numerical solver to perform parameter inference  
   - Accelerates computation using a Random Forest surrogate  

---

## Direct Numerical Power Spectrum

This component computes the full numerical power spectrum for a given WI model.

---

### Defining Your WI Model

The model must be defined in:

```
model_calc.cpp
```

As in the GQ module, define:

- $V(\phi)$  
- $V'(\phi)$  
- $V''(\phi)$  
- $\Upsilon(\phi, T)$ *(if not of the form $T^p \phi^c$)* 
- Derivatives of $\Upsilon$ *(if not of the form $T^p \phi^c$)*  

```{important}
Any modification to `model_calc.cpp` requires recompilation.
```

---

### Running the Solver

Navigate to:

```bash
cd PS_Calculator/Power_Spectrum
```

Run:

```bash
python -u ps_script.py
```

---

### Input Parameters (`ps_script.py`)

#### Numerical Controls

- `kp` *(float, default = 0.05)*  
  Pivot scale in $\mathrm{Mpc}^{-1}$

- `em_step` *(float, default = 1e-5)*  
  Minimum step size used for SDE solver

- `Nrealz` *(int, default = 2048)*  
  Number of stochastic realizations  

- `kmin`, `kmax` *(float, default = -6.0, 2.0)*  
  Logarithmic range of $k$ (actual range is $10^{k_{\min}}$ to $10^{k_{\max}}$)

- `points_k` *(int, default = 50)*  
  Number of $k$ points to sample in the defined $k$-range 

- `Np_autocalc` *(int, default = 1)*  
  Controls pivot scale exit calculation  

  ```{note}
  Automatic computation is recommended unless precise control over the pivot scale is required.
  ```

- `verbosity` *(int, default = 1)*  
  Controls logging  

---

#### Output Control

- `write_bg` *(bool, default = True)*  
  Save background evolution  

- `fname_bg`, `fname_ps`  
  Output file names. If changed here then should be modified in the plotting notebook as well consistently.  

---

#### Model Parameters

Example:

```python
V0 = ...
gst = ...
Q0 = ...
phi0 = ...
Np = ...
therm = ...
rad_noise = ...
p = ...
c = ...
```

```{important}
These must match the function signature in `model_calc.cpp`.
```

---

### Output Files

- `bg.dat`

  $N \quad \phi \quad \phi' \quad T$

- `ps.dat`

  $k \quad  P_{\mathcal{R}}(k)$

- `PT_kp.dat`  
  - Tensor amplitude at pivot scale  

---

### Analysis and Plotting

The notebook:

```
Plotting_NB.ipynb
```

- Plots the raw power spectrum  
- Fits the spectrum to extract:
  - $A_s$, $n_s$, $\alpha_s$, $\beta_s$  
- Computes $r$  
- Plots background evolution (if enabled)

---

#### Power Spectrum Fitting

The following functional form is used:

```python
def fitting_fn(lnk, lnAs, ns, alphs, betas):
    return lnAs + (lnk - np.log(kp)) * (
        ns - 1
        + 0.5 * alphs * (lnk - np.log(kp))
        + (1/6) * betas * (lnk - np.log(kp))**2
    )
```

This corresponds to a power-law spectrum with running.

---

#### `functions_bg_diag.py`

This file is used **only for plotting**.

```{note}
The WI model must be redefined here so that Python can evaluate background quantities for visualization. This does not affect the numerical solver.
```

---

## Emulator-Based Parameter Inference

This module enables parameter inference using the numerical solver.

---

### RF_Acc_Cobaya (Accelerated Inference)

- Uses full numerical solver during MCMC  
- Trains a Random Forest model on-the-fly  
- Gradually replaces expensive solver evaluations  

```{note}
This mode must be run with a single chain.
```

---

#### Key Features

- Computes and fits numerical $P(k)$ → extracts observables  
- Compares with observational constraints using Gaussian likelihood  
- Incorporates solver uncertainty into likelihood  

---

#### Solver Uncertainty

The script:

```
std_gen.py
```

estimates the intrinsic stochastic variance of the solver by repeated evaluations at the same parameters.

This contributes to the total uncertainty in:

- $A_s$, $n_s$, $\alpha_s$, $\beta_s$

```{note}
The result of this script should be used to modify `yerr_solver` in `llihood_Observables.py`.
```

---

### RF Training and Control Parameters

The Random Forest (RF) emulator is trained dynamically during MCMC to approximate the mapping:

$$\{\text{model parameters}\} \rightarrow \{A_s, n_s, \alpha_s, \beta_s\}$$

This section lists all configurable parameters controlling training, usage, and reliability of the emulator.

---

#### Model Loading

```python
load_previous_rf = False
```

- If `True`, loads an existing trained model (`rf_model.pkl`)
- Useful for:
  - restarting runs  
  - continuing training  
  - using pre-trained emulator  

---

### Training Schedule

```python
update_frequency = 100
min_points_before_rf = 100
max_training_points = 3000
```

- `min_points_before_rf`  
  Minimum number of **true solver evaluations** before the RF is first trained  

- `update_frequency`  
  RF is retrained every `update_frequency` true evaluations  

- `max_training_points`  
  Maximum number of training samples retained  

---

#### Sliding Window Training

When the number of training samples exceeds `max_training_points`, only the **most recent samples** are retained.

This is referred to as a *sliding window*:

- Old data points are discarded  
- New data continuously replaces them  
- Keeps training focused on the **high-likelihood region**  

```{note}
Sliding window training improves local accuracy near the posterior peak and avoids memory growth.
```

---

### RF Usage Control

```python
rf_uncertainty_tol = 3.0
forced_true_fraction = 0.05
```

- `rf_uncertainty_tol`  
  Controls when the RF prediction is trusted  

  The RF prediction is used if:

  $$
  \frac{\text{RF variance}}{\text{observational variance}} < \text{rf_uncertainty_tol}
  $$

  - Lower values → more conservative (more true evaluations)  
  - Higher values → more aggressive RF usage (faster, less accurate)  

- `forced_true_fraction`  
  Fraction of samples where the **true numerical solver is forced**, even if RF is trusted. For example, if `forced_true_fraction = 0.05`, then the true solver is used in approximately 5% of evaluations regardless of RF confidence.

  Purpose:
  - Prevent RF drift  
  - Maintain training quality  

```{note}
Typical values: 0.01–0.1.  
Higher values improve robustness but increase runtime.
```

---

### Numerical Solver Parameters

These must also be set consistently (same as direct PS computation):

```python
kp = 0.05
em_step = 1e-5
Nrealz = 2048
kmax = np.log10(1e2)
kmin = np.log10(1e-6)
points_k = 50
Np_autocalc = 1
verbosity = 0
```
---

### Practical Recommendations

#### For Initial Runs

```python
min_points_before_rf = 100
update_frequency = 100
rf_uncertainty_tol = 3-4
forced_true_fraction = 0.05
```

---

#### For Faster Runs 

```python
rf_uncertainty_tol = 5
forced_true_fraction = 0.01
```

---

#### For Maximum Accuracy (increased runtimes)

```python
rf_uncertainty_tol = 1-2
forced_true_fraction = 0.1
```

---

### Summary

- RF training and usage starts after sufficient true evaluations  
- Uses sliding window to focus on relevant parameter space  
- Dynamically balances:
  - accuracy (true solver)  
  - speed (RF emulator)  
- Controlled via uncertainty threshold `rf_uncertainty_tol` and forced sampling `forced_true_fraction`  
---

### RF_Only_Cobaya (Emulator Only Inference)

- Uses pre-trained model (`rf_model.pkl`)  
- Does not call numerical solver  
- Extremely fast  

```{note}
This mode is useful for validating emulator performance and running extended analyses (full CMB likelihood).
```

---

### Limitations

- Full CMB likelihood not yet implemented  
- RF_Acc_Cobaya supports only single-chain runs
- Parameter inference is still computationally expensive and initial burn-in can take a long time

```{Important}
It is recommended to run the numerical solver based inference on HPCs with many CPU threads enabled.
```

```{note}
Full CMB inference can be implemented by adapting the approach in `SA_PS_Calculator/An_CAMB.py`.
```

---

### Workflow

1. Define WI model in `model_calc.cpp`  
2. (Optional) Run `ps_script.py` to inspect spectrum  
3. Run RF_Acc_Cobaya to train emulator  
4. Use RF_Only_Cobaya for fast inference  
5. Analyze chains using `getdist`  

---

### Common Pitfalls

- Incorrect function signature  
- Too few realizations → noisy spectrum  
- Poor RF training → inaccurate inference  
- Not accounting for solver uncertainty  
- Using RF_Acc_Cobaya with multiple chains  

---

## Summary

- Computes full numerical WI power spectrum  
- Most accurate but computationally expensive  
- Emulator enables efficient inference  
- Modular design 