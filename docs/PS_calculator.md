# `PS` Calculator

This module computes the full Warm Inflation (WI) primordial power spectrum numerically using either stochastic perturbation evolution or the deterministic perturbation solver (DSWIM). It also provides parameter-inference pipelines ranging from emulator-accelerated analyses to direct full CMB likelihood constraints.

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
├── Emulator/
    ├── RF_Acc_Cobaya/
    └── RF_Only_Cobaya/
└── CMB_Const_FP/
    ├── An_CAMB.py
    └── Input.yaml
```

---

## Overview

This module has three main components:

1. **Direct Numerical Power Spectrum**
   - Computes $P_{\mathcal{R}}(k)$ for a given WI model and parameters  

2. **Emulator-Based Parameter Inference**
   - Uses the numerical solver to perform parameter inference  
   - Supports both Random Forest emulator–accelerated inference and direct deterministic (DSWIM) likelihood evaluation
         
3. **Direct CMB Parameter Inference (DSWIM)**
   - Computes the primordial power spectrum deterministically
   - Performs parameter inference through direct comparison of the resulting CMB power spectra with full observational likelihoods

---

## Direct Numerical Power Spectrum

This component computes the full numerical primordial power spectrum for a given WI model using either stochastic realizations or deterministic correlation-matrix evolution.

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

```{note}
For a detailed guide on implementing WI models, including custom dissipation coefficients, refer to the `GQ_Calculator` documentation.
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

#### Solver Selection

- `want_FP` *(int: 0 or 1, default: 1)*

  Selects the perturbation solver:

  - `0` → stochastic perturbation evolution (default SWIM solver)
  - `1` → deterministic perturbation evolution (DSWIM)

#### Numerical Controls

- `kp` *(float, default: 0.05)*  
  Pivot scale in $\mathrm{Mpc}^{-1}$

- `em_step` *(float, default: 1e-5)*  
  Minimum step size used for SDE solver

- `Nrealz` *(int, default: 2048)*  
  Number of stochastic realizations. Ignored when `want_FP = 1`  

- `kmin`, `kmax` *(float, default: -6.0, 2.0)*  
  Logarithmic range of $k$ (actual range is $10^{k_{\min}}$ to $10^{k_{\max}}$)

- `points_k` *(int, default: 50)*  
  Number of $k$ points to sample in the defined $k$-range.

    ```{note}
    For stochastic calculations, increasing `points_k` increases the total runtime approximately linearly, since each $k$-mode requires averaging over many realizations.

    For deterministic calculations (`want_FP = 1`), independent $k$-modes are evolved in parallel. Consequently, substantially larger values of `points_k` can be used with a modest impact on runtime. Values of $\mathcal{O}(10^2$--$10^3)$ are typically practical and may be desirable when constructing high-resolution primordial power spectra.
    ```

- `Np_autocalc` *(int, default: 1)*  
  Controls pivot scale exit calculation  

  ```{note}
  Automatic computation is recommended unless precise control over the pivot scale is required.
  ```

- `verbosity` *(int, default: 1)*  
  Controls logging  

---

#### Output Control

- `write_bg` *(bool, default: True)*  
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

The same notebook is applicable to both stochastic and deterministic approaches.

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
- Supports both stochastic (`want_FP = 0`) and deterministic (`want_FP = 1`) perturbation solvers.

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

This section lists all configurable parameters controlling training, usage, and reliability of the emulator. Edit the `llihood_Observables.py` script to implement them.

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

#### Training Schedule

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

#### RF Usage Control

```python
rf_uncertainty_tol = 3.0
forced_true_fraction = 0.05
```

- `rf_uncertainty_tol`  
  Controls when the RF prediction is trusted  

  The RF prediction is used if:

  $$\frac{\text{RF variance}}{\text{observational variance}} < \text{rf\_uncertainty\_tol}$$

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

#### Numerical Solver Parameters

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

The emulator framework also supports the deterministic perturbation solver (**DSWIM**) through the parameter

```python
want_FP = 1
```

When deterministic mode is enabled, the primordial power spectrum is computed using correlation-matrix evolution rather than stochastic realizations. Owing to the computational efficiency of the deterministic solver, emulator training is bypassed and parameter inference is performed through direct likelihood evaluation.

```{note}
The Random Forest emulator is used only when `want_FP = 0`. In deterministic mode (`want_FP = 1`), all likelihood evaluations are performed using DSWIM directly.
```

```{important}
The emulator-related settings are relevant only when `want_FP = 0`. They are ignored when using the deterministic solver (`want_FP = 1`).
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

```{important}
The `RF_Only_Cobaya` workflow is intended exclusively for emulator-based inference. When using DSWIM (`want_FP = 1`), `RF_Acc_Cobaya` performs direct likelihood evaluation and bypasses emulator training. As a result, no `rf_model.pkl` file is created and `RF_Only_Cobaya` cannot be used.
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
3. Choose an inference pipeline:
   - **Stochastic mode (`want_FP = 0`)**
     - Run `RF_Acc_Cobaya` to train the emulator and perform parameter inference
     - Optionally use `RF_Only_Cobaya` for subsequent fast inference using the trained emulator
   - **Deterministic mode (`want_FP = 1`)**
     - Run `RF_Acc_Cobaya` to perform direct parameter inference using DSWIM
4. Analyze chains using `getdist`  

---

### Common Pitfalls

- Incorrect function signature  
- Too few realizations (`Nrealz`) when using the stochastic solver, leading to noisy power spectra
- Poor RF training or insufficient training samples when using the emulator workflow 
- Not accounting for solver uncertainty  
- Using RF_Acc_Cobaya with multiple chains
- Attempting to use `RF_Only_Cobaya` after running in deterministic mode (`want_FP = 1`), since no emulator is trained 

---

## CMB_Const_FP — Full CMB Inference using DSWIM

This subdirectory contains a direct CMB parameter-inference pipeline based on the deterministic perturbation solver (DSWIM). The primordial power spectrum is computed deterministically and supplied to `CAMB` for comparison with observational likelihoods through `Cobaya`.

The overall workflow, input YAML configuration, likelihood installation, and chain analysis are largely identical to those described for `SA_PS_Calculator/An_CAMB.py`. As with the other submodules of `PS_Calculator`, the warm inflation model is defined in `model_calc.cpp` and inherited automatically by the inference pipeline.

### Key Differences

- Uses the full numerical primordial power spectrum rather than the semi-analytical WI spectrum.
- Uses deterministic perturbation evolution (`want_FP = 1`) by default.
- Does not require a precomputed $G(Q)$ correction function.
- Does not require emulator training.

```{important}
Before running this module, ensure that the required CMB likelihoods have been installed and configured in Cobaya, as described in `SA_PS_Calculator/An_CAMB.py`.
```

```{note}
Although DSWIM enables full CMB inference at a fraction of the cost of the stochastic approach, the semi-analytical pipeline remains the faster alternative. For warm inflation models where $G(Q)$ depends non-trivially on the model parameters, however, direct deterministic inference is generally the preferred approach. We therefore recommend running this module on HPC systems with multiple CPU threads.
```
---

## Summary

- Computes full numerical WI power spectrum
- Most accurate but computationally expensive  
- Emulator enables efficient inference
- DSWIM provides alternative efficient and accurate approach 
- Modular design 