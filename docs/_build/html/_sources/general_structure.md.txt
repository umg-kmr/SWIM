# General Structure

SWIM is organized into three core modules:

1. $G(Q)$ Calculator  
2. Semi-Analytical Power Spectrum Calculator  
3. Numerical Power Spectrum Calculator  

These modules can be used independently or combined within a workflow, depending on the use case and objectives.

---

## Overview

- The **$G(Q)$ calculator module** computes the correction factor that enters the semi-analytical power spectrum in Warm Inflation.  
- The **semi-analytical module** uses this function to efficiently evaluate the power spectrum and perform parameter inference via `Cobaya`.  
- The **numerical module** computes the full primordial power spectrum directly and also provides emulator-based acceleration for parameter inference via `Cobaya`.   

---

## $G(Q)$ Calculator

The correction function $G(Q)$ accounts for the effect of the coupling between inflaton and radiation fluctuations on the primordial power spectrum in Warm Inflation. It is obtained by numerically solving the full-set of perturbation equations.

This module is conceptually similar to existing tools such as [WI2Easy](https://github.com/RudneiRamos/WI2easy) and [WarmSpy](https://github.com/GabrieleMonte/WarmSPy).

Once computed, $G(Q)$ can be used in the semi-analytical expression for the power spectrum:

$$
P_{\mathcal{R}}(k) =
\left(\dfrac{H^2}{2\pi \dot{\phi}}\right)^2
\left[
\dfrac{T}{H}\dfrac{2\sqrt{3}\pi Q}{\sqrt{3+4\pi Q}} + 1 + 2\mathcal{N}
\right]
G(Q)
$$

where all quantities are evaluated at horizon crossing $(k = aH)$, and $\mathcal{N}$ denotes the thermal distribution of inflaton fluctuations.

### Structure

This module is implemented in:

```
GQ_Calculator/
├── bg/      # Background evolution
├── pert/    # Perturbation equations
```

It also includes Python scripts and notebooks to:
- compute initial conditions  
- evaluate $G(Q)$  
- visualize and smooth the results  

---

## Semi-Analytical Power Spectrum Calculator

This module uses the computed $G(Q)$ to evaluate the Warm Inflation (WI) power spectrum within the semi-analytical framework. It provides a computationally efficient pipeline for constraining model parameters using cosmological observations.

Two modes of operation are available, depending on the desired level of constraints and computational cost:

##### Constraints on $A_s$, $n_s$, and $r$:

In this mode, constraints are placed directly on the predicted values of $A_s$, $n_s$, and $r$. The likelihood assumes that $A_s$ and $n_s$ are uncorrelated and is implemented in `llihood.py` using a Gaussian log-likelihood.

Since this approach does not require computation of the full CMB power spectrum, it is significantly faster and well-suited for rapid parameter exploration.


##### Full CMB Constraints:

In this mode, parameter inference is performed by computing the full CMB power spectrum using `CAMB`, with the WI primordial power spectrum supplied as an external input.

This approach provides more robust and accurate constraints, as it directly compares the full CMB power spectrum with observational data. However, it is computationally more expensive than the $A_s$, $n_s$ constraints mode.

This functionality is implemented in `An_CAMB.py`.

### Structure

This module is implemented in:

```
SA_PS_Calculator/
```

It includes example configuration files for use with `Cobaya`.

## Numerical Power Spectrum Calculator

This module directly computes the WI primordial power spectrum $P_{\mathcal{R}}(k)$ by solving the full set of stochastic perturbation equations. It can be used independently of the other modules of SWIM.

Key features:

- Full stochastic-Langevin evolution of perturbations  
- Direct computation of $P_{\mathcal{R}}(k)$ as a function of $k$  
- No reliance on the $G(Q)$ approximation  

### Structure

This module is implemented in:

```
PS_Calculator/
```

It includes two sub-directories- 1) `Power_Spectrum` and 2) `Emulator`

---

### Power Spectrum

This subdirectory contains Python scripts for computing the numerical power spectrum of a WI model, along with Jupyter notebooks for visualizing and analysing the results. The outputs are stored as `.dat` files, including both the primordial power spectrum and the corresponding background evolution. These files can be used for further analysis and for extracting observables such as $A_s$, $n_s$, and $r$. This mode requires the WI model parameters to be known and explicitly specified.

### Emulator

To reduce the computational cost of parameter inference, the numerical power spectrum module includes an emulator based on Random Forest Regression (RFR).

- The emulator is trained on outputs of the numerical power spectrum solver  
- It maps model parameters to power spectrum parameters ($A_s$, $n_s$, $\alpha_s$ and $\beta_s$) 
- It is trained and used dynamically during inference to replace expensive numerical evaluations  

The `Emulator` subdirectory is organized into two components:

1. **`RF_Acc_Cobaya`**  
   Performs parameter inference using `Cobaya` with the full numerical solver. During the run, the RFR emulator is trained on-the-fly and progressively replaces calls to the computationally expensive numerical solver.

2. **`RF_Only_Cobaya`**  
   Performs parameter inference using the pre-trained emulator. This mode is highly efficient and is suitable for more detailed or extended analyses once the emulator has been trained.

```{note}
The trained emulator is specific to the Warm Inflation model under consideration. Any change in the model, including modifications to the perturbation setup (e.g., thermalization assumptions or radiation noise), requires retraining the emulator.
```

---