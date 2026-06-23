# Quick Start

This guide walks through the core workflow of SWIM using a built-in example Warm Inflation (WI) model:

$$
V(\phi) = \dfrac{1}{4} V_0 \phi^4, \qquad
\Upsilon(\phi,T) = C_{\Upsilon} T^3
$$

---

## Overview

SWIM provides two complementary workflows:

### Semi-Analytical Workflow

1. Compute $G(Q)$
2. Perform parameter inference using the semi-analytical primordial power spectrum

### Fully Numerical Workflow

1. Compute the primordial power spectrum using stochastic or deterministic perturbation evolution
2. Perform parameter inference using either emulator-accelerated (stochastic) or direct likelihood evaluation (deterministic)
3. Optionally perform full CMB likelihood inference using deterministic framework

---

## 1. Compute $G(Q)$

Navigate to the module:

```bash
cd GQ_Calculator
```

Clean previous outputs:

```bash
rm ics.dat GQ.dat GQ_smooth.dat
```

Compute initial conditions:

```bash
python -u find_ICs.py
```

Then compute $G(Q)$:

```bash
python -u find_GQ.py
```

**Deterministic mode (DSWIM):** To compute $G(Q)$ using the deterministic perturbation solver, set

```python
want_FP = 1
```

in `find_GQ.py` before execution. In this mode, the perturbation correlation matrix is evolved deterministically using the Fokker--Planck formalism. The default stochastic implementation is recovered by setting `want_FP = 0`.

Outputs:
- `ics.dat` — initial conditions  
- `GQ.dat` — raw $G(Q)$  

Use the notebook `GQ_Plotting_NB.ipynb` to visualize and smooth the output (`GQ_smooth.dat`).

---

## 2. Semi-Analytical Inference 

Requires `Cobaya` and cosmological likelihoods.

Navigate to:

```bash
cd ../SA_PS_Calculator
```

Remove old chains:

```bash
rm -rf chains_Asns chains_fullCMB
```

Run inference (example with 8 chains):

```bash
mpirun -n 8 cobaya-run Input_asns.yaml
```

For full CMB likelihoods:

```bash
mpirun -n 8 cobaya-run Input.yaml
```

```{note}
When running multiple chains, limit CPU usage per chain:

    export OMP_NUM_THREADS=$(( $(nproc --all) / 8 ))
```

```{important}
Make sure to set the `packages_path:` in `Input.yaml` to the location of your Cobaya installation. Otherwise, Cobaya will not be able to locate external codes (e.g. CAMB) and likelihoods.
```

---

## 3. Numerical Power Spectrum

Compute the full numerical WI power spectrum:

```bash
cd ../PS_Calculator/Power_Spectrum
```

Clean outputs:

```bash
rm bg.dat ps.dat PT_kp.dat
```

use all CPU threads:

```bash
export OMP_NUM_THREADS=$(nproc --all)
```

Run:

```bash
python -u ps_script.py
```

**Deterministic mode (DSWIM):** To use the deterministic perturbation solver, set

```python
want_FP = 1
```

in `ps_script.py` before execution. 

Outputs:
- `bg.dat` — background evolution  
- `ps.dat` — power spectrum  
- `PT_kp.dat` — tensor amplitude  

Use `Plotting_NB.ipynb` to visualize results and extract $(A_s, n_s, r)$.

---

## 4. Emulator-Based Inference 

Accelerated inference using a Random Forest emulator:

```bash
cd ../../Emulator/RF_Acc_Cobaya
```

Remove previous chains:

```bash
rm -rf chains*
```

Run:

```bash
cobaya-run Input_asns.yaml
```

After ~100 valid samples, the emulator is trained and saved as: `rf_model.pkl`

```{note}
The emulator is trained only for the chosen WI model. Any change in the model or perturbation settings requires retraining.
```

**Deterministic mode (DSWIM):** To use the deterministic perturbation solver, set

```python
want_FP = 1
```

in `llihood_Observables.py` before running Cobaya. Owing to the computational efficiency of the deterministic solver, parameter inference is performed directly and the Random Forest emulator is bypassed.

---

## 5. CMB Inference using DSWIM

Full CMB parameter inference using the deterministic perturbation solver:

```bash
cd ../../PS_Calculator/CMB_Const_FP
```

Remove previous chains:

```bash
rm -rf chains*
```

Run Cobaya (example with 8 chains):

```bash
mpirun -n 8 cobaya-run Input_asns.yaml
```

```{note}
When running multiple chains, limit CPU usage per chain:

    export OMP_NUM_THREADS=$(( $(nproc --all) / 8 ))
```

This module uses the deterministic perturbation solver (**DSWIM**) by default. Although the stochastic solver can also be used, it is generally not recommended for parameter inference owing to its significantly higher computational cost. This module performs direct likelihood evaluation using the full CMB angular power spectra and corresponding observational likelihoods.

Before running this module, ensure that all required CMB likelihoods have been installed and configured in Cobaya.


---

## Summary

- Use **GQ module** → computes correction factor using stochastic or deterministic perturbation evolution 
- Use **SA module** → fast semi-analytical inference  
- Use **PS module - Power_Spectrum** → full numerical spectrum using stochastic or deterministic approaches 
- Use **PS module - Emulator** → efficient parameter inference using numerical power spectrum (stochastic/deterministic) 
- Use **PS module - CMB_Const_FP** → full CMB likelihood parameter inference using DSWIM (no emulator) 

---

For full details of each module, see the corresponding sections in the documentation.