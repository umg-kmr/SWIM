# Quick Start

This guide walks through the core workflow of SWIM using a built-in example Warm Inflation (WI) model:

$$
V(\phi) = \dfrac{1}{4} V_0 \phi^4, \qquad
\Upsilon(\phi,T) = C_{\Upsilon} T^3
$$

---

## Overview

The typical workflow in SWIM is:

1. Compute the correction function $G(Q)$  
2. Perform semi-analytical parameter inference using the computed $G(Q)$
   
**or**

1. Bypass the semi-analytical approach and compute the full numerical power spectrum 
2. Perform parameter inference using full numerical power spectrum accelerated using a RandomForest regression emulator  

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

Compute the full numerica WI power spectrum:

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
rm -rf chains
```

Run:

```bash
cobaya-run Input_asns.yaml
```

After ~100 valid samples, the emulator is trained and saved as: `rf_model.pkl`

```{note}
The emulator is trained only for the chosen WI model. Any change in the model or perturbation settings requires retraining.
```

---

## Summary

- Use **GQ module** → compute correction factor  
- Use **SA module** → fast semi-analytical inference  
- Use **PS module** → full numerical spectrum  
- Use **PS module - Emulator** → efficient parameter inference using numerical power spectrum  

---

For full details of each module, see the corresponding sections in the documentation.