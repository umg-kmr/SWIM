# `SA_PS` Calculator

This module provides an efficient framework for performing parameter inference in Warm Inflation (WI) models using a semi-analytical power spectrum.

It relies on the correction function $G(Q)$ and avoids solving the full stochastic perturbation equations, making it significantly faster than the numerical module.

---

## Directory Structure

```
SA_PS_Calculator/
├── Bg.cpp
├── model_calc.cpp
├── llihood.py
├── An_CAMB.py
├── Input_asns.yaml
└── Input.yaml
```

---

## Overview

This module computes the WI power spectrum using:

- Numerical background evolution  
- Semi-analytical expressions for the scalar power spectrum  
- Precomputed or analytical forms of $G(Q)$  

It is designed primarily for **parameter inference using Cobaya**.

---

## Defining Your WI Model

The WI model is defined in:

```
model_calc.cpp
```

The structure is identical to the $G(Q)$ module, with the following simplifications:

- No need to define $V''(\phi)$  
- No need to define derivatives of $\Upsilon(\phi,T)$  

Only the following functions are required:

- $V(\phi)$  
- $V'(\phi)$  
- $\Upsilon(\phi, T)$ *(only if not of the form $T^p \phi^c$)*  

If the dissipation coefficient is of the form $\Upsilon \propto T^p \phi^c$, it is already implemented internally and the user only needs to specify the parameters `p` and `c` in the Python scripts. For more general forms of $\Upsilon(\phi, T)$, the user must explicitly define it in `model_calc.cpp`. 

```{important}
Any modification to `model_calc.cpp` requires recompilation.
```

---

## Role of Core Files

### `Bg.cpp`

- Solves background evolution  
- Determines pivot scale exit (if enabled)  
- Computes the semi-analytical power spectrum  

---

### `llihood.py`

- Implements likelihood for:
  - $A_s$, $n_s$  
  - optionally $\alpha_s$ (*running of the spectral index*) 
- The tensor-to-scalar ratio $r$ is constrained internally (upper bound) in `Bg.cpp`. 

This script is used with `Input_asns.yaml`.

---

### `An_CAMB.py`

- Injects the WI primordial power spectrum into CAMB  
- Used for computing the CMB angular power spectra using CAMB 
- Performs full CMB likelihood analysis with Cobaya

---

## $G(Q)$ Handling

The semi-analytical power spectrum requires $G(Q)$, which can be supplied in multiple ways.

---

### Option 1: External $G(Q)$ (Recommended)

Set:

```python
read_GQ_from_file = 1
```

- Uses precomputed $G(Q)$ data  
- File format:
  ```
  Q   G(Q)
  ```
- Interpolation is handled internally in C++  

$G(Q)$ can be generated using:
- SWIM GQ Calculator  
- WI2Easy  

```{note}
If using SWIM’s GQ Calculator, it is strongly recommended to smooth $G(Q)$ using the notebook provided in the GQ module before using it here.

If using a $G(Q)$ file computed with WI2Easy, update the path to the file in `model_calc.cpp` by modifying the parameter `GQfname` of type `string`.
```

---

### Option 2: Internal Analytical Fits

Set:

```python
read_GQ_from_file = 0
```

The code uses analytical fitting functions defined in `model_calc.cpp`:

```cpp
auto GQ = [p](double Q) -> double {
    if (GQ_from_file == 1) {
        return exp10(logGasQ(log10(Q)));
    }
    else {
        // internal fits depending on p
    }
};
```

- Several commonly used fits are already implemented  
- These are available only for specific cases (e.g. $p = 3, 1, -1$)  

Additional fitting functions can be:
- computed using WarmSPy  
- implemented by modifying the `auto GQ` function  


---

```{warning}
Ensure that the $G(Q)$ used corresponds to the same WI model being analyzed.
```

---

## C++ – Python Interface

As in the $G(Q)$ module, the function signature must match between C++ and Python.

Example:

```cpp
void model(double phi_ini,double gst,double Q_ini,double V0,double Np,int c,int p,int therm);
```

This must be consistent in:

- `model_calc.cpp`  
- `llihood.py`  (*modify the function `logp`*)
- `An_CAMB.py`  (*modify the function `feature_power_spectrum`*)

```{important}
Any mismatch in function signature will lead to runtime errors or incorrect results.
```
### Full CMB Parameter Interface (`An_CAMB.py`)

When using full CMB constraints, the WI model parameters must also be defined in the `FeaturePrimordialPk` class inside `An_CAMB.py`.

This step is required for Cobaya to correctly recognize and pass the model parameters.

---

#### What Needs to Be Modified

Update the following section to include your WI model parameters:

```python
class FeaturePrimordialPk(Theory):
    # define parameter names here
    params = {"phi0":None,"gst":None,"Q0":None,"V0":None,"Np":None,"c":None,"p":None,"therm":None}

    def calculate(self, state, want_derived=True, **params_values_dict):
        # extract parameters here
        phi0,gst,Q0,V0,Np,c,p,therm = \
        [params_values_dict[itr] for itr in ["phi0","gst","Q0","V0","Np","c","p","therm"]] 

        ks, Pks = feature_power_spectrum(phi0,gst,Q0,V0,Np,c,p,therm,kp=self.kp) #match function signature

```

---

#### Key Points

- The `params` dictionary must include **all WI model parameters**  
- The order of parameters extracted in `calculate()` must match:
  - the `params` dictionary  
  - the C++ function signature  
- The call to `feature_power_spectrum(...)` must also follow the same ordering  

```{important}
If you add or remove model parameters, you must update:
- `feature_power_spectrum` arguments list
- the `params` dictionary  
- the parameter extraction list  
- the function call to `feature_power_spectrum(...)`
```

---

## Input Parameters

All runtime parameters are configured in the Python scripts.

---

## `llihood.py` — $A_s$, $n_s$ Inference

### Pivot Scale Handling

- `Np_autocalc` *(int, default = 1)*  
  Controls how the pivot scale exit is determined.

  - `1` → Automatically compute the e-fold $N_*$ at which the pivot scale exits the horizon  
  - `0` → Use user-specified value of `Np`  

  ```{note}
  The integration starts from $N = 0$, so `Np` must be specified consistently in this convention when `Np_autocalc = 0`.

  Automatic computation is recommended unless precise control over the pivot scale is required.
  ```

- `kp` *(float, default = 0.05)*  
  Pivot scale in $\mathrm{Mpc}^{-1}$.

---

### $G(Q)$ Handling

- `read_GQ_from_file` *(int, default = 1)*  
  - `1` → Use precomputed $G(Q)$ from file  
  - `0` → Use internal analytical fits  

---

### Spectrum Control

- `want_full_spectrum` *(int, do not modify default = 0)*  
  - `0` → Compute only $A_s$, $n_s$  
  - `1` → Compute full power spectrum (not used in this mode)

---

### Numerical Parameters

- `verbosity` *(int, default = 0)*  
  Controls logging output. 

---

### Auxiliary Parameters

*(Required by the interface but not used in this mode)*

- `kmin` *(float, default = 1e-6)*  
- `kmax` *(float, default = 100.0)*  

---

## `An_CAMB.py` — Full CMB Inference

### Pivot Scale Handling

- `Np_autocalc` *(int, default = 1)*  
  Same meaning as in `llihood.py`.

---

### $G(Q)$ Handling

- `read_GQ_from_file` *(int, default = 1)*  
  Same meaning as in `llihood.py`.

---

### Spectrum Control

- `want_full_spectrum` *(int, do not modify default = 1)*  
  - `1` → Required for full CMB computation  
  - `0` → Not suitable for this mode  

---

### Numerical Parameters

- `verbosity` *(int, default = 0)*  
  Controls logging output.

---

## Running Parameter Inference

### Mode 1: $A_s$, $n_s$ (fast)

```bash
cobaya-run Input_asns.yaml
```

- Uses `llihood.py`  
- No full CMB computation  
- Fast and efficient  

---

### Mode 2: Full CMB Constraints

```bash
cobaya-run Input.yaml
```

- Uses `An_CAMB.py`  
- Requires:
  - CAMB  
  - likelihoods (Planck, ACT, SPT, DESI)  

Refer to the [Cobaya documentation](https://cobaya.readthedocs.io/en/latest/likelihood_external.html) for installation.

```{important}
Make sure to set the `packages_path:` in `Input.yaml` to the location of your Cobaya installation. Otherwise, Cobaya will not be able to locate external codes (e.g. CAMB) and likelihoods.
```

---

## Workflow

1. Define WI model in `model_calc.cpp`  
2. Ensure function signatures match in Python files  
3. Provide or compute $G(Q)$  
4. Configure input YAML files (priors, parameters)  
5. Run Cobaya  
6. Analyze chains using `getdist`  

---

## Common Pitfalls

- Using incorrect or inconsistent $G(Q)$  
- Not smoothing $G(Q)$ computed by SWIM (can introduce noise in inference)  
- Insufficient $Q$ range in $G(Q)$ data  
- Function signature mismatch between C++ and Python  
- CAMB overusing threads in multi-chain runs  

```{Important}
When running multiple chains (`mpirun`), limit CAMB threads using:


    export OMP_NUM_THREADS=$(( $(nproc --all) / <number-of-chains> ))
```

---

## Summary

- Fast parameter inference using semi-analytical power spectrum  
- Requires $G(Q)$ as input  
- Supports both simple ($A_s$, $n_s$) and full CMB analyses  
- Designed for efficiency and flexibility  