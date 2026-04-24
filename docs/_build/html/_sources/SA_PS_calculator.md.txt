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
  - optionally $\alpha_s$  
- The tensor-to-scalar ratio $r$ is constrained internally (upper bound)  

This script is used with `Input_asns.yaml`.

---

### `An_CAMB.py`

- Injects the WI primordial power spectrum into CAMB  
- Computes CMB angular power spectra  
- Used for full CMB likelihood analysis  

---

## $G(Q)$ Handling

The semi-analytical power spectrum requires $G(Q)$, which can be supplied in multiple ways.

---

### Option 1: Precomputed $G(Q)$ (Recommended)

Set:

```python
read_GQ_from_file = 1
```

- Uses precomputed data (e.g. from SWIM GQ calculator)  
- File format:
  ```
  Q   G(Q)
  ```
- Interpolation is handled internally in C++  

```{note}
It is strongly recommended to smooth $G(Q)$ using the notebook provided in the GQ module.
```

---

### Option 2: Internal Analytical Fits

Set:

```python
read_GQ_from_file = 0
```

The code uses built-in fitting functions in `model_calc.cpp`:

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

```{note}
These fits are available only for specific cases (e.g. $p=3,1,-1$).  
Users must modify `model_calc.cpp` to add new fits if required.
```

---

### External Sources

$G(Q)$ can also be generated using:

- SWIM GQ Calculator  
- WI2Easy  
- WarmSPy  

```{warning}
Ensure that the $G(Q)$ used corresponds to the same WI model being analyzed.
```

---

## C++–Python Interface

As in the $G(Q)$ module, the function signature must match between C++ and Python.

Example:

```cpp
void model(double phi_ini,double gst,double Q_ini,double V0,double Np,int c,int p,int therm);
```

This must be consistent in:

- `model_calc.cpp`  
- `llihood.py`  
- `An_CAMB.py`  

```{important}
Any mismatch in function signature will lead to runtime errors or incorrect results.
```

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

Refer to the Cobaya documentation for installation.

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
- Not smoothing $G(Q)$ (can introduce noise in inference)  
- Insufficient $Q$ range in $G(Q)$ data  
- Function signature mismatch between C++ and Python  
- CAMB overusing threads in multi-chain runs  

```{note}
When running multiple chains, limit CAMB threads using:
```

```bash
export OMP_NUM_THREADS=...
```

---

## Summary

- Fast parameter inference using semi-analytical power spectrum  
- Requires $G(Q)$ as input  
- Supports both simple ($A_s$, $n_s$) and full CMB analyses  
- Designed for efficiency and flexibility  