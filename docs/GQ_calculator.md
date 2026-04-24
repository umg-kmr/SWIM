# `GQ` Calculator

This module computes the correction function $G(Q)$ for a given Warm Inflation (WI) model by solving the background and perturbation equations numerically.

It is organized into two main components `bg` and `pert` with helper Python scripts:

```
GQ_Calculator/
├── bg/                # Background evolution
├── pert/              # Perturbations + power spectrum
├── find_ICs.py        # Initial condition finder
├── find_GQ.py         # Computes G(Q)
└── GQ_Plotting_NB.ipynb  # Plotting and Smoothing G(Q)
```

---

## Overview

The computation proceeds in two stages:

1. **Initial condition search (`find_ICs.py`)**  
   Determines the field value $\phi_{\text{initial}}$ corresponding to a given $Q_{\text{initial}}$ such that a desired duration of inflation is obtained.

2. **$G(Q)$ computation (`find_GQ.py`)**  
   Uses these initial conditions to solve the background as well as perturbation equations and computes:

$$G(Q_*) = \left.\frac{P_{\mathrm{numerical}}}{P_{\mathrm{analytical}}}\right|_{Q = Q_*}$$

---

## Defining Your WI Model

The WI model must be defined in C++ in:

```
bg/model_calc.cpp
pert/model_calc.cpp
```

### Required Functions

The user must define the following functions using C++ lambda expressions:

- Potential $V(\phi)$: `V(phi)`
- First derivative $V_{,\phi}(\phi)$: `Vd(phi)`
- (Perturbation module only) Second derivative $V_{,\phi\phi}(\phi)$: `Vdd(phi)`
- Dissipation coefficient $\Upsilon(\phi,T)$: `Ups_wo_Cy(phi,T)`

Lambda functions are used throughout SWIM for flexibility, allowing model definitions to be modified without changing the solver structure.

### Lambda Function Template

In SWIM, these functions are implemented using **C++ lambda functions**, which allow model parameters to be directly captured inside the function definition.

General template:

```cpp
auto function_name = [captured_parameters](input_variables) -> return_type {
    // function definition
};
```

- `auto` lets the compiler infer the function type  
- `[captured_parameters]` allows passing model parameters (e.g. $V_0$) into the function  
- `(input_variables)` are the function arguments (e.g. $\phi$, $T$)  

---

### Example: Defining the Potential

```cpp
auto V = [V0](double phi) -> double {
    return (V0/4.0)*pow(phi,4.0);
};
```

This defines the quartic potential:

$$V(\phi) = \frac{V_0}{4}\phi^4$$

---

### Example: First Derivative

```cpp
auto Vd = [V0](double phi) -> double {
    return V0*phi*phi*phi;
};
```

which corresponds to:

$$V'(\phi) = V_0 \phi^3$$

### Dissipation Coefficient

If the dissipation coefficient has the form: 

$$\Upsilon(\phi, T) \propto T^p \phi^c$$

then no modification in the C++ code is required. The functional form is already implemented internally, and the user only needs to specify the parameters `p` and `c` in the Python scripts.

The constant $C_{\Upsilon}$ is computed internally using the initial conditions $(\phi_\text{initial}, Q_\text{initial})$.

---

For more general forms of $\Upsilon(\phi, T)$, the user must modify the function:

```cpp
auto Ups_wo_Cy = [custom_parameters](double phi,double T) -> double {
    // define custom form here
};
```

In this case, the corresponding derivatives of $\Upsilon$ must also be specified in the `pert/model_calc.cpp` file.

---

```{important}
The model must be defined consistently in both `bg/model_calc.cpp` and `pert/model_calc.cpp`.
```
---

## C++ – Python Interface

The WI model is exposed from C++ to Python through a C-compatible function:

```cpp
extern "C" {
    void model(double phi_ini, double Q_ini, double gst,
               double V0, int p, int c, int hybrid_inf);
}
```

---

### “Matching The Function Signature”

The *function signature* refers to:

- Function name (`model`)  
- Number of arguments  
- Order of arguments  
- Type of each argument (`double`, `int`, etc.)  

All of these must be **identical everywhere** the function is defined or called.

---

### Where Must It Match?

The same function signature must appear in:

- `bg/model_calc.cpp`  
- `pert/model_calc.cpp`  
- `find_ICs.py` (via `ffi.cdef`)  
- `find_GQ.py` (via `ffi.cdef`)  

If any of these differ, the code will either:
- fail at runtime, or  
- produce incorrect results  

---

### When Do You Need to Modify It?

You must update the function signature if your WI model introduces new parameters.

#### Example

Suppose your model potential $V(\phi)$ depends on an additional parameter $\alpha$.

---

#### Step 1: Modify C++ (`model_calc.cpp`)

```cpp
extern "C" {
    void model(double phi_ini, double Q_ini, double gst,
               double V0, double alpha, int p, int c, int hybrid_inf);
}
```

---

#### Step 2: Use $\alpha$ in the Model Definition

Update your lambda functions (`V(phi),Vd(phi)`) to capture the new parameter:

```cpp
auto V = [V0, alpha](double phi) -> double {
    return V0 * exp(-alpha * phi);
};
```

Here:
- `alpha` is included in the capture list `[V0, alpha]`
- This makes it accessible inside the function

```{important}
If a parameter is not included in the capture list `[ ... ]`, it will not be available inside the lambda function and will lead to a compilation error.
```

---

```{note}
In the `pert` module as well, you must update the lambda functions (`V(phi),Vd(phi)` and `Vdd(phi)`) consistently with the modified potential.
```

---

#### Step 3: Update Python Binding (`find_ICs.py`)

```python
ffi.cdef("void model(double phi_ini,double Q_ini,double gst,double V0,double alpha,int p,int c,int hybrid_inf);")
```

---

#### Step 4: Update Python Binding (`find_GQ.py`)

```python
ffi.cdef("void model(double phi_ini,double Q_ini,double gst,double V0,double alpha,int p,int c,int therm,int rad_noise,int hybrid_inf);")
```

---

#### Step 5: Update Function Calls

```python
lib.model(phi0, Q0, gst, V0, alpha, p, c, hybrid_inf)
```

---

```{important}
The number, order, and type of arguments must be identical in all locations. Even a small mismatch (e.g., missing parameter or wrong order) will break the interface.
```

---

### Key Rules (Safe Checklist)

- Always copy the argument list from `model_calc.cpp` into Python  
- Keep the same order of parameters everywhere  
- Use the same data types (`double`, `int`)  
- Update both `find_ICs.py` and `find_GQ.py`
- Ensure new parameters are captured in lambda functions  
---

```{important}
Any modification to `.cpp` files requires recompilation of SWIM.
```

---

## Input Parameters

All runtime parameters are configured in the Python scripts.  
The parameters are grouped below by functionality.

---

## `find_ICs.py` — Initial Condition Setup

### Model Parameters

These define the underlying Warm Inflation model.

- `V0` *(float)*  
  Overall normalization of the potential ($V_0$).

- `gst` *(float)*  
  Effective number of relativistic degrees of freedom ($g_*$).

- `p, c` *(int, optional)*  
  Exponents for $\Upsilon \propto T^p \phi^c$.  
  Only used if the dissipation coefficient follows this form.

- `hybrid_inf` *(int: 0 or 1)*  
  Disables (0)/Enables (1) hybrid inflation mode.  
  If set to 1, the critical field value must also be specified.

---

### Inflation Control

- `dur_N` *(float, default: 60.0)*  
  Desired duration of inflation in e-folds.

---

### $Q$ Sampling

- `Qlow`, `Qup` *(float)*  
  Lower and upper bounds of $Q$ range to be explored.

- `npts` *(int, default: 150)*  
  Number of points used to sample the $Q$ range.

---

### Root-Finding Parameters

Used in the brute-force search for initial field values.

- `lower_bound`, `upper_bound` *(float)*  
  Defines the search interval for $\phi$ (in log-space) and should be set according to the WI model being studied.  
  Wider ranges increase robustness but also runtime.

---

### Parallelization

- `Nprocs` *(int)*  
  Number of CPU threads used for parallel root finding. Set depending upon your system, `nprocs --all` lists all available threads.

---

## `find_GQ.py` — $G(Q)$ Computation

### Model Parameters

Same as in `find_ICs.py`.  
These must be kept consistent.

---

### Physical Options

- `therm` *(int: 0 or 1)*  
  Disables (0)/ Enables(1) thermalization effects in perturbations.

- `rad_noise` *(int: 0 or 1)*  
  Excludes (0)/Includes (1) stochastic thermal noise in radiation perturbations.

---

### Numerical Controls

- `Nrealz` *(int, default: 2048)*  
  Number of stochastic realizations used for averaging.

- `Nstar` *(float, default: 7)*  
  E-fold at which the power spectrum is evaluated. Set it such that enough e-folds are available for the mode to initialize and freeze.

- `verbosity` *(int, default: 0)*  
  Controls logging output (0 = silent).

```{note}
The number of CPU threads used to parallelize over independent stochastic realizations is controlled by the environment variable `OMP_NUM_THREADS`. By default, all available threads are used. Set this environment variable to limit the threads used by the `GQ_Calculator`.

The runtime scales approximately inversely with the number of threads (up to hardware limits), since realizations are computed independently. 
```

---

```{note}
The parameters in `find_ICs.py` and `find_GQ.py` must be consistent.  
Mismatch in model parameters will lead to incorrect results.
```

## Running the Pipeline

```bash
python -u find_ICs.py
python -u find_GQ.py
```

---

## Output Files

- `ics.dat`: $(\phi_{\text{initial}}, Q_\text{initial})$  
- `GQ.dat`: $(Q_*, G(Q_*))$  
- `GQ_smooth.dat`: smoothed version (via notebook)

---

## Methodology

### Initial Conditions

The initial field value $\phi_\text{initial}$ is obtained via a **brute-force root-finding procedure** over a user-defined interval, ensuring a fixed duration of inflation (`dur_N`).

### $G(Q)$ Computation

For each $(\phi_\text{initial}, Q_\text{initial})$:

- Background evolution is solved
- Analytical power spectrum is computed
- Stochastic perturbation equations are solved and numerical power spectrum obtained after averaging
- The ratio gives $G(Q)$

---

## Module Roles

- `bg/`: Background evolution only  
- `pert/`: Background + analytical PS + perturbations  

```{note}
The analytical power spectrum is always computed inside the `pert` module.
```

---

## Smoothing

Smoothing of $G(Q)$ is performed in `GQ_Plotting_NB.ipynb` using 1-D splines with knots.

---

## Common Pitfalls

- Incorrect parameter ranges → invalid or incorrect initial conditions  
- Function signature mismatch → runtime errors or incorrect results  
- Forgetting to recompile after modifying `.cpp` files  

```{warning}
If no valid initial condition is found, the code may still write incorrect values to `ics.dat`. Users should verify the output of `find_ICs.py` by checking that the reported "function value at" is close to zero for all $Q$ values.

The function value corresponds to:

$$N_{\mathrm{end}} - \texttt{dur_N}$$

Typical acceptable values are $\sim 10^{-3}$ or smaller. Significant deviation from zero indicates that the desired inflation duration was not achieved.
```

---

## Summary

- Define your WI model in `model_calc.cpp`  
- Ensure consistent function signatures across C++ and Python  
- Use `find_ICs.py` → compute initial conditions  
- Use `find_GQ.py` → compute $G(Q)$  
- Validate and smooth results  

This module provides a flexible and powerful framework for computing $G(Q)$ for a wide range of Warm Inflation models.