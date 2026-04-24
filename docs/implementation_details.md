# Implementation Details

SWIM is designed with a hybrid architecture that combines the computational efficiency of C++ with the flexibility and ease of use of Python.

The computationally intensive components—such as background evolution and stochastic perturbation dynamics—are implemented in C++, while Python provides a lightweight interface for running pipelines, managing data, and performing analysis.

---

## Architecture Overview

- **C++ backend**: Handles all heavy numerical computations  
- **Python interface**: Provides user-facing scripts, plotting, and inference tools  
- **Interfacing**: Python interacts with compiled C++ shared libraries via `cffi`  

Each module compiles into shared libraries (`.so` files), which are dynamically loaded and executed from Python.

---

## C++ Backend

### ODE Solver (Background Evolution)

The background dynamics are solved using the Boost ODEInt library with an adaptive step-size integrator:

- **Method**: Runge–Kutta–Fehlberg 7(8) (`runge_kutta_fehlberg78`)  
- **Variables evolved**: $(\phi, \phi', T)$  
- **Time variable**: Number of e-folds $N$ (with $' \equiv d/dN$)  

All quantities are expressed in reduced Planck units ($M_{\mathrm{Pl}} = 1$).

The end of different phases is determined using:

- **End of inflation**: $\epsilon_H = 1$  
- **End of reheating / start of radiation domination**: $\epsilon_H = 2$  

---

### SDE Solver (Perturbations)

The stochastic evolution of perturbations is implemented using a custom scheme:

- **Method**: RI1W1 scheme (Rößler) for Itô stochastic differential equations  
- **Step-size**: Controlled via an internal empirical prescription (neither fixed nor fully adaptive)  
- **Noise**: Gaussian with zero mean, scaled with the e-fold step  

The curvature perturbation is obtained by averaging over multiple stochastic realizations.

---

### Parallelization

Parallelization is implemented using **OpenMP**:

- Independent stochastic realizations are computed in parallel  
- The final power spectrum is obtained by averaging over these realizations  

For parameter inference, parallelization across chains is handled externally via `Cobaya` using MPI.

---

### Boost Libraries

SWIM makes extensive use of Boost for numerical routines:

- ODE integration (`boost::numeric::odeint`)  
- Root finding  (`boost::math::tools`)
- Interpolation  (`boost::math::interpolators`)

---

## Python Interface

Python provides a lightweight interface to the C++ backend:

- Loads shared libraries using `cffi`  
- Passes model parameters to the C++ solvers  
- Handles execution pipelines  
- Performs plotting and data analysis  
- Implements the Random Forest Regression (RFR) emulator  

No computationally intensive tasks are performed in Python. All heavy calculations are delegated to the C++ layer.

---

## Data Flow

The typical workflow is:

```
Python → C++ (computation) → output (.dat files) → Python (analysis/plotting)
```

- C++ solvers write outputs (e.g., background evolution, power spectrum) to `.dat` files  
- Python scripts and notebooks read these files for visualization and post-processing  

---
```{admonition} Design Philosophy
:class: tip

SWIM is designed to balance computational efficiency with usability:

- C++ is used for performance-critical numerical computations  
- Python provides a flexible interface for workflow management and analysis  
- Modular design allows independent use of different components  
```
