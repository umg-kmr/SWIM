# SWIM: Stochastic Warm Inflation Module

[![Docs](https://img.shields.io/readthedocs/swim?label=docs)](https://swim.readthedocs.io/en/latest/)
[![License](https://img.shields.io/badge/license-GPL--3.0-blue)](https://github.com/umg-kmr/SWIM/blob/main/LICENSE)
![Platform](https://img.shields.io/badge/platform-Linux-lightgrey)
![Language](https://img.shields.io/badge/language-C++%20%26%20Python-blue)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-black?logo=github)](https://github.com/umg-kmr/SWIM) 
[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)]()

SWIM is a numerical framework for computing the primordial power spectrum in warm inflation models.  
It combines a full stochastic solver, semi-analytical approximations, and fast inference tools within a modular C++/Python architecture.

SWIM is intended for both detailed numerical studies of warm inflation dynamics and efficient parameter inference.

---

## Key Features

- Full stochastic warm inflation solver (background + perturbations)
- Semi-analytical power spectrum using the $G(Q)$ correction function
- Random Forest emulator for accelerated parameter inference
- Integration with Cobaya and CAMB for cosmological analysis
- Parallelized implementation (OpenMP for realizations, MPI via Cobaya)

---

## Module Overview

- **GQ Calculator**  
  Computes the correction function $G(Q)$ for semi-analytical warm inflation power spectra

- **SA_PS Calculator**  
  Fast parameter inference using semi-analytical power spectra with precomputed $G(Q)$, integrated with Cobaya and CAMB

- **PS Calculator**  
  Full stochastic numerical solver for the power spectrum, with optional Random Forest emulator for accelerated parameter inference

---

## Quick Workflow

A typical workflow in SWIM consists of:

1. Define the model in `model_calc.cpp`

2. Choose the computation submodule:

   - **Full numerical (PS Calculator):**  
     Directly compute the numerical power spectrum, with optional parameter inference using the Random Forest emulator

   - **Semi-analytical (SA_PS Calculator):**  
     Compute the correction function $G(Q)$ using the GQ calculator, then perform parameter inference using the semi-analytical module

See {doc}`quick_start` for a minimal working example.

---


## Design and Implementation

SWIM follows a hybrid design:

- **C++ backend** for numerical evolution (ODE/SDE solvers, parallel execution)
- **Python interface** for workflow control, inference, and analysis 

The framework is designed to be:

- **modular**, with separate components for background evolution, perturbations, and inference, while abstracting heavy numerical computation from the user  
- **extensible**, allowing custom inflationary potentials $V(\phi)$ and dissipative coefficients $\Upsilon(\phi, T)$  
- **efficient**, with computationally intensive routines implemented in C++ with parallelization, and optional acceleration through a Random Forest emulator
- **accessible**, with most interactions handled through Python scripts and only limited modification of C++ required for model definition  

---

## Performance

SWIM uses parallel computation (OpenMP) and supports emulation-based acceleration.  
In practice, this can provide improved performance compared to existing warm inflation solvers, particularly for inference workflows.

---

## Tested Environment

The following configuration was used for testing:

- **CPU:** AMD Ryzen 9 7900X (12 cores / 24 threads, 4.7 GHz base)
- **Memory:** 32 GB DDR5 (6000 MT/s)
- **OS:** Fedora Linux 37  
- **Kernel:** 6.5.12  

**Software stack:**

- Python 3.13  
- g++ 14.2.0  
- Boost 1.87.0  
- Cobaya 3.5.7  

---

## Citation

If you use SWIM in your work, please cite:

> *[Reference to SWIM paper — to be added]*

---

## Acknowledgments

Thanks to Gabriele Montefalcone, Alejandro Perez Rodriguez, Dipankar Bhattacharya for useful discussions and support during the development of SWIM.

Thanks also to many others who contributed through testing, feedback, and discussions at various stages of development.

## Contents

```{toctree}
:maxdepth: 1
:caption: Documentation

installation
quick_start
general_structure
implementation_details
GQ_calculator
SA_PS_calculator
PS_calculator