# SWIM: Stochastic Warm Inflation Module

![Platform](https://img.shields.io/badge/platform-Linux-blue)
![Python](https://img.shields.io/badge/python-3.10+-blue)
![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b)

SWIM is a numerical framework for computing the power spectrum in Warm Inflation (WI) models and performing parameter inference using `Cobaya`. The code provides a flexible and modular pipeline for both semi-analytical and fully numerical approaches.

The main components of SWIM include:

- **G(Q) correction function calculator**
- **Semi-analytical power spectrum module**
- **Full numerical stochastic solver**
- **Random Forest emulator for accelerated inference**

---

# Installation Instructions:
First, clone the repository to a convenient location on your system:

```bash
git clone https://github.com/umg-kmr/SWIM.git
cd SWIM
```

We recommend creating a dedicated `conda` environment to manage dependencies. After installing Anaconda or Miniconda, create a new environment:

```bash
conda create -n SWIM python=3.13
```

Then activate the environment:

```bash
conda activate SWIM
```

### Dependencies

Install the required packages for SWIM using `conda`:

```bash
conda install -c conda-forge gcc gxx gfortran numpy matplotlib scipy scikit-learn jupyterlab cffi joblib
```

SWIM relies on the Boost C++ libraries for numerical routines. Download the latest release from:

https://www.boost.org/releases/latest/

and extract it to a convenient location using:

```bash
tar -xvf boost_*.tar.bz2
```

After extraction, ensure that the Boost library paths are correctly specified in the `compile_SWIM.sh` script.

### Multiprocessing support

SWIM supports parallel computation using OpenMP and MPI (useful for Cobaya).

Install the required packages using `conda`:

```bash
conda install -c conda-forge openmp openmpi 
```

### Install Cobaya (Optional)

SWIM can interface with `Cobaya` for Bayesian analysis and parameter inference. While not required for basic functionality, installing Cobaya is recommended for performing parameter estimation.

Follow the official installation guide:

https://cobaya.readthedocs.io/en/latest/installation.html

All dependencies required by Cobaya should be installed within the same `conda` environment created for SWIM.

### Compilation

After installing all the required dependencies, SWIM can be compiled using the `compile_SWIM.sh` script.

First, edit the paths to the SWIM directory and Boost libraries in the script:

```bash
export BOOST_PATH=/path/to/boost
export SWIM_PATH=/path/to/SWIM
```

Then make the script executable and run it:

```bash
chmod +x compile_SWIM.sh
./compile_SWIM.sh
```

SWIM is now ready to be used.

---

# Quick Start Guide

As an illustrative example, the code includes an implementation of the following Warm Inflation (WI) model:

$$
V(\phi) = \dfrac{1}{4} V_0 \phi^4
$$

$$
\Upsilon(\phi,T) = C_{\Upsilon} T^3
$$

### $G(Q)$ Calculator

To compute the $G(Q)$ correction factor, use this module.

Navigate to the directory:

```bash
cd GQ_Calculator
```

Remove any existing data files to generate fresh outputs:

```bash
rm ics.dat GQ.dat GQ_smooth.dat
```

Run the Python script to compute the initial conditions:

```bash
python -u find_ICs.py
```

This step may take some time depending on system specifications. Upon completion, a file named `ics.dat` will be generated.

---

Next, run the script to compute the $G(Q)$ correction factor:

```bash
python -u find_GQ.py
```

After completion, the output will be saved in the file `GQ.dat`. The computed $G(Q)$ can be visualized using the Jupyter notebook:`GQ_Plotting_NB.ipynb`

The notebook also includes utilities to smooth the raw output and save it as `GQ_smooth.dat`.

### Semi-Analytical Power Spectrum Calculator

This module is used to perform parameter inference for WI models using the semi-analytical power spectrum, which incorporates the $G(Q)$ correction factor.

To use this module, `Cobaya` must be installed along with the required likelihoods (Planck, ACT DR6, DESI DR2, and SPT). Follow the official Cobaya guide to install external likelihoods:

https://cobaya.readthedocs.io/en/latest/likelihood_external.html#list-of-external-packages

---

From the main SWIM directory, navigate to:

```bash
cd SA_PS_Calculator
```

Remove any existing chains:

```bash
rm -rf chains_Asns chains_fullCMB
```

To perform parameter inference using constraints only on $A_s$, $n_s$, and $r$, run:

```bash
mpirun -n 8 cobaya-run Input_asns.yaml
```

This launches 8 parallel chains.

For a single chain, run:

```bash
cobaya-run Input_asns.yaml
```
To perform parameter inference including full CMB power spectrum constraints (requires `CAMB`), run:

```bash
mpirun -n 8 cobaya-run Input.yaml
```

By default, CAMB uses all available CPU threads. When running multiple chains in parallel, it is recommended to limit the number of threads per chain.

You can check the total number of available CPU threads using:

```bash
nproc --all
```

Set the number of OpenMP threads as:

```bash
export OMP_NUM_THREADS=<total_threads / number_of_chains>
```

For example, if your system has 32 threads and you run 8 chains:

```bash
export OMP_NUM_THREADS=4
```
The resulting chains can be analyzed using the `getdist` utility provided by Cobaya.

### Numerical Power Spectrum Calculator

This module computes the Warm Inflation (WI) primordial power spectrum $P_{\mathcal{R}}(k)$ numerically, without relying on the $G(Q)$ correction factor. It operates independently of the `GQ_Calculator` and `SA_PS_Calculator` modules.

The module can be used to:
- compute the full WI power spectrum for a given set of model parameters, and  
- perform parameter inference using the numerical power spectrum by interfacing with `Cobaya`.

From the main SWIM directory, navigate to:

```bash
cd PS_Calculator
```
To compute the full WI power spectrum for the model:

$$
V(\phi) = \dfrac{1}{4} V_0 \phi^4
$$

$$
\Upsilon(\phi,T) = C_{\Upsilon} T^3
$$

The WI model and corresponding parameter values are already implemented in the code.

Navigate to the power spectrum directory and remove existing data files:

```bash
cd Power_Spectrum
rm bg.dat ps.dat PT_kp.dat
```

(Optional) If OpenMP threads were previously limited, reset them to use all available CPUs:

```bash
export OMP_NUM_THREADS=$(nproc --all)
```

Then run the Python script:

```bash
python -u ps_script.py
```

The script generates the following files:

- `bg.dat` — background evolution of $\phi$, $\phi'$, and $T$ as a function of e-folds $N$  
- `ps.dat` — numerical WI power spectrum as a function of $k$  
- `PT_kp.dat` — tensor power spectrum amplitude at the pivot scale

The notebook `Plotting_NB.ipynb` can then be used to:

- visualize the background evolution  
- plot the raw power spectrum $P_{\mathcal{R}}(k)$  
- compute derived quantities such as $A_s$, $n_s$, and $r$

---

To perform parameter inference using `Cobaya` (requires Cobaya to be installed), navigate to:

```bash
cd ../Emulator/RF_Acc_Cobaya
```

Remove any existing chains:

```bash
rm -rf chains
```

Then run Cobaya:

```bash
cobaya-run Input_asns.yaml
```

This runs a single MCMC chain and uses a Random Forest Regression (RFR) emulator to accelerate the inference.

During the initial phase, the code evaluates the full numerical solver and stores valid samples. Once a sufficient number of points (100) have been accumulated, the emulator is trained and saved as: `rf_model.pkl`

This trained model can then be used to perform parameter inference more efficiently.

#### Notes

- This module is the most computationally expensive component of SWIM and is best suited for execution on a high-performance computing (HPC) system with multiple CPU threads.  
- After training, the emulator provides a fast approximation to the full numerical solver in the high-likelihood region of parameter space.

---

## Platform Compatibility

SWIM has currently been tested only on Linux-based systems. Compatibility with macOS or Windows has not been extensively verified.

---

## Repository Structure

```
SWIM/
├── GQ_Calculator/
│   ├── bg/                     # Background evolution
│   ├── pert/                   # Perturbation module 
│   ├── find_ICs.py            # Initial condition solver
│   ├── find_GQ.py             # Computes G(Q)
│   └── GQ_Plotting_NB.ipynb   # Visualization & smoothing

├── SA_PS_Calculator/
│   ├── model_calc.cpp         # WI model definition
│   ├── Bg.cpp                 # Solver
│   ├── llihood.py             # Likelihood definition (A_s, n_s, r constraints)
│   ├── An_CAMB.py             # CAMB external primordial power spectrum
│   ├── Input_asns.yaml        # Cobaya config (A_s, n_s, r)
│   └── Input.yaml             # Full CMB inference

├── PS_Calculator/
│   ├── Emulator/              # RF emulator based inference
│   ├── Power_Spectrum/        # Numerical power spectrum calculation
│   ├── model_calc.cpp         # Model definition
│   └── Bg.cpp                 # Solver

├── compile_SWIM.sh            # Compilation script
├── arXiv_plots_data/          # Figures and data from paper
└── README.md
```

## Contact

For questions, feedback, or issues related to SWIM, please contact:

- **Umang Kumar**  
  Email: umang.kumar_phd21@ashoka.edu.in  

Alternatively, feel free to open an issue on this repository.