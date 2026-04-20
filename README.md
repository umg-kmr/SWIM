# SWIM: Stochastic Warm Inflation Module
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

## Quick Start Guide

As an illustrative example, the code includes an implementation of the following Warm Inflation (WI) model:

$$
V(\phi) = \dfrac{1}{4} V_0 \phi^4
$$

$$
\Upsilon(\phi,T) = C_{\Upsilon} T^3
$$

