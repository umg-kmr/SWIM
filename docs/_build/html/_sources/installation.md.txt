# Installation

```{note}
SWIM has currently been tested only on Linux-based systems.
```

## Clone Repository

Clone the SWIM repository to a convenient location:

```bash
git clone https://github.com/umg-kmr/SWIM.git
cd SWIM
```

---

## Python Environment Setup

We recommend using a dedicated `conda` environment.

```bash
conda create -n SWIM python=3.13
conda activate SWIM
```

---

## Dependencies

Install the required Python and compiler dependencies:

```bash
conda install -c conda-forge gcc gxx gfortran numpy matplotlib scipy scikit-learn jupyterlab cffi joblib
```

---

## Boost C++ Libraries

SWIM relies on Boost libraries for numerical routines (e.g., ODE solvers).

Download the latest version from:  
https://www.boost.org/releases/latest/

Extract the archive:

```bash
tar -xvf boost_*.tar.bz2
```

```{important}
Downloading Boost is required for SWIM to function.
```

```{note}
SWIM has been tested with Boost v1_87_0, but should work with later versions.
```

---

## Multiprocessing Support

SWIM supports parallel computation via OpenMP and MPI.

Install:

```bash
conda install -c conda-forge openmp openmpi
```

---

## Cobaya (Optional)

SWIM can interface with `Cobaya` for Bayesian parameter inference.

Installation guide:  
https://cobaya.readthedocs.io/en/latest/installation.html

```{note}
Ensure all Cobaya dependencies are installed within the same SWIM environment.
```

---

## Compilation

No installation is required beyond compilation.

### Step 1: Configure paths

Edit `compile_SWIM.sh` and set:

```bash
export BOOST_PATH=/path/to/boost
export SWIM_PATH=/path/to/SWIM
```

---

### Step 2: Compile

```bash
chmod +x compile_SWIM.sh
./compile_SWIM.sh
```

This will compile all required modules with optimized flags.

---

## Manual Compilation (Optional)

Alternatively, compile modules manually:

```bash
g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/GQ_Calculator/bg/libbg.so <path-to-SWIM>/GQ_Calculator/bg/model_calc.cpp -lm -fopenmp

g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/GQ_Calculator/pert/libpert.so <path-to-SWIM>/GQ_Calculator/pert/model_calc.cpp -lm -fopenmp

g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/PS_Calculator/libmodel.so <path-to-SWIM>/PS_Calculator/model_calc.cpp -lm -fopenmp

g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/SA_PS_Calculator/libmodel.so <path-to-SWIM>/SA_PS_Calculator/model_calc.cpp -lm -fopenmp
```
```{note}
The compilation flags provided above correspond to those used during testing and are optimized for performance. Users may modify these flags depending on their system architecture and their desired balance between performance and numerical accuracy.
```

---
## Ready to Use

SWIM is now ready to be used.