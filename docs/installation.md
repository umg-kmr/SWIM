# Installation

```{note}
SWIM is supported natively on Linux-based systems, with macOS and Windows support is provided through Docker.
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

### For Mac/Windows users

Install Docker from the official website: https://www.docker.com/ and pull the latest image of SWIM:

```bash
docker pull umgkmr/swim:latest
```

Run the container in interactive mode and expose port 8888 for JupyterLab:

```bash
docker run -it -p 8888:8888 --name swim_container umgkmr/swim:latest
```

The image includes all required dependencies to run and compile SWIM, including Cobaya. To keep the image size manageable, external likelihoods and cosmological codes are not preinstalled. These can be installed following the Cobaya documentation: https://cobaya.readthedocs.io/en/latest/installation_cosmo.html

Once the container starts, you are placed in a bash shell as the default user `swim`. From this shell, SWIM can be compiled using:

```bash
./compile_SWIM.sh
```

To run JupyterLab:

```bash
jupyter-lab --port=8888 --ip=0.0.0.0 --no-browser --ServerApp.token=''
```

Then open in your browser: `http://localhost:8888`

The following commands should be run on your host system (not inside the container).
To copy files from the container to your host system:

```bash
docker cp swim_container:/path/in/container /path/on/host
```

To restart the container:

```bash
docker start -ai swim_container
```

```{note} 
GUI applications from Cobaya (e.g. `getdist-gui`) are not supported within the Docker container.
```

A Dockerfile is available in the `docker/` directory and can be modified to suit individual use cases.