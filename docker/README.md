# SWIM Docker Environment

This Docker image provides a ready-to-use environment for running and developing **SWIM (Stochastic Warm Inflation Module)**.

---

## Contents of the Image

The container includes:

- **Base system**
  - Debian (slim) with essential development tools
  - `gcc`, `g++`, `gfortran`, `cmake`, `make`
  - `openmpi` for parallel computing
  - `libopenblas` 

- **Utilities**
  - `git`, `curl`, `wget`
  - `vim`, `nano`, `htop`
  - `sudo` (passwordless for user `swim`)

- **Python environment**
  - Python 3.14
  - Scientific stack:
    - `numpy`, `scipy`, `matplotlib`
    - `scikit-learn`, `joblib`, `numba`
  - Cosmology / inference:
    - `cobaya`
    - `candl-like` for SPT datasets
  - Other dependencies:
    - `cffi`, `cython`, `pyside6`, `jupyterlab`, `mpi4py`

- **C++ dependencies**
  - Boost (v1.91.0, locally installed at `/home/swim/boost`)

- **SWIM source code**
  - Cloned into:
    ```
    /home/swim/SWIM
    ```

---