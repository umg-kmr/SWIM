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

