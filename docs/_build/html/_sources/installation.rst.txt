Installation
============

.. note::

   Currently we have only tested the module on linux systems.

Recommended python environment setup
_________________________________________

1. Create a fresh python environment using conda:

   ::
   
      conda create -n SWIM
      
2. Activate the environment and install the following packages:
   NumPy, SciPy, CFFI, matplotlib, JupyterLab, numba (optional).
   
   ::
      
      conda activate SWIM
      conda install -c <channel-name (eg. conda-forge)> <package-name (eg. numpy)>
      
   Additionally, if you wish to use the Artificial Neural Network (ANN) emulator for SWIM, please also install PyTorch. Refer to the pytorch `documentation <https://pytorch.org/get-started/locally/>`_ and install for your specific system.  
    
3. It is also recommended to have a working installation of open-mpi (for parallelization) and C++ compiler (g++ for linux) if not already installed:

  ::
   
     conda install -c conda-forge openmpi
     conda install -c conda-forge gcc
     conda install -c conda-forge gxx (or cxx-compiler)
     
4. The C++ implementation of the code uses the boost libraries to implement for example ODE solvers. No need to install it, just download (and extract) the package from `boost.org <https://www.boost.org/releases/latest/>`_ at a convenient location which can then be pointed to the compiler during compilation. **Download of this library is crucial for the code to work**.

.. note::

   Currently only tested with boost v1_87_0. But should work with later versions as well.

Compilation
____________

No installation is required to use SWIM. To obtain and use SWIM follow the following steps:

1. Download the package from GitHub:

   ::
   
      git clone https://github.com/umg-kmr/SWIM
      
2. Enter the cloned repository and edit the file: "compile_SWIM.sh". You need to change the following lines:

   ::
   
     export $BOOST_PATH=<change-to-boost-location>
     export $SWIM_PATH=<change-to-SWIM-package-location>
     
3. Make the script executable and then execute:
   
   ::
     
      chmod +x compile_SWIM.sh   
      ./compile_SWIM.sh
      
   The script will automatically compile the required modules with the compilation flags specified (these flags were used for testing).  
   
4. Alternate method: It is possible to compile the modules without using the convenience script by executing the following commands:

   ::
   
      g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/GQ_Calculator/bg/libbg.so <path-to-SWIM>/GQ_Calculator/bg/model_calc.cpp -lm -fopenmp
      g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/GQ_Calculator/pert/libpert.so <path-to-SWIM>/GQ_Calculator/bg/model_calc.cpp -lm -fopenmp
      g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/PS_Calculator/libmodel.so <path-to-SWIM>/PS_Calculator/model_calc.cpp -lm -fopenmp
      g++ -shared -fPIC -I <path-to-boost> -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o <path-to-SWIM>/SA_PS_Calculator/libmodel.so <path-to-SWIM>/SA_PS_Calculator/model_calc.cpp -lm -fopenmp
  
      
Now SWIM is ready to be used!
