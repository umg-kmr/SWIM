:math:`G(Q)` Calculator
========================

This module is subdivided into ``bg`` and ``pert`` modules.

``bg`` submodule
__________________

The main function of this module is to compute the background dynamics of a specified WI model. As a part of :math:`G(Q)` calculator this submodule is used to calculate the initial conditions for different :math:`Q` values.

The directory itself contains 3 files:

1. ``Bg.cpp`` is the program that numerically integrates the background and computes the required parameters.
2. ``model_calc.cpp`` is where the WI model is specified. It is only required to enter the model potential and its derivative if :math:`\Upsilon = C_{\Upsilon} T^p \phi^c`. Otherwise ``Ups_wo_Cy`` function needs to be modified with the form of :math:`\Upsilon` being considered. 
3. ``libbg.so`` is the compiled library that interfaces with python helper scripts in the :math:`G(Q)` calculator module.

Some examples of potentials commonly seen in literature are already available in the code itself. To modify the WI model in ``model_calc.cpp``:

1. Modify the potential function ``V`` as such:

.. code-block:: cpp

   auto V = [m] (double phi) -> double {
            return 0.5*m*m*phi*phi; //V = 1/2 m^2 \phi^2
        };

Most of the functions in SWIM are ``lambda`` functions that can take arguments as well as capture fixed parameters from the top-level function. The square brackets are where the captured parameters (fixed) specified. It should be modified accordingly to the potential specified. 

2. Modify the potential derivative function ``Vd`` as such:

.. code-block:: cpp

   auto Vd = [m] (double phi) -> double {
            return m*m*phi; //V' =  m^2 \phi
        };
        
3. Depending on the model, you can also modify the ``Ups_wo_Cy`` function. The default form is :math:`\Upsilon \propto T^p \phi^c` .

.. code-block:: cpp

   auto Ups_wo_Cy = [p,c] (double phi,double T) -> double {  //Form of Upsilon without the constant
            return pow(T,p) * pow(phi,c);
        };

The constant ``Cy`` (:math:`C_{\Upsilon}`) is internally calculated using :math:`\Upsilon = 3HQ`. 

4. Finally, also modify the signature of the top-level function ``model`` to include your model parameters with their data types.

.. code-block:: cpp
   
   void model (double phi_ini,double Q_ini, double gst, double m, int p, int c, int hybrid_inf)
   
Generally, you can keep all the parameters (``phi_ini``, ``Q_ini``, ``gst``, ``p``, ``c``, ``hybrid_inf``) and add the other model parameters specific to your use case.

Apart from these functions, no other function needs to be modified for WI model definition. Once done with the modification make sure to run the compilation script again.

``pert`` submodule
___________________

This submodule performs the perturbation calculations for WI using the stochastic formalism. It is primarily used to calculate the growth function :math:`G(Q)` using the python helper script ``find_GQ.py``. The directory structure of this submodule is similar to the ``bg`` submodule with code files having exactly the same name. ``Bg.cpp`` contains functions to integrate the stochastic differential equations and finally compute :math:`G(Q)`. ``model_calc.cpp`` is where the model has to be specified again so that the perturbation module uses the correct model for calculations.

The set-up procedure is similar to ``bg`` module:

1. Specify the model potential and its derivative in the functions ``V`` and ``Vd``.

2. Additionally, also specify the second derivative of the potential in ``Vdd``.

.. code-block:: cpp

   auto Vdd = [m] (double phi) -> double {
            return m*m; //V'' = m^2
        };
        
The second derivative of the potential is required for the perturbation equations.

3. If the model :math:`\Upsilon` differs from the default form of the code, modify the ``Ups_wo_Cy`` function similar to the ``bg`` module. 

4. Skip this step if no modification is done to :math:`\Upsilon`. Specify the partial derivatives of :math:`\Upsilon(\phi,T)` in the functions ``pT_Ups`` (temperature derivative) and ``pph_Ups`` (phi derivative).

.. code-block:: cpp

   auto pT_Ups = [Cy,p,c] (double phi, double T) -> double {
            return p * Cy * pow(T,p-1.0) * pow(phi,c);
        };

   auto pph_Ups = [Cy,p,c] (double phi, double T) -> double {
            return c * Cy * pow(T,p) * pow(phi,c-1.0);
        };
        
5. Similar to the ``bg`` module, change the signature of the top-level ``model`` function to accept the newly added model parameters.

.. code-block:: cpp

   void model (double phi_ini,double Q_ini, double gst, double m, int p, int c, int therm, int rad_noise, int hybrid_inf) 
   
Keep the parameters (``phi_ini``, ``Q_ini``, ``gst``, ``p``, ``c``, ``therm``, ``rad_noise``, ``hybrid_inf``) as it is and add your model parameters with their data types.

With this the model specification in the ``pert`` module is done and the code should be re-compiled.

Now we can move on to compute :math:`G(Q)` by using the helper scripts ``find_ICs.py`` and ``find_GQ.py``. These scripts are a convenient way to interface with the ``C++`` code and thus is the preferred way of using the module. It it possible to directly interact with the ``C++`` code by compiling the code with the inclusion of a ``main`` function in the ``model_calc.cpp`` files and removing the ``-shared`` and ``-fPIC`` compilation flags.

``find_ICs.py`` script
_______________________

The main function of this script is to generate initial values of :math:`\phi_i` such that the duration of inflation is equal to the required number of e-folds specified by ``dur_N`` parameter. The python script sets the global parameter values and repeatedly calls the compiled ``C++`` library from the ``bg`` submodule for different :math:`Q` values. Finally the computed values of :math:`\phi_i` and :math:`Q` are stored in a file which can then be used by the ``find_GQ.py`` script.


To use the helper script:

1. Specify the model parameters at the top level in the script.

   .. data:: m
      :type: double
      :value: 1e-6
      
      Potential parameter
      
   .. data:: gst
      :type: double
      :value: 106.75
      
      Relativistic degrees of freedom (:math:`g_*`)
   
2. Specify the :math:`\Upsilon` parameters.

   .. data:: p
      :type: int
      :value: 3
      
      :math:`T^p`
      
   .. data:: c
      :type: int
      :value: -2
      
      :math:`\phi^c`
      
3. If working with hybrid inflation, set these two parameters.
      
   .. data:: hybrid_inf
      :type: int
      :value: 0 or 1
      
      Parameter indicating inflation type. Can be either ``0`` (disabled) or ``1`` (enabled). 
      
   .. data:: ph_crit
      :type: double
      :value: 0.2
      
      Parameter to set the critical value of the field. To be used with ``hybrid_inf=1``.
   
4. Finally, set the other code related parameters.

   .. data:: Qlow
      :type: double
      :value: 1e-9
      
      Lower limit for the range of :math:`Q` values.  
      
   .. data:: Qup
      :type: double
      :value: 1e4
      
      Upper limit for the range of :math:`Q` values.
      
   .. data:: npts
      :type: int
      :value: 100
      
      Number of points to be generated between ``Qlow`` and ``Qup``, logarithmically spaced.
      
   .. data:: dur_N
      :type: double
      :value: 60.0
      
      Parameter to specify the duration of inflation in e-folds.
       
   .. data:: Nprocs
      :type: int
      :value: 12
      
      Number of CPUs to be used for multiprocessing. To utilize all available CPU resources, set this number to the number of logical processors (threads/core :math:`\times` cores) in your system.
      
5. Finding the value of :math:`\phi` such that inflation ends at ``dur_N`` is effectively a minimization problem and a versatile albeit inefficient ``brute_force`` algorithm is used to achieve this in this code. The algorithm finds a local minima of the objective function specified by the parameters:

   .. data:: lower_bound
      :type: double
      :value: -1.0
      
      Parameter to specify the lower bound of the input parameter (:math:`\phi`) in :math:`\log_{10}`.
      
   .. data:: upper_bound
      :type: double
      :value: 1.0
      
      Parameter to specify the upper bound of the input parameter (:math:`\phi`) in :math:`\log_{10}`.
      
 As the algorithm evaluates the objective function multiple times, it is an inefficient way to minimize the function. The algorithm used can be replaced by a more efficient algorithm better suited for your use case by modifying this line of the code:
  
 .. code-block:: python

    soln = optimize.brute(objfn,ranges=ranges, args=(i,), Ns=1000, full_output=1, disp=False, workers=Nprocs,finish=optimize.fmin)
    
6. Finally, modify the ``C++`` library definition and also the signature of the function being called from the loaded library.

   .. code-block:: python
   
      ffi.cdef("void model (double phi_ini,double Q_ini, double gst, double m,int p, int c, int hybrid_inf);void clear_Nend();extern double Nend;void set_phi_crit (double x);",override=True)
      
   The command above specifies the functions (and their signatures) that are to be loaded from the ``C++`` shared library of the ``bg`` module. You only need to modify the signature of ``void model`` function to match with the definition in ``model_calc.cpp`` of the ``bg`` module.
   
   .. code-block:: python
   
      lib.model(phi0, Q0, gst, m, p, c,hybrid_inf)
      
   This is where we pass the python parameter values to the shared library and hence should match the signature of ``void model`` exactly. The name of the python variables can be different from the ``C++`` variable names but the logical quantity held by these parameters should match. For instance, ``phi0 = phi_ini``.
   
Now the helper script is ready to compute the initial conditions. Run the command ``python find_ICs.py`` in the conda environment setup during the installation stage. The script will create a file ``ics.dat`` in the same directory with the format: ``phi_initial``, ``Q_initial``.

.. note::
 
   If the helper script is re-run, the new computation results will keep on appending to the same ``ics.dat`` file.
 
``find_GQ.py`` script
_______________________

This script computes the growth function :math:`G(Q)` utilizing the initial conditions calculated by the ``find_ICs.py`` script in the previous step. The shared ``C++`` library of the ``pert`` module is repeatedly called for each pair of initial conditions in the ``ics.dat`` file. The ``pert`` shared library will internally calculate both the numerical scalar power spectrum and analytical power spectrum and finally save the ratio of the two in a file ``GQ.dat``.

The usage of the helper script is similar to the ``find_ICs.py`` script with a few additional parameters:

1. Specify the model and :math:`\Upsilon` parameters as before in ``find_ICs.py`` script.

2. Include the parameters for hybrid inflation if required.

3. Set the other parameters:
   
   .. data:: therm
      :type: int
      :value: 0 or 1
      
      Parameter indicating the state of thermalization of inflaton with the radiation bath. ``0`` for no thermalization and ``1`` inflaton is thermalized and we assume a Bose-Einstein distribution.
      
   .. data:: rad_noise
      :type: int
      :value: 0 or 1
      
      Set to ``1`` to include the thermal noise term in the radiation perturbation equation and ``0`` to not.
      
   .. data:: Nrealz
      :type: int
      :value: 1024
      
      Parameter that sets the number of realiazations of the stochastic noise term to average over. Higher values will increase computation time.
      
      
   .. data:: Nstar
      :type: double
      :value: 7.0
      
      Parameter that sets the e-folding value at which the power spectrum should be computed. It should be set in a way such that the mode under study should have enough e-folds before to initialize and after to freeze.
      
4. Like before, modify the ``C++`` library definition and also the signature of the function being called from the loaded library.

   .. code-block:: python
   
      ffi.cdef("void model (double phi_ini,double Q_ini, double gst, double m, int p, int c, int therm, int rad_noise,int hybrid_inf);void clear_vars();void set_phi_crit (double x);extern double Nend;extern double GQ_val;extern double Q_val; void set_globals (int N_realizations, double Nstar, int verbosity);",override=True)
      
   You only need to modify the signature of ``void model`` function to match with the definition in ``model_calc.cpp`` of the ``pert`` module.
   
   .. code-block:: python
   
      lib_pert.model(i,j,gst,m,p,c,therm,rad_noise,hybrid_inf)      
      
   where ``i`` and ``j`` are dummy variables looping over the pair of :math:`\phi` and :math:`Q` values from ``ics.dat``.
   
Now the helper script can be run using the command ``python find_GQ.py`` in the conda environment setup during the installation stage. The script will create a file ``GQ.dat`` in the same directory with the format: :math:`Q_*`, :math:`G(Q_*)`.

.. note:: 
  
   1. If the helper script is re-run, the new computation results will keep on appending to the same ``GQ.dat`` file.
   
   2. The multiprocessing method in the ``find_GQ.py`` script is different from ``find_ICs.py``. In the ``find_GQ.py`` case, multiprocessing is implemented in the ``C++`` library itself using ``OpenMP`` and by default it uses all the available CPUs in the system. To change the number of CPUs available to ``OpenMP`` you can update the flag:
   
      .. code-block:: bash
   
         export OMP_NUM_THREADS=12
      
  This change is temporary and is only valid untill the current shell is closed. 
  
  
``GQ_Plotting_NB.ipynb`` notebook
__________________________________

This is ``Jupyter`` notebook that includes codes to visualize the computation of the ``GQ_Calculator``. It is optional to use this notebook and you can use your own preferred plotting facilities. However, it has an important code-block that smooths the raw output of the ``GQ_Calculator``. This is important if you wish to use the computed :math:`G(Q)` in the semi-analytical approach to compute the WI power-spectrum. ``SWIM`` provides another submodule ``SA_PS_Calculator`` that implements this approach, can read an external ``GQ.dat`` file and perform statistical analysis using ``Cobaya`` and recent cosmological likelihoods. It is important to smooth the raw output before using it in the ``SA_PS_Calculator`` module.

The usage of the ``Jupyter`` notebook is mostly self-explanatory,

1. The notebook loads the raw output ``GQ.dat`` and then plots it.

2. In the next step smoothing is done using ``UnivariateSpline`` in log scale and then plots both the raw and smooth output to compare. Here, the parameter `s` can be tuned to achieve varying degree of smoothness. For more information on ``s`` and other paramters of ``UnivariateSpline`` visit `here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html>`_.

3. The notebook also provides a code to compare the derivatives of raw and smooth :math:`G(Q)`. It can be used as a visual tool to check the smoothing effect of the previous step.

4. Once you are satisfied with the smoothing of the raw data, the final code block can be run to save the smoothed data to a file ``GQ_smooth.dat``. This file can then be used by the ``SA_PS_Calculator``.

Optional functionality
_______________________

In the ``Bg.cpp`` file of the ``pert`` module a commented out version of the analytical power spectrum can be found that doesn't drop the slow-roll parameters :math:`\eta_V` and :math:`\beta_V`. The implementation in the code directly follows from `R.O.Ramos and L.A. da Silva <https://inspirehep.net/literature/1219309>`_.

To use this form of the analytical power spectrum just comment out the other implementation and uncomment:

.. code-block:: cpp
  
   auto Pp = [Nask,phiasN,phpasN,TasN,H,Q,therm,V,Vd,Vdd,Ups,pph_Ups] (double k) -> double {
   
   
