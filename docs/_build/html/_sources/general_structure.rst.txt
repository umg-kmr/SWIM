General Structure of SWIM
=========================

SWIM has been divided into following sub-modules:

1. :math:`G(Q)` calculator 
2. Numerical power-spectrum calculator
3. Semi-analytical power-spectrum calculator

:math:`G(Q)` Calculator
________________________
The function :math:`G(Q)` accounts for the enhancement in the primordial power spectrum of Warm Inflation (WI) arising due to the coupling between inflaton and radiation fluctuations. This function has to be computed numerically by integrating WI fluctuation equations which is what this submodule does, similar to codes like `WI2Easy <https://github.com/RudneiRamos/WI2easy>`_ and `WarmSpy <https://github.com/GabrieleMonte/WarmSPy>`_.  Once obtained, it can be used to calculate the semi-analytical primordial power-spectrum :math:`P_{\mathcal{R}}(k)` for WI.

.. math::

   P_{\mathcal{R}}(k) = \left(\dfrac{H^2}{2\pi \dot{\phi}}\right)^2\left[\dfrac{T}{H}\dfrac{2\sqrt{3}\pi Q}{\sqrt{3+4\pi Q}} + 1 + 2\mathcal{N}\right] G(Q) \ ,
   
where all the quantities are evaluated at horizon crossing :math:`(k=aH)` and :math:`\mathcal{N}` denotes the thermal distribution of the inflaton field due to the presence of radiation.

This submodule contains all the functionality required to compute the growth function :math:`G(Q)` for a given WI model, further divided into: ``bg`` and ``pert`` modules. It also includes python scripts and a Jupyter Notebook to interact with the numerical solver.

Numerical Power Spectrum Calculator
___________________________________

Labelled as ``PS_Calculator`` in the code. It is the heart of SWIM as it directly computes the WI power spectrum numerically using a Stochastic-Langevin approach (similar to ``WarmSpy``). This submodule can be run independently of all the other submodules within SWIM, as the power spectrum is directly calculated as a function of the comoving wavenumber :math:`k`.  



