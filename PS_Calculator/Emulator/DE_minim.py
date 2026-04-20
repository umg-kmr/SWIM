import numpy as np
from scipy.optimize import differential_evolution
import llihood as L


def loglk (x):
    x = np.asarray(x)
    phi0,gst,V0,Cy = 10**x
    p = int(3)
    c = int(0)
    f = 5.0
    therm = int(0)
    rad_noise = int(0)
    Np = 1.0
    return L.logp(phi0,gst,V0,f,Cy,Np,p,c,therm,rad_noise)
def print_status(xk, convergence):
    """
    xk: The best parameter set of the current generation
    convergence: The fractional convergence (0 to 1)
    """
    # 1. Convert log parameters back to linear for readability
    #    (Wrap in np.array just in case xk comes as a list)
    # real_params = np.array(xk)
    
    # 2. Print the details
    print(f"--- Generation Best ---")
    print(f"Log Params:  {xk}")
    
    # Optional: Print the score if you want to track it
    # Note: This adds a small overhead of 1 extra calculation per step
    # score = objective_function(xk)
    # print(f"Score:       {score}")
    
    # print("-" * 50)

bounds = [(-7.0,np.log10(np.pi*5.0)),(0.0,4.0),(-80.0,-1.0),(-10.0,50.0)] #phi_ini, gst, V0, Cy

result = differential_evolution(
    loglk, 
    bounds, 
    maxiter=1000,        # Maximum number of generations
    popsize=100,          # Population size multiplier
    disp=True,
    workers=1,
    #x0=x0,
    callback=print_status,
)

# 3. Display the results
print("--- Optimization Results ---")
print(f"Success: {result.success}")
print(f"Message: {result.message}")
print(f"Optimal Parameters: {result.x}")
print(f"Objective Value at Minimum: {result.fun}")
print(f"Number of Evaluations: {result.nfev}")
