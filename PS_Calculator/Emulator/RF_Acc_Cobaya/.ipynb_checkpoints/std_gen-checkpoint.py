import numpy as np
import llihood_Observables as L

x = [23.6895,538.664,0.311553,2.05237e-15,1.0,3,0,0,0] #should match the function signature

ln10As = []
ns = []
alphs = []
betas = []
for i in range(50):
    z = L.true_solver_obs(x)
    print(z)
    ln10As.append(z[0])
    ns.append(z[1])
    alphs.append(z[2])
    betas.append(z[3])

ln10As = np.array(ln10As)
ns = np.array(ns)
alphs = np.array(alphs)
betas = np.array(betas)

print("ln10As std: ", np.std(ln10As))
print("ns std: ", np.std(ns))
print("alphs std: ", np.std(alphs))
print("betas std: ", np.std(betas))