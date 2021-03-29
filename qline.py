
from qutiepy import *
import numpy as np
from scipy.linalg import expm

"""

Ax = b

"""

D = 2  #equation dimension
A = np.eye(D) # ????????

t = 4   # bits in b

b = register(t)  # 2 bits per value in b, value of B is now [8, 0]?
T = 2 ** t    # states in b
amps = np.sqrt(2/T) * np.array([np.sin((np.pi*(tau+0.5)/T)) for tau in range(T)])
phi0 = register(t)
phi0.setAmps(amps)



phi0b = prod(phi0, b)

t0 = 1
hamMat = []

for tau in range(T):
    tau_ = np.zeros((T, 1))
    tau_[:,0] = np.array([int(i == tau) for i in range(T)])
    print("tau_: ", tau_.shape)
    tau_T = tau_.T
    print("tau_T: ", tau_T.shape)
    tautau = np.kron(tau_, tau_T)
    print("tautau: ", tautau.shape)
    oper = expm(1j*tau*t0*A/T)
    print("oper: \n", oper)
    term = np.kron(tautau, oper)
    print("TERM: ", term.shape)
    hamMat.append(term)

print(np.array(hamMat).shape)
hamMat = np.sum(hamMat, axis=0) 
print(hamMat.shape)
print(phi0b.NStates)

