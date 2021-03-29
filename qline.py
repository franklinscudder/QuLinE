
from qutiepy import *
import numpy as np
from scipy.linalg import expm

"""

Ax = b

"""


A = np.kron(np.array([[1,0],[0,1]]),np.array([[1+0j, 1j], [-1j, 0j]])) # ????????

t = 4   # bits in b
b = register(t)  
T = 2 ** t    # states in b
amps = np.sqrt(2/T) * np.array([np.sin((np.pi*(tau+0.5)/T)) for tau in range(T)])
phi0 = register(t)
phi0.setAmps(amps)

phi0b = prod(phi0, b)
t0 = 1

hamMatTerms = []

for tau in range(T):                    #construct hamilton operator
    tautau = np.zeros((T, T))
    tautau[tau, tau] = 1
    print("tautau: ", tautau.shape)
    oper = expm(1j*tau*t0*A/T)
    print("oper: \n", oper)
    term = np.kron(tautau, oper)
    print("TERM: ", term.shape)
    hamMatTerms.append(term)

print()
hamMat = np.sum(hamMatTerms, axis=0)   
print(hamMat.shape)                    # these should match...
print(phi0b.NStates)                   # ... but they dont
ham = genericGate(phi0b.NBits)        #make it a gate
ham.matrix = hamMat

