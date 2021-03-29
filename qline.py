
from qutiepy import *
import numpy as np
from scipy.linalg import expm

"""

Ax = b

"""


# ????????

t = 4   # bits in b
A = np.eye(t)
b = register(4)  
T = 2 ** t    # states in b = 16
amps = np.sqrt(2/T) * np.array([np.sin((np.pi*(tau+0.5)/T)) for tau in range(T)])
phi0 = register(t)
phi0.setAmps(amps)
#print(len(phi0.amps))  ## 16
print(b)
input()
phi0b = prod(phi0, b)
#print(phi0b.NStates)  ## 256??? should be tT???
#print(len(np.kron(phi0.amps, b.amps)))
#input()
t0 = 1

hamMatTerms = []

for tau in range(T):                    #construct hamilton operator
    tautau = np.zeros((T, T))
    tautau[tau, tau] = 1             # t x t
    print("tautau: ", tautau.shape)
    oper = expm(1j*tau*t0*A/T)      # t x t
    print("oper: \n", oper.shape)
    term = np.kron(tautau, oper)    # tT x tT
    print("TERM: ", term.shape)
    hamMatTerms.append(term)

print()
hamMat = np.sum(hamMatTerms, axis=0)   
print(hamMat.shape)                    # these should match...
print(phi0b.NStates)                   # ... but they dont
ham = genericGate(phi0b.NBits)        #make it a gate
ham.matrix = hamMat

