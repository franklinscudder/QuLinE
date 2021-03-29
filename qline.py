
from qutiepy import *
import numpy as np
from scipy.linalg import expm

"""

Ax = b

"""

def QFTmatrix(N, omg):
    matrix = np.zeros((N,N), dtype=complex)
    for x in range(N):
        for y in range(N):
            matrix[x,y] = omg ** (x*y)
    
    return matrix / np.sqrt(N)

bAmps = [3, 4]

A = np.eye(len(bAmps))
b = register(1)
b.setAmps(bAmps)
t = 4  # bits in phi
T = 2 ** t    # states in phi
amps = np.sqrt(2/T) * np.array([np.sin((np.pi*(tau+0.5)/T)) for tau in range(T)])
phi0 = register(t)
phi0.setAmps(amps)
phi0b = prod(phi0, b)

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


hamMat = np.sum(hamMatTerms, axis=0)                   
ham = genericGate(phi0b.NBits)        #make it a gate
ham.matrix = hamMat

QFT = genericGate(phi0b.NBits)    # only phi gets qft'd
QFT.matrix = QFTmatrix(phi0b.NStates, 1j)

phib = QFT(ham(phi0b))

ancilla = register(1)    
UinvBinit = prod(phib, ancilla)




















