
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

A = np.array([[1,-0.333],[-0.333,1]])
bBits = 1
b = register(bBits)
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
    print("oper: ", oper.shape)
    term = np.kron(tautau, oper)    # tT x tT
    print("TERM: ", term.shape)
    hamMatTerms.append(term)

print("==============================")

hamMat = np.sum(hamMatTerms, axis=0)                   
ham = genericGate(bBits+t)        #make it a gate
ham.matrix = hamMat
print(hamMat.shape)

QFT = genericGate(t)    ###### only phi gets qft'd
QFTMat = QFTmatrix((2**bBits)*T, 1j)
print(QFTMat.shape)
QFT.matrix = np.kron(QFTMat, np.eye(bBits))

phi1b = QFT(ham(phi0b))
print(phi1b)

ancilla = register(1)    
phi1ba = prod(phi1b, ancilla)

#first ? bits are lambdaK, next ? bits are Uj, finally ancilla.

### [C C C ... C C 0 0 ... 0 0 Rx(?)]

toKron = [np.eye(2)] * t + [np.zeros((2,2))] * bBits + [[0,-1],[1,0]] #?????

res = toKron[0]
for m in toKron[1:]:
    res = np.kron(res, m)

contRotGate = genericGate(t + bBits + 1)
contRotGate.matrix = res

phi1ba = contRotGate(phi1ba)

iQFT = genericGate(t)
iQFT.matrix = np.array(np.asmatrix(QFT.matrix).H)





















