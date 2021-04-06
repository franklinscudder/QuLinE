
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

### HAMILTONIAN SIMULATION

hamMatTerms = []

for tau in range(T):                    #construct hamilton operator
    tautau = np.zeros((T, T))
    tautau[tau, tau] = 1             # t x t
    oper = expm(1j*tau*t0*A/T)      # t x t
    term = np.kron(tautau, oper)    # tT x tT
    hamMatTerms.append(term)

hamMat = np.sum(hamMatTerms, axis=0)                   
ham = genericGate(bBits+t)        #make it a gate
ham.matrix = hamMat

phib = ham(phi0b)

### QFT

QFTGate = QFT(t)    ###### only phi gets qft'd
QFTGate = parallelGate([QFTGate, identity(bBits)]) 
phib = QFTGate(phib)

### ADD ANCILLA

ancilla = register(1)    
phiba = prod(phib, ancilla)

### CONTOLLED U ROTATION
# want to make a less hacky implementation here - will need to rethink controlled
#  gates in qutiepy...

toKron = [np.eye(2)] * t + [np.array([[1,0],[0,0]])] * bBits + [np.array([[0,-1],[1,0]])]  #?????

res = toKron[0]
for m in toKron[1:]:
    res = np.kron(res, m)
    
contRotGate = genericGate(t + bBits + 1)
contRotGate.matrix = res

phiba = contRotGate(phiba)

### iQFT

iQFTGate = QFTGate.H()
iQFTGate = parallelGate([iQFTGate, identity(1)])

phiba = iQFTGate(phiba)

### INVERSE HAMILTONIAN

iham = ham.H()
iham = parallelGate([iham, identity(1)])

phiba = iham(phiba)

### OBSERVE ANCILLA

print(phiba.observe(bit=t+bBits))




















