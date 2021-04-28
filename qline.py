
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


def main(debug=False):
    A = np.array([[1,0],
                  [0,1]])
    
    k = np.linalg.cond(A)
    print("k = ", k)
    
    bBits = int(np.log2(A.shape[0]))
    bAmps = [1, -1]
    
    b = register(bBits)
    b.setAmps(bAmps)
    t = 8  # bits in phi
    T = 2 ** t    # states in phi
    amps = np.sqrt(2/T) * np.array([np.sin((np.pi*(tau+0.5)/T)) for tau in range(T)])
    phi0 = register(t)
    phi0.setAmps(amps)
    phi0b = prod(phi0, b)
    
    t0 = 1
    
    ###     HAMILTONIAN SIMULATION
    
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
    
    if debug:
        print("Hamiltonian:\n", phib)
        input()
    
    ### QFT
    
    QFTGate = QFT(t)    ###### only phi gets qft'd
    QFTGate = parallelGate([QFTGate, identity(bBits)]) 
    phib = QFTGate(phib)
    
    if debug:
        print("QFT:\n", phib)
        input()
    
    ### ADD ANCILLA
    
    ancilla = register(1)    
    phiba = prod(phib, ancilla)
    
    if debug:
        print("Add Ancilla:\n", phiba)
        input()
    
    ### CONTOLLED U ROTATION
    
    # if phi = x, apply Rx(f(x)) on ancilla
    
    gatesInSeries = []
    for x in range(1, T):
        xBin = f'{x:0{t}b}'
        preGates = []
        for bit in xBin:
            if bit == "0":
                preGates.append(pauliX(1))
            else:
                preGates.append(identity(1))
        
        preGates.append(identity(bBits + 1))
        preGate = parallelGate(preGates)
        
        theta = np.arccos(0.0035/((x/T))) ## I think...?
        controlledGate = genericGate(1)
        s = np.sin(theta/2)
        c = np.cos(theta/2)
        controlledGate.matrix = np.array([[c,-s*1j],[-s*1j,c]])
        offsets = list(range(-bBits-1, -(bBits+t)-1, -1))
        controlledGate = controlledGate.addControlBits(offsets)
        
        postGate = preGate.H()
        
        gatesInSeries.append(postGate(controlledGate(preGate)))
    
    contRot = gatesInSeries[0]
    for gate in gatesInSeries[1:]:
        contRot = gate(contRot)
    
    phiba = contRot(phiba)
    
    if debug:
        print("Controlled Rotation:\n", phiba)
        input()
    
    ### iQFT
    
    iQFTGate = QFTGate.H()
    iQFTGate = parallelGate([iQFTGate, identity(1)])
    
    phiba = iQFTGate(phiba)
    
    if debug:
        print("iQFT:\n", phiba)
        input()
    
    ### INVERSE HAMILTONIAN
    
    iham = ham.H()
    iham = parallelGate([iham, identity(1)])
    
    phiba = iham(phiba)
    
    if debug:
        print("Inv. Hamiltonian:\n", phiba)
        input()
    
    ### OBSERVE ANCILLA
    
    ancilla = phiba.observe(bit=phiba.NBits-1)
    print(ancilla)
    
    if ancilla:
        print(f'{0:0{phiba.NBits-1}0b}'.format(phiba.observe()))
        return True
    
    else:
        return False
        
if __name__ == "__main__":
    while not main(False):
        pass



















