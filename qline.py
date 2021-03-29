
from qutiepy import *
import numpy as np

"""

Cx = b

"""

C = np.eye(2)

t = 4

b = register(t)  # 2 bits per value in b, value of B is now [8, 0]?
T = 2 ** t
amps = np.sqrt(2/T) * np.array([np.sin((np.pi*(tau+0.5)/T)) for tau in range(T)])
phi0 = register(t)
phi0.setAmps(amps)


print(phi0)

phi0b = prod(phi0, b)

t0 = 1
for tau in range(T):
    amps = np.zeros(T, dtype=np.dtype(complex))
    amps[tau] = phi0b.amps[tau]
    print(amps)
    kettau = register(t)
    kettau.setAmps(amps)
    tautau = prod(kettau, kettau)
    
