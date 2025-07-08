import numpy as np

#from CompositionalModel import Rachford_Rice
from CompositionalModel import Flash_Algorithm

# Exercício Equação de estado cúbica
T = 220 # K
P = 30*10**5 # Pa
R = 8.314
x = np.array([0.2457,0.7443])
y = np.array([0.7682,0.2318])
Pc = np.array([45.99*10**5, 73.74*10**5]) # Pa
Tc = np.array([190.56,304.12]) # K
w = np.array([0.011,0.225])
M = np.array([16.043,44.01]) # g/mol
kij = 0.0259
z = np.array([0.5,0.5])
K = np.array([3.1266,0.3073])


x, y, _, _, _, _ = Flash_Algorithm(T=T, Tc=Tc, P=P, Pc=Pc, x=x, y=y, z=z, EEC='PR', w=w, K=K, R=R, kij=kij, M=M).execute()


#Flash_Algorithm(T=T, Tc=Tc, P=P, Pc=Pc, x=x, y=y, z=z, EEC='PR', w=w, K=K, R=R, kij=kij, M=M).execute()
