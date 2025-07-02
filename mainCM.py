import numpy as np

#from CompositionalModel import Rachford_Rice
from CompositionalModel import Flash_Algorithm

# Exercício Rachford-Rice
"""
x = np.array([0.2457,0.7443])
y = np.array([0.7682,0.2318])
K = np.array([3.1266,0.3073])
z = np.array([0.5,0.5])

V, L = Flash_Algorithm(z=z, K=K, EEC='PR').calculate_LV()


print('O V calculado é: ',V)
print('O L calculado é: ',L)
"""
# Exercício Equação de estado cúbica
T = 220
P = 30
x = np.array([0.2457,0.7443])
y = np.array([0.7682,0.2318])
Pc = np.array([3.1266,0.3073])
Tc = np.array([190.56,304.12])
w = np.array([0.011,0.225])
M = np.array([16.043,44.01])
kij = np.array([0.0259])
z = np.array([0.5,0.5])
K = np.array([3.1266,0.3073])

Flash_Algorithm(T=T, Tc=Tc, P=P, Pc=Pc, x=x, y=y, z=z, EEC='PR', w=w, K=K, kij=kij, M=M).execute()
