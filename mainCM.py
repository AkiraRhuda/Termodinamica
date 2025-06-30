import numpy as np

from CompositionalModel import CompositionalModel

# Exercício Rachford-Rice

x = np.array([0.2457,0.7443])
y = np.array([0.7682,0.2318])
K = np.array([3.1266,0.3073])
z = np.array([0.5,0.5])

CompositionalModel(EEC='PR', z=z, K=K, calculateonly='Rachford-Rice')

# Exercício Equação de estado cúbica
T = 220
P = 30
x = np.array([0.2457,0.7443])
y = np.array([0.7682,0.2318])
Pc = np.array([3.1266,0.3073])
Tc = np.array([190.56,304.12])
w = np.array([0.011,0.225])
Tc = np.array([16.043,44.01])

CompositionalModel(T=T, Tc=Tc, P=P, Pc=Pc, xc=x, yc=y, EEC='PR', w=w)