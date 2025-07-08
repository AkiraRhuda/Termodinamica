import numpy as np

#from CompositionalModel import Rachford_Rice
from CompositionalModel import Flash_Algorithm
from CompositionalModel import PlotIsotermico
from CompositionalModel import PlotIsobarico
# Exercício Equação de estado cúbica
dict = ['C1', 'CO2']
T = 220 # K
#P = np.array([5*10**5,15*10**5,30*10**5, 60*10**5, 100*10**5, 150*10**5, 200*10**5, 300*10**5]) # Pa
#P = np.array([15*10**5,20*10**5, 25*10**5, 30*10**5, 35*10**5, 40*10**5, 45*10**5, 50*10**5, 55*10**5, 60*10**5]) # Pa
P = np.array([15*10**5,20*10**5, 25*10**5, 30*10**5, 35*10**5, 40*10**5, 45*10**5, 50*10**5]) # Pa
#P = np.array([30*10**5])
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

X = []
Y = []

for i in range(len(P)):
    xi, yi, l, o = Flash_Algorithm(T=T, Tc=Tc, P=P[i], Pc=Pc, x=x, y=y, z=z, EEC='PR', w=w, K=K, R=R, kij=kij, M=M).execute()
    X.append(xi)
    Y.append(yi)


PlotIsotermico(X, Y, P, T, dict)
P = 30*10**5
T = np.array([200, 205, 210, 215, 220, 225, 230, 235]) # K
X = []
Y = []
for i in range(len(T)):
    xi, yi, l, o = Flash_Algorithm(T=T[i], Tc=Tc, P=P, Pc=Pc, x=x, y=y, z=z, EEC='PR', w=w, K=K, R=R, kij=kij, M=M).execute()
    X.append(xi)
    Y.append(yi)

PlotIsobarico(X, Y, P, T, dict)
print(X)