from BlackOil import GasPhase_Correlations
from BlackOil import OilPhase_Correlations

# Dados do problema
T = 122
TR = 122 + 460
P = 3626
Pb = 5000
do = 0.86
dg = 0.84

print('### Gas ###')
GasPhase_Correlations(dg,P,TR).output()
print('')

print('### Oil ###')
OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).output()