from BlackOil import GasPhase_Correlations
from BlackOil import OilPhase_Correlations

T = 122
TR = 122 + 460
P = 3626
Pb = 5000
do = 0.86
dg = 0.84

print('### Gas ###')

print('Massa específica do gás: ',GasPhase_Correlations(dg,P,TR).rho_g()) # bateu
print('Fator Z: ', GasPhase_Correlations(dg, P, TR).Z()) # bateu
print('Compressibilidade do gás: ', GasPhase_Correlations(dg, P, TR).Cg()) # bateu
print('Viscosidade do gás: ', GasPhase_Correlations(dg, P, TR).mu_lee1966()) # N bateu
print('Bg: ', GasPhase_Correlations(dg, P, TR).Bg()) # bateu

print('')
print('### Oil ###')

print('Razão de solubilidade gás-óleo: ', OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).Rs()) # bateu
print('Fator volume- formação do óleo: ', OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).Bo()) # bateu
print('Massa específica do óleo: ', OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).rho_o()) # bateu
print('Compressibilidade do óleo: ', OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).Co()) # será q bateu?
print('Viscosidade do óleo morto: ', OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).mu_od()) # bateu
print('Viscosidade do óleo saturado: ', OilPhase_Correlations(do=do, dg=dg, API=None, P=P, Pb=Pb, T=T).mu_ob()) # bateu