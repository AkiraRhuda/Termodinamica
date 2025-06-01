from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
from barril.units import Scalar
from barril.units import ChangeScalars

densidaderelativagas = 0.84
densidaderelativaoleo = 0.86
tempseparador = 80 # °F
Pressaoseparador = 100 # psia
Pressaodebolha = 5000 # psia

pressaolinha = 3626 # psia
temperaturalinha = 122 + 460 # R
#temperaturalinha = Scalar(122, "F") # °F
#temperaturalinhaconvertida = ChangeScalars(temperaturalinha, "R")

ρg, Cg, μ = GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=densidaderelativagas, P=pressaolinha, T=temperaturalinha).output()

print('Massa específica do gás: ', ρg)
print('Compressibilidade do gás: ',Cg) # bateu
print('Viscosidade do gás: ', μ) # bateu

temperaturalinha = 122 # F

ρo, Co, μod, μob = OilPhase("Standing", "Standing", "Beggs",densidaderelativaoleo, densidaderelativagas, pressaolinha, Pressaodebolha, temperaturalinha).output()

print('Massa específica do óleo: ', ρo) # bateu
print('Compressibilidade do óleo: ', Co) # será q bateu?
print('Viscosidade do óleo morto: ', μod) # bateu
print('Viscosidade do óleo saturado: ', μob) # bateu