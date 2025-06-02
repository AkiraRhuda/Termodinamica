from BlackOilModel import GasPhase
from BlackOilModel import OilPhase
from barril.units import Scalar
from barril.units import ChangeScalars


# Dados iniciais
densidaderelativagas = 0.84
densidaderelativaoleo = 0.86
tempseparador = 80 # °F
Pressaoseparador = 100 # psia
Pressaodebolha = 5000 # psia
pressaolinha = 3626 # psia
temperaturalinhaR = 122 + 460 # R
temperaturalinha = 122 # F
#temperaturalinha = Scalar(122, "F") # °F
#temperaturalinhaconvertida = ChangeScalars(temperaturalinha, "R")

GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=densidaderelativagas, P=pressaolinha, T=temperaturalinhaR).output()
OilPhase("Standing", "Standing", "Beggs",densidaderelativaoleo, densidaderelativagas, pressaolinha, Pressaodebolha, temperaturalinha).output()

