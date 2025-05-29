import correlations
from barril.units import Scalar
from barril.units import ChangeScalars

densidaderelativagas = 0.84
densidaderelativaoleo = 0.86
tempseparador = 80 # °F
Pressaoseparador = 100 # psia
Pressaodebolha = 5000 # psia

pressaolinha = 3626 # psia
temperaturalinha = 122
temperaturalinha = 122 + 460 # R
#temperaturalinha = Scalar(122, "F") # °F
#temperaturalinhaconvertida = ChangeScalars(temperaturalinha, "R")

Z, Cg, μ = correlations.FaseGas(correlation="Papay", dg=densidaderelativagas, P=pressaolinha, T=temperaturalinha).output()

print('Fator Z: ',Z)
print('Compressibilidade do gás: ',Cg)
print('Viscosidade do gás: ', μ)
