import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class FaseGas:
    """
    Calcula o fator de compressibilidade térmico

    Parametros
    ----------
    correlation : str
        Informa a correlação a ser usada
    P : float
        Pressao absoluta
    Ppc : float
        Pressao pseudocritica
    Ppr : float
        Pressao pseudoreduzida
    T : float
        Temperatura
    Tpc : float
        Temperatura pseudocritica
    Tpr : float
        Temperatura pseudoreduzida
    """
    def __init__(self, correlation, dg=None,P=None, Ppc=None, Ppr=None, T=None, Tpc=None, Tpr=None):
        self.dg = dg
        self.P, self.Ppc, self.Ppr, self.T, self.Tpc, self.Tpr = P, Ppc, Ppr, T, Tpc, Tpr
        self.correlation = correlation
        self.Cpr, self.Cg, self.Z = None, None, None
        self.initialproperties()
        self.functionselector()
        self.compressisoterm()
        self.dempsey()


    def initialproperties(self):
        if self.dg is not None and self.Ppc is None and self.Tpc is None:
            if self.dg < 0.75:
                self.Ppc = 677 + 15*self.dg - 37.5*self.dg**2
                self.Tpc = 168 + 325*self.dg - 12.5*self.dg**2
            else:
                self.Ppc = 706 - 51.7*self.dg - 11.1*self.dg**2
                self.Tpc = 187 + 330*self.dg - 71.5*self.dg**2
        if self.Ppr is None and self.Tpr is None:
            self.Ppr = self.P/self.Ppc
            self.Tpr = self.T/self.Tpc
        if self.dg is None:
            pass
    def functionselector(self):
        if self.correlation == "Papay":
            self.papay()

    def papay(self, Ppr=None):
        if Ppr is None:
            self.Z = 1 - (3.53*self.Ppr)/(10**(0.9813*self.Tpr)) + (0.274*self.Ppr**2)/(10**(0.8157*self.Tpr))
        else:
            return 1 - (3.53*Ppr)/(10**(0.9813*self.Tpr)) + (0.274*Ppr**2)/(10**(0.8157*self.Tpr))

    def dempsey(self):
        self.μΔgN2, self.μΔgCO2, self.μΔgH2S = 0, 0, 0
        μgncorrigida = 1.709 * 10**(-5) - 2.062*10**(-6)*self.dg*self.T + 8.188*10**(-3)-6.15*10**(-3)*np.log10(self.dg)
        self.μg = μgncorrigida + self.μΔgN2 + self.μΔgCO2 + self.μΔgH2S
        a0, a1, a2, a3, a4, a5 = -2.4621, 2.9705, -2.8626e-1, 8.0542e-3, 2.8086, -3.4980
        a6, a7, a8, a9, a10, a11 = 3.6037e-1, -1.0443e-2, -7.9339e-1, 1.3964, 1.4914e-1, 4.4102e-3
        a12, a13, a14, a15 = 8.3939e-2, -1.8641e-1, 2.0336e-2, -6.0958e-4

        lnexpression = a0 + a1*self.Ppr + a2*self.Ppr**2 + a3*self.Ppr**3 + self.Tpr*(a4 + a5*self.Ppr + a6*self.Ppr**2
                    + a7*self.Ppr**3) + self.Tpr**2*(a8 + a9*self.Ppr + a10*self.Ppr**2 +a11*self.Ppr**3) \
                    + self.Tpr**3*(a12 + a13*self.Ppr + a14*self.Ppr**2 + a15*self.Ppr**3)
        self.μ = np.exp(lnexpression)*self.μg*1/self.Tpr

    def compressisoterm(self):
        Ppr = self.Ppr
        dx = 10**-10
        self.Cpr = 1/self.Ppr - 1/self.Z * (self.papay(Ppr+dx)-self.papay(Ppr))/dx
        self.Cg = self.Cpr/self.Ppc

    def output(self):
        return self.Z, self.Cg, self.μ
class FaseOleo:
    def __init__(self):
        pass