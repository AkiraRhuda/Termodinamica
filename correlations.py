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
        self.prints()


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
            self.papay(self)

    def papay(self, Ppr=None):
        if Ppr is None:
            self.Z = 1 - (3.53*self.Ppr)/(10**(0.9813*self.Tpr)) + (0.274*self.Ppr**2)/(10**(0.8157*self.Tpr))
        else:
            return 1 - (3.53*Ppr)/(10**(0.9813*self.Tpr)) + (0.274*Ppr**2)/(10**(0.8157*self.Tpr))

    def dempsey(self):
        #μgncorrigida = 1.709 * 10 **(-5) - 2.062*10**(-6)*self.dg*T + 8.188*10**(-3)-6.15*10**(-3)*np.log10(self.dg)
        #self.μg = μgncorrigida + μ
        pass

    def compressisoterm(self):
        Ppr = self.Tpr
        dx = 0.0001
        self.Cpr = 1/self.Ppr - 1/self.Z * (self.papay(Ppr+dx)-self.papay(Ppr))/dx
        self.Cg = self.Cpr/self.Ppc

    def prints(self):
        print(self.Z)
        print(self.Cg)
