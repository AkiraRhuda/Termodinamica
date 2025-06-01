import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class GasPhase_Correlations:

    def __init__(self, dg, P, T):
        
        if dg == None:
            self.dg = self.rho_g()/self.rho_ar
        else:
            self.dg = dg

        if self.dg < 0.75:
            self.Ppc = 677 + 15.0*self.dg - 37.5*(self.dg**2)
            self.Tpc = 168 + 325*self.dg - 12.5*(self.dg**2)
        else:
            self.Ppc = 706 - 51.7*self.dg - 11.1*(self.dg**2)
            self.Tpc = 187 + 330*self.dg - 71.5*(self.dg**2)
        self.standardconditions()
        self.T = T
        self.P = P
        self.Ppr = P/self.Ppc
        self.Tpr = T/self.Tpc
        self.Mg = self.dg * 28.96 #M ar

    def standardconditions(self):
        self.Psc = 14.7 #psia
        self.Tsc = 60 #F
    
    def rho_g(self):
        R = 10.73 # psia·ft³/ (lb·mol·°R)
        return (self.P * self.Mg)/(self.Z() * R * self.T)
    
    def Bg(self):
        return self.Psc/(self.Tsc) * self.Z() * (self.T-460)/self.P

   ## Papay (1985)

    def Z(self, Ppr=None, Tpr=None, dg=None):
        if Ppr is not None:
            self.Ppr = Ppr
        return 1 - (3.53 * self.Ppr) / (10 ** (0.9813 * self.Tpr)) + (0.274 * self.Ppr ** 2) / (10 ** (0.8157 * self.Tpr))

    def dZdP(self, Tpr=None, Ppr=None, h=0.5*10**(-12)):

        frwrd = self.Z(self.Ppr + h, self.Tpr)
        bckwrd = self.Z(self.Ppr - h, self.Tpr)
        return (frwrd - bckwrd) / (h)


    def mu_lee1966(self):
        
            x_v = 3.448 + (986.4 / self.T) + (0.01009 * self.Mg)
            y_v = 2.4 - (0.2 * x_v)
            k_v = ((9.379 + 0.0160 * self.Mg) * (self.T ** 1.5)) / (209.2 + (19.26 * self.Mg) + self.T)

            exp = x_v * ((self.rho_g() / 62.4) ** y_v)

            return 1e-4 * k_v * np.exp(exp)
## Dempsey (1965)

    def mu_dempsey(self, dg=None, Tpr=None, Ppr=None):

        mu_g = (1.709e-5 - 2.062e-6 * self.dg) * (self.T-460) + 8.188e-3 - 6.15e-3 * np.log10(self.dg)
        print('mu_g: ', mu_g)
        a0 = -2.4621
        a1 = 2.9705
        a2 = -2.8626e-1
        a3 = 8.0542e-3
        a4 = 2.8086
        a5 = -3.4980
        a6 = 3.6037e-1
        a7 = -1.0443e-2
        a8 = -7.9339e-1
        a9 = 1.3964
        a10 = 1.4914e-1
        a11 = 4.4102e-3
        a12 = 8.3939e-2
        a13 = -1.8641e-1
        a14 = 2.0336e-2
        a15 = -6.0958e-4

        poly = (
                a0 +
                (a1 * self.Ppr) +
                (a2 * (self.Ppr ** 2)) +
                (a3 * (self.Ppr ** 3)) +
                (self.Tpr * (a4 + (a5 * self.Ppr) + (a6 * (self.Ppr ** 2)) + (a7 * (self.Ppr ** 3)))) +
                (self.Tpr ** 2) * (a8 + (a9 * self.Ppr )+ (a10 * (self.Ppr ** 2)) + (a11 * (self.Ppr ** 3))) +
                (self.Tpr ** 3) * (a12 + (a13 * self.Ppr) + (a14 * (self.Ppr ** 2)) + (a15 * (self.Ppr ** 3)))
                 )


        mu = mu_g * np.exp(poly) / self.Tpr
        return mu

    def rho(self, P, Mg, Z, T):
        R = 0
        self.ρ = P * self.Mg / Z * R * T

    def Cpr(self, Ppr=None, dZdP=None):
        return 1/self.Ppr - (1/self.Z() * self.dZdP())
    def Cg(self):
        return self.Cpr()/self.Ppc

class OilPhase_Correlations:

    def __init__(self, do, dg, API, P, Pb, T):

        self.Pb = Pb
        self.P = P
        self.T = T
        self.dg = dg

        if do == None:
            self.API = API
            self.do = 141.5/(API + 131.5)
        else:
            self.do = do
            self.API = (141.5/do) - 131.5

        if Pb == None:
            self.Pb = self.bubblepressure()
        else:
            self.Rs()
    
    ## Standing
    def bubblepressure(self):

        a = 0.00091 * self.T - 0.0125 * self.API

        return 18.2 * ((self.Rs())/self.dg) * (10**a) - 1.4


    def Rs(self):

        if self.P > self.Pb:
            a = 0.0125*self.API - 0.00091*self.T
            return self.dg * (((self.Pb/18.2) + 1.4) * 10**a)**(1/0.83)
        else:

            a = 0.0125*self.API - 0.00091*self.T
            return self.dg * (((self.P/18.2) + 1.4) * 10**a)**(1/0.83)

    def dRsdP(self):
        a = 0.0125*self.API - 0.00091*self.T
        return self.dg * (((1/18.2) + 1.4) * 10**a)**(1/0.83)
    
    def dBodP(self):
        return 0.00012*(self.dRsdP() * (1.2)*(self.dg/self.do)**(0.5) + 1.25*self.T)**(0.2)
    
    def Co(self):
        if self.P >= self.Pb:
            a = self.rho_ob() + 0.004347 * (self.P - self.Pb) - 79.1
            b = 0.0007141 * (self.P - self.Pb) - 12.938
            return 1e-6 * np.exp(a/b)
        else:
            return -1/self.Bo() * self.dBodP() + GasPhase_Correlations(self.dg, self.P, self.T+460).Bg()/self.Bo() * self.dRsdP()
        
    def Bo(self):
        if self.P > self.Pb:
            return self.Bob * np.exp(-self.Co() * (self.P - self.Pb))

        else:
            return 0.9759 + 0.00012*(self.Rs() * (self.dg/self.do)**(0.5) + 1.25*self.T)**(1.2)
            
    # Fim das correlações Standing
    def rho_o(self):
        if self.P > self.Pb:
            return self.rho_ob() * np.exp(self.Co()*(self.P-self.Pb))

        else:
            return (62.4 * self.do + 0.0136 * self.Rs() * self.dg)/self.Bo()
        
    def rho_ob(self):
        return Exception('Não implementado')
        
    # Correlação de Beggs e Robinson
    def mu_od(self):
        A = 10 ** (3.0324 - 0.02023 * self.API)
        return (10 ** (A * self.T ** (-1.163))) - 1
    
    def mu_ob(self):
        a = 10.715 * (self.Rs() + 100)**(-0.515)
        b = 5.44 * (self.Rs() + 150)**(-0.338)
        return a * self.mu_od()**b