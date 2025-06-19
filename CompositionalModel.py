import numpy as np

def scarvorought1966(n):
    return 0.5*(10**(2-n)); ## Critério de parada


def f_epest(v_new,v_old):
    return abs((v_new - v_old)/v_new)*100; ## Erro Percentual Estimado


def newtonraphson(xo, f, df):
    Epest = 100
    x_new = 0
    x_old = xo
    Eppara = scarvorought1966(6)
    n = 0
    erro = []
    roots = []
    while Epest > Eppara and n <= 100:
        x_new = x_old - f(x_old)/ df(x_old)
        Epest = f_epest(x_new, x_old)
        erro.append(Epest)
        roots.append(x_new)
        x_old = x_new
        n+=1

    
    return roots, erro


class CompositionalModel:

    def __init__(self, T, Tc, P, Pc, z, EEC, PM, Vc, Zc, w):
        self.T, self.Tc, self.P, self.Pc, self.z = T, Tc, P, Pc, z
        self.EEC, self.PM, self.Vc, self.Zc, self.w = EEC, PM, Vc, Zc, w

    def K(self):
        self.K = self.Pc/self.P*np.exp(5.37*(1+self.w)*(1-self.Tc/self.T))

    def Rachford_Rice(self):
        soma = 0
        for i in range(len(Tc)):
            soma += z[i]*(K[i]-1)/(1-V+V*K[i])
        
        return soma

    def dRachford_Rice(self):
        soma = 0
        for i in range(nc):
            soma += z[i]*(K[i]-1)**2/(1+V*(K[i]-1))**2
    
        return -soma
    
    def LV(self):
        self.V, _ = newtonraphson(self.K, self.Rachford_Rice(), self.dRachford_Rice())
        self.L = 1 - self.V

    def xiyi(self):
        self.x = self.z/(1+self.V*(self.K-1))
        self.y = self.K*self.z/(1+self.V*(self.K-1))

    def coefficients(self):
        omegaa = 0.45724 
        omegab = 0.0778
        
    
        b = omegab
        #Bliq += self.B[i]*self.y[i] for i in range(len(self.y))
    
    def Gibbs(self, Z):
        return (Z-1)*np.log(Z-B) - A/((delta1-delta2)*B)*np.log((Z+delta1*B)/(Z+delta2*B))

    @staticmethod
    def Zcoef(A, B):
        delta1 = 1+ np.sqrt(2)
        delta2 = 1 - np.sqrt(2)
        C1= 1 # Z**3
        C2 = (B-1)*(delta1 + delta2 - 1) # Z**2
        C3 = ((A + delta1*delta2 * B**2) - (delta1+delta2)*B*(B+1))
        C4 = (A*B + delta1*delta2*(B**2)*(B+1))
        return C1, C2, C3, C4
    
    @staticmethod
    def Zfunc(Z, C1, C2, C3, C4):
        return C1*Z**3 + C2*Z**2 + C3*Z + C4
    
    @staticmethod
    def dZfunc(Z,C1, C2, C3):
        return 3*C1*Z**2 + 2*C2*Z + C3
    
    def Z(Z,A,B):
        xo = [0, 0, 0]
        C1, C2, C3, C4 = Zcoef(A,B)
        Z, _ = newtonraphson(xo, f=Zfunc(Z, C1, C2, C3, C4), df=dZfunc(Z, C1, C2, C3, C4))

        """if isinstance(Z[i], complex) for i in Z:
            
        else:
            zmin = min(Z)
            zmax = max(Z)
            if self.Gibbs(zmin)<self.Gibbs(zmax):
                Zfase = zmin
            else:
                Zfase = zmax
        return Zfase"""

    def fugacidade(self):
        self.fL = self.phiL*self.x*self.P
        self.fV = self.phiV*self.y*self.P

    def teste(self):
        for i in range(len(self.y)):
            soma = (self.fL[i]/self.fV[i] - 1)**2
        if soma < 10**-8:
            K = self.y/self.x
        print('Resultados: ', self.V, self.L, self.x, self.y, self.ZL, self.ZV, self.K)
        else:
            # REPETICAO
