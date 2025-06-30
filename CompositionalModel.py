import numpy as np

def scarvorought1966(n):
    return 0.5*(10**(2-n)); ## Critério de parada


def f_epest(v_new,v_old):
    return abs(v_old-v_new)
    #return abs((v_new - v_old)/v_new)*100 ## Erro Percentual Estimado


def newtonraphson(xo, f, df, Eppara):
    Epest = 100
    x_new = 0
    x_old = xo
    n = 0
    erro = []
    roots = []
    while Epest > Eppara:
        x_new = x_old - f(x_old)/ df(x_old)
        Epest = f_epest(x_new, x_old)
        erro.append(Epest)
        roots.append(x_new)
        x_old = x_new
        n+=1

    
    return roots[-1], erro


class CompositionalModel:

    def __init__(self, T=None, Tc=None, P=None, Pc=None, z=None, xc=None, yc=None, EEC=None, PM=None, Vc=None, Zc=None, w=None, K=None, Vchute=None, calculateonly=None):
        self.T, self.Tc, self.P, self.Pc, self.z = T, Tc, P, Pc, z
        self.eec, self.PM, self.Vc, self.Zc, self.w = EEC, PM, Vc, Zc, w
        self.K = K
        self.R = 8.314
        self.xc, self.yc = xc, yc
        if EEC == 'PR' or EEC == 'Peng-Robinson':
            self.omegaa = 0.45724
            self.omegab = 0.0778
            self.delta1 = 1 + np.sqrt(2)
            self.delta2 = 1 - np.sqrt(2)
        elif EEC == 'SRK' or EEC == 'Soave-Redlich-Kwong':
            self.omegaa = 0.42748
            self.omegab = 0.08664
            self.delta1 = 0
            self.delta2 = 0
        else:
            raise Exception('EEC deve ser SRK (Soave-Redlich-Kwong) ou PR (Peng-Robinson)!')

        if Vchute is not None:
            self.Vchute = Vchute
        else:
            self.Vchute = 0.5

        if calculateonly == 'Rachford-Rice':
            if K is None:
                self.K()
            self.calculate_LV()
            print('O V calculado é: ',self.V)
            print('O L calculado é: ',self.L)
        elif xc is not None and yc is not None:
            if len(self.xc) != len(self.yc):
                raise Exception('xc e yc devem ter o mesmo tamanho!')
            else:
                self.EEC()

        else:
            self.execute()

    def execute(self):
        self.K()
        self.calculate_LV()
        self.calculate_molar_fractions()
        self.EEC()

    def K(self):
        self.K = np.zeros(len(self.w))
        for i in range(len(self.w)):
            self.K[i] = self.Pc[i]/self.P*np.exp(5.37*(1+self.w[i])*(1-self.Tc[i]/self.T))

    def Rachford_Rice(self, V):
        soma = 0
        for i in range(len(self.K)):
            soma += self.z[i]*(self.K[i]-1)/(1-V+V*self.K[i])
        
        return soma

    def dRachford_Rice(self, V):
        soma = 0
        for i in range(len(self.K)):
            soma += self.z[i]*(self.K[i]-1)**2/(1+V*(self.K[i]-1))**2
    
        return -soma
    
    def calculate_LV(self):
        self.V, _ = newtonraphson(self.Vchute, self.Rachford_Rice, self.dRachford_Rice, 10**-3)
        self.L = 1 - self.V

    def calculate_molar_fractions(self):
        self.xc, self.yc = np.zeros(len(self.z)), np.zeros(len(self.z))
        for i in range(len(self.z)):
            self.xc[i] = self.z[i]/(1+self.V*(self.K[i]-1))
            self.yc[i] = self.K[i]*self.z[i]/(1+self.V*(self.K[i]-1))

    ### EEC ###

    def M(self,w):
        if self.eec == 'PR' or self.eec == 'Peng-Robinson':
            if w<0.49:
                return 0.37464 + 1.54226*w - 0.26992*w**2
            else:
                return 0.379642 + 1.48503*w - 0.1644236*w**2 + 0.16667*w**3
        else:
            return 0.480 + 1.574*w - 0.176*w**2

    def mixture(self):
        kij = 0 # Coeficiente de interação química
        prod_x = 1
        bgas, bliq = 0, 0
        agas, aliq = 0, 0
        a = np.zeros(len(self.xc))
        b = np.zeros(len(self.yc))

        for i in range(len(self.xc)):
            prod_x = prod_x*self.xc[i]
        Aij = np.sqrt(prod_x)*(1-kij)

        ####### QUEM É pri, Tri????
        for i in range(len(self.yc)):
            b[i] = self.omegab*self.Pr[i]/self.Tr[i]
        for i in range(len(self.xc)):
            a[i] = self.omegaa*self.Pr[i]/self.Tr[i]**2 * (1+self.M(self.w[i])*(1-np.sqrt(self.Tr[i])))**2

        for i in range(len(self.yc)):
            for j in range(len(self.xc)):
                if i == j:
                    aliq += self.xc[i]*self.xc[j]
                    agas += self.yc[i]*self.yc[j]
                else:
                    aliq += self.xc[i]*self.xc[j]*Aij
                    agas += self.yc[i]*self.yc[j]*Aij

        for i in range(len(self.yc)):
            bliq += self.xc[i] * b[i]
            bgas += self.yc[i] * b[i]

        self.B[0] = bliq * self.P / (self.R * self.T)
        self.B[1] = bgas * self.P / (self.R * self.T)
        self.A[0] = aliq * self.P / (self.R * self.T)**2
        self.A[1] = agas * self.P / (self.R * self.T)**2

    def Gibbs(self, Z):
        self.gibbs = ((Z-1)-np.log(Z-self.B)
                      -self.A/((self.delta1-self.delta2)*self.B)*np.log((Z+self.delta1*self.B)/(Z+self.delta2*self.B)))

    def Zcoef(self, A, B, delta1, delta2):
        self.C1 = 1 # Z**3
        self.C2 = (delta1 + delta2 - 1)*B-1 # Z**2
        self.C3 = ((A + delta1*delta2 * B**2) - (delta1+delta2)*B*(B+1)) # Z**1
        self.C4 = (A*B + delta1*delta2*(B**2)*(B+1))

    def Zfunction(self, Z):
        return self.C1*Z**3 + self.C2*Z**2 + self.C3*Z + self.C4

    def dZfunction(self, Z):
        return 3*self.C1*Z**2 + 2*self.C2*Z + self.C3
    
    def EEC(self):
        A, B = self.mixture()
        xo = 0
        self.Zcoef(A,B, self.delta1, self.delta2)
        Z1, _ = newtonraphson(xo, f=self.Zfunction, df=self.dZfunction, Eppara=10**-6)
        # Regra de Briot-Ruffini #
        a = self.C1
        b = a*Z1+self.C2
        c = b*Z1+self.C3
        #d = c*Z1+C3
        delta = b**2 - 4*a*c
        if delta == 0:
            Z2 = -b/(2*a)
            if self.gibbs(Z2)>self.gibbs(Z1):
                self.Z = self.Z1
            else:
                self.Z = Z2
        elif delta > 0:
            Z2 = (-b+np.sqrt(delta))/(2*a)
            Z3 = (-b-np.sqrt(delta))/(2*a)
            Zmin = min([Z1, Z2, Z3])
            Zmax = max([Z1, Z2, Z3])
            if self.gibbs(Zmin)>self.gibbs(Zmax):
                self.Z = self.Zmax
            else:
                self.Z = Zmin
        else:
            self.Z = Z1

    def calculateZ(self):
        pass

    def fugacidade(self):
        for i in range(len(self.x)):
            self.fL, self.self.fV = np.zeros(len(self.x)), np.zeros(len(self.y))
            self.fL[i] = self.phiL[i]*self.x[i]*self.P
            self.fV[i] = self.phiV[i]*self.y[i]*self.P

    def teste(self):
        soma = 0
        for i in range(len(self.y)):
            soma += (self.fL[i]/self.fV[i] - 1)**2
        if soma < 10**-8:
            for i in range (len(self.y)):
                self.K[i] = self.y[i]/self.x[i]
                print('Resultados: ', self.V, self.L, self.x, self.y, self.ZL, self.ZV, self.K)
        else:
            for i in range(len(self.y)):
                self.K[i] = self.K[i]*(self.fL[i]/self.fV[i])