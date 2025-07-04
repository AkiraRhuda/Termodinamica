import numpy as np
import matplotlib.pyplot as plt

def scarvorought1966(n):
    return 0.5*(10**(2-n)); ## Critério de parada

def f_epest(v_new,v_old):
    return abs((v_new - v_old)/v_new)*100 ## Erro Percentual Estimado


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


class Rachford_Rice:
    def __init__(self, z, K, Vchute=None):
        self.K = K
        self.z = z
        if Vchute is not None:
            self.Vchute = Vchute
        else:
            self.Vchute = 0.5

    def rachford_rice(self, V):
        soma = 0
        for i in range(len(self.K)):
            soma += self.z[i] * (self.K[i] - 1) / (1 - V + V * self.K[i])

        return soma

    def drachford_rice(self, V):
        soma = 0
        for i in range(len(self.K)):
            soma += self.z[i] * (self.K[i] - 1) ** 2 / (1 + V * (self.K[i] - 1)) ** 2

        return -soma

    def calculate_LV(self):
        V, _ = newtonraphson(self.Vchute, self.rachford_rice, self.drachford_rice, 10 ** -3)
        L = 1 - V
        return V, L


class Flash_Algorithm:

    def __init__(self, T=None, Tc=None, P=None, Pc=None, z=None, x=None, y=None, EEC=None, kij=None, PM=None, Vc=None, Zc=None, w=None, M=None, K=None, Vchute=None, a=None, b=None, param=None):
        self.T, self.Tc, self.P, self.Pc, self.z = T, Tc, P, Pc, z
        self.eec, self.PM, self.Vc, self.Zc, self.w = EEC, PM, Vc, Zc, w
        self.M, self.kij, self.K, self.Vchute = M, kij, K, Vchute
        self.R = 8.314
        self.x, self.y = x, y
        self.a, self.b = a, b
        self.it = 0
        if kij is not None:
            self.kij = kij
        else:
            self.kij = 0

        if Vchute is None:
            self.Vchute = 0.5

        if self.Tc is not None and self.Pc is not None:
            if len(self.Pc) == len(self.Tc):
                self.Pr = np.zeros(len(self.Pc))
                self.Tr = np.zeros(len(self.Tc))
                for i in range(len(self.Pc)):
                    self.Pr[i] = self.P/self.Pc[i]
                    self.Tr[i] = self.T/self.Tc[i]

            else:
                raise Exception('Faltam dados de temperatura ou pressão!')

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

        if param == 'mixture':
            self.mixture()

        if K is None:
            self.wilson_correlation()

        if x is not None and y is not None:
            if len(self.x) != len(self.y):
                raise Exception('x e y devem ter o mesmo tamanho!')


    def execute(self):
        self.calculate_LV()
        if self.it > 0:
            self.calculate_molar_fractions()
        self.EEC()
        self.fugacidade()
        self.verification()

    def wilson_correlation(self):
        self.K = np.zeros(len(self.w))
        for i in range(len(self.w)):
            self.K[i] = self.Pc[i]/self.P*np.exp(5.37*(1+self.w[i])*(1-self.Tc[i]/self.T))

    def Rachford_Rice(self, V):
        soma = 0
        for i in range(len(self.z)):
            soma += self.z[i]*(self.K[i]-1)/(1-V+V*self.K[i])
        
        return soma

    def dRachford_Rice(self, V):
        soma = 0
        for i in range(len(self.z)):
            soma += self.z[i]*(self.K[i]-1)**2/(1+V*(self.K[i]-1))**2
    
        return -soma
    
    def calculate_LV(self):

        if self.it == 0:
            self.V, _ = newtonraphson(self.Vchute, self.Rachford_Rice, self.dRachford_Rice, 10**-6)
        else:
            self.V, _ = newtonraphson(self.V, self.Rachford_Rice, self.dRachford_Rice, 10 ** -6)
        self.L = 1 - self.V
        return self.V, self.L

    def calculate_molar_fractions(self):
        self.x, self.y = np.zeros(len(self.z)), np.zeros(len(self.z))
        for i in range(len(self.z)):
            self.x[i] = self.z[i]/(1+self.V*(self.K[i]-1))
            self.y[i] = self.K[i]*self.z[i]/(1+self.V*(self.K[i]-1))

    ### EEC ###

    def m(self,w):
        if self.eec == 'PR' or self.eec == 'Peng-Robinson':
            if w<0.49:
                return 0.37464 + 1.54226*w - 0.26992*w**2
            else:
                return 0.379642 + 1.48503*w - 0.1644236*w**2 + 0.16667*w**3
        else:
            return 0.480 + 1.574*w - 0.176*w**2


    def calculate_a_and_b(self):
        a = np.zeros(len(self.x))
        b = np.zeros(len(self.y))
        for i in range(len(self.y)):
            b[i] = self.omegab*self.Pr[i]/self.Tr[i]
        for i in range(len(self.x)):
            a[i] = self.omegaa*self.Pr[i]/self.Tr[i]**2 * (1+self.m(self.w[i])*(1-np.sqrt(self.Tr[i])))**2
        return a, b

    def mixture(self):
        prod_a = 1
        bgas, bliq = 0, 0
        agas, aliq = 0, 0
        self.B = np.zeros(len(self.x))
        self.A = np.zeros(len(self.x))
        self.Psi = np.zeros(len(2)) # número de fases
        if self.a is None and self.b is None:
            a, b = self.calculate_a_and_b()
        else:
            a, b = self.a, self.b
        for i in range(len(self.x)):
            prod_a = prod_a*a[i]
        Aij = np.sqrt(prod_a)*(1-self.kij)

        for i in range(len(self.y)):
            for j in range(len(self.x)):
                if i == j:
                    aliq += self.x[i]*self.x[j]*a[i]
                    agas += self.y[i]*self.y[j]*a[i]
                else:
                    aliq += self.x[i]*self.x[j]*Aij
                    agas += self.y[i]*self.y[j]*Aij

        for i in range(len(self.y)):
            bliq += self.x[i] * b[i]
            bgas += self.y[i] * b[i]
        self.B[0], self.B[1] = bliq, bgas
        self.A[0], self.A[1] = aliq, agas

        for i in range(len(2)):
            for j in range(len(self.x)):
                if i == j:
                    self.Psi[i] += Aij*self.y[i]
                else:
                    self.Psi[i] += Aij*self.y[i]




    def gibbs_energy(self, Z, i):
        self.gibbs = ((Z-1)-np.log(Z-self.B[i])
                      -self.A[i]/((self.delta1-self.delta2)*self.B[i])*np.log((Z+self.delta1*self.B[i])/(Z+self.delta2*self.B[i])))
        return self.gibbs

    def Zcoef(self, A, B, delta1, delta2):
        self.C1 = 1 # Z**3
        self.C2 = (delta1 + delta2 - 1)*B-1 # Z**2
        self.C3 = ((A + delta1*delta2 * B**2) - (delta1+delta2)*B*(B+1)) # Z**1
        self.C4 = (A*B + delta1*delta2*(B**2)*(B+1))

    def Zfunction(self, Z):
        return self.C1*Z**3 + self.C2*Z**2 + self.C3*Z - self.C4

    def dZfunction(self, Z):
        return 3*self.C1*Z**2 + 2*self.C2*Z + self.C3
    
    def EEC(self):
        self.mixture()
        Z1 = np.zeros(len(self.x))
        Z2 = np.zeros(len(self.x))
        Z3 = np.zeros(len(self.x))
        self.Z = np.zeros(len(self.x))
        for i in range(len(self.x)):
            xo = 1
            self.Zcoef(self.A[i],self.B[i], self.delta1, self.delta2)
            Z1[i], _ = newtonraphson(xo, f=self.Zfunction, df=self.dZfunction, Eppara=10**-12)
            # Regra de Briot-Ruffini #
            a = self.C1
            b = a*Z1[i]+self.C2
            c = b*Z1[i]+self.C3
            #d = c*Z1+C3
            delta = b**2 - 4*a*c
            if delta == 0:
                Z2[i] = -b/(2*a)
                if Z1[i] < 0:
                    self.Z[i] = Z2[i]
                elif Z2[i] < 0:
                    self.Z[i] = Z1[i]
                else:
                    if self.gibbs_energy(Z2[i], i)>self.gibbs_energy(Z1[i], i):
                        self.Z[i] = Z1[i]
                    else:
                        self.Z[i] = Z2[i]
            elif delta > 0:
                Z2[i] = (-b+np.sqrt(delta))/(2*a)
                Z3[i] = (-b-np.sqrt(delta))/(2*a)
                Zmin = min([Z1[i], Z2[i], Z3[i]])
                Zmax = max([Z1[i], Z2[i], Z3[i]])
                roots_energy = np.array([self.gibbs_energy(Zmin, i), self.gibbs_energy(Zmax, i)])
                if roots_energy[0]>roots_energy[1]:
                    self.Z[i] = Zmax
                else:
                    self.Z[i] = Zmin
            else:
                self.Z[i] = Z1[i]

    def calculate_phi(self):


    def fugacidade(self):
        # Phi-phi
        self.fL, self.fV = np.zeros(len(self.x)), np.zeros(len(self.y))
        for i in range(len(self.x)):
            self.phiL, self.phiV =
            self.fL[i] = self.phiL*self.x[i]*self.P
            self.fV[i] = self.phiV*self.y[i]*self.P

    def verification(self):
        soma = 0
        for i in range(len(self.y)):
            soma += (self.fL[i]/self.fV[i] - 1)**2
        if soma < 10**-8 or self.it<10:
            for i in range (len(self.y)):
                self.K[i] = self.y[i]/self.x[i]
                print('Z', self.Z)
                #print('Resultados: ', self.V, self.L, self.Z,self.x, self.y, self.K)
        else:
            for i in range(len(self.y)):
                self.K[i] = self.K[i]*(self.fL[i]/self.fV[i])
                self.wilson_correlation()
                self.it += 1
                self.execute()