import math

# Dados fornecidos
d_g = 0.84  # densidade do gás
T = 122  # temperatura em °F
P = 3626  # pressão em psia

# 1. Calcular propriedades pseudocríticas (Standing, 1977)
P_pc = 677 + 15 * d_g - 37.5 * d_g**2  # psia
T_pc = 168 + 325 * d_g - 12.5 * d_g**2  # °R (Rankine)

# 2. Converter temperatura de °F para °R e calcular P_pr e T_pr
T_rankine = T + 459.67
P_pr = P / P_pc
T_pr = T_rankine / T_pc

# 3. Coeficientes para correção de viscosidade (Lee-Gonzalez-Eakin)
a = [2.46211820, 2.97054714, -0.28626405, 0.00805420,
     2.80860949, -3.49803305, 0.36037302, -0.01044324,
     -0.79338568, 1.39643306, -0.14914493, 0.00441016,
     0.08393872, -0.18640885, 0.02033679, -0.00060958]

# 4. Calcular viscosidade do gás a pressão atmosférica (µg*)
def calcular_mu_g_estrela(d_g, T):
    termo1 = (1.709e-5 - 2.062e-6 * d_g) * T
    termo2 = 8.188e-3 - 6.15e-3 * math.log10(d_g)
    return termo1 + termo2

# 5. Calcular viscosidade corrigida para alta pressão (µg)
def calcular_mu_g(mu_g_estrela, P_pr, T_pr, a):
    X = (a[0] + a[1]*P_pr + a[2]*P_pr**2 + a[3]*P_pr**3 +
         T_pr * (a[4] + a[5]*P_pr + a[6]*P_pr**2 + a[7]*P_pr**3) +
         T_pr**2 * (a[8] + a[9]*P_pr + a[10]*P_pr**2 + a[11]*P_pr**3) +
         T_pr**3 * (a[12] + a[13]*P_pr + a[14]*P_pr**2 + a[15]*P_pr**3))
    return (math.exp(X) * mu_g_estrela) / T_pr

# 6. Executar cálculos
mu_g_estrela = calcular_mu_g_estrela(d_g, T)
mu_g = calcular_mu_g(mu_g_estrela, P_pr, T_pr, a)

# 7. Exibir resultados
print("=== Resultados ===")
print(f"Pressão pseudocrítica (P_pc): {P_pc:.2f} psia")
print(f"Temperatura pseudocrítica (T_pc): {T_pc:.2f} °R")
print(f"Pressão pseudoreduzida (P_pr): {P_pr:.4f}")
print(f"Temperatura pseudoreduzida (T_pr): {T_pr:.4f}")
print(f"Viscosidade do gás a pressão atmosférica (µg*): {mu_g_estrela:.6f} cP")
print(f"Viscosidade do gás corrigida (µg): {mu_g:.6f} cP")