#%%
# Load libraries
import numpy as np
import matplotlib.pyplot as plt
import gibbs

#%%

T_list = np.linspace(373, 1473, 100)
G_H2O = np.zeros_like(T_list)
G_H2 = np.zeros_like(T_list)
G_Fe2O3 = np.zeros_like(T_list)
G_Fe3O4 = np.zeros_like(T_list)
G_FeO = np.zeros_like(T_list)
G_Fe = np.zeros_like(T_list)
print(gibbs.G_Fe2O3(298))
print(gibbs.G_Fe3O4(298))
print(gibbs.G_FeO(298))
print(gibbs.G_Fe(298))
print(gibbs.G_H2O(298, phase='liquid'))
print(gibbs.G_H2(298))
for i, T in enumerate(T_list):
    G_H2O[i] = gibbs.G_H2O(T, phase='gas')
    G_H2[i] = gibbs.G_H2(T)
    G_Fe2O3[i] = gibbs.G_Fe2O3(T)
    G_Fe3O4[i] = gibbs.G_Fe3O4(T)
    G_FeO[i] = gibbs.G_FeO(T)
    G_Fe[i] = gibbs.G_Fe(T)

fig, ax = plt.subplots()
ax.plot(T_list, G_H2O, label='H2O')
ax.plot(T_list, G_H2, label='H2')
ax.plot(T_list, G_Fe2O3, label='Fe2O3')
ax.plot(T_list, G_Fe3O4, label='Fe3O4')
ax.plot(T_list, G_FeO, label='FeO')
ax.plot(T_list, G_Fe, label='Fe')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Gibbs free energy (kJ/mol)')
ax.legend()
# plt.show()

#%%
# Calculate the equilibrium p-T diagram
# The equilibrium constant K is given by:
# G(T) = -RTln(K)
# K = p_H2O/p_H2

# T_list = np.arange(500, 1500)
R = 8.314/1000 # kJ/(mol K)
# For the reaction: (1) 3Fe2O3(s) + H2(g) = 2Fe3O4(s) + H2O(g)
G_1 = 2*G_Fe3O4 + G_H2O - 3*G_Fe2O3 - G_H2
lnK_1 = -G_1/(R*T_list)
# For the reaction: (2) Fe3O4(s) + H2(g) = 3FeO(s) + H2O(g)
G_2 = 3*G_FeO + G_H2O - G_Fe3O4 - G_H2
lnK_2 = -G_2/(R*T_list)
# For the reaction: (3) FeO(s) + H2(g) = Fe(s) + H2O(g)
G_3 = G_Fe + G_H2O - G_FeO - G_H2
lnK_3 = -G_3/(R*T_list)
# For the reaction: (4) 1/4 Fe3O4(s) + H2(g) = 3/4 Fe(s) + H2O(g)
G_4 = (3*G_Fe + 4*G_H2O - G_Fe3O4 - 4*G_H2)/4
lnK_4 = -G_4/(R*T_list)

fig, ax = plt.subplots()
ax.plot(T_list, -lnK_1, label='K1')
ax.plot(T_list, -lnK_2, label='K2')
ax.plot(T_list, -lnK_3, label='K3')
ax.plot(T_list, -lnK_4, label='K4')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel(r'$ln[p(H_{2})/p(H_{2}O)]$')
ax.legend()

#%%
# Find the intersection of the curves

from scipy.optimize import fsolve

# For the reaction: (5) 4FeO -> Fe + Fe3O4
G_5 = lambda T: gibbs.G_Fe(T) + gibbs.G_Fe3O4(T) - 4*gibbs.G_FeO(T)
res = fsolve(G_5, 500)
Tc = res[0]
print('Equilibrium temperature for reaction 5:', Tc)

T_list_1 = T_list[T_list < Tc]
T_list_2 = T_list[T_list >= Tc]
# lnK_1_Tc = lnK_1[T_list >= Tc]
lnK_2_Tc = lnK_2[T_list >= Tc]
lnK_3_Tc = lnK_3[T_list >= Tc]
lnK_4_Tc = lnK_4[T_list < Tc]

fig, ax = plt.subplots()
ax.plot(T_list, -lnK_1, label='K1')
ax.plot(T_list_2, -lnK_2_Tc, label='K2')
ax.plot(T_list_2, -lnK_3_Tc, label='K3')
ax.plot(T_list_1, -lnK_4_Tc, label='K4')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel(r'$ln[p(H_{2})/p(H_{2}O)]$')
ax.legend()
plt.show()

