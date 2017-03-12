import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches


LL_dielectric = []
inverse_LL_dielectric = []
LL_G = []

kerker_dielectric = []

k=10000

#fermi_wavevector = 1.0
a0 = 1.0
delta = 0.08

kerker_a = 0.8
kerker_b = 1.5

#LL_prefactor = 2.0 / (fermi_wavevector * a0 * 3.14159265)


j=0

LL_G.append(10)
LL_dielectric.append(1)
kerker_dielectric.append(1)



inverse_LL_dielectric.append(1/LL_dielectric[j])

for j in range(1,k):

    LL_G.append( LL_G[j-1] - 10/float(k) )

    fermi_wavevector = 0.98*(0.2*(LL_G[j]**2) + 2.25)
    LL_prefactor = 2.0 / (fermi_wavevector * a0 * 3.14159265)

    LL_dielectric.append(1 + LL_prefactor * ( 1/((LL_G[j]/fermi_wavevector)**2) - (delta/(2*((LL_G[j]/fermi_wavevector)**3)))*
                        ( np.arctan( (2*(LL_G[j]/fermi_wavevector) - ((LL_G[j]/fermi_wavevector))**2) / delta )
                        + np.arctan( (2*(LL_G[j]/fermi_wavevector) + ((LL_G[j]/fermi_wavevector))**2) / delta )  ) +
                        ( delta**2/(8*(LL_G[j]/fermi_wavevector)**5) +
                        1 /(2.0*(LL_G[j]/fermi_wavevector)**3) - 1/(8*(LL_G[j]/fermi_wavevector)) ) *
                        ( np.log( delta**2 + (2*(LL_G[j]/fermi_wavevector) + (LL_G[j]/fermi_wavevector)**2)**2 )
                        - np.log( delta**2 + (2*(LL_G[j]/fermi_wavevector) - (LL_G[j]/fermi_wavevector)**2)**2 ) ) ))

    inverse_LL_dielectric.append(1/LL_dielectric[j])


    if LL_G[j]>8.0:

        kerker_dielectric.append(1.0)
        
    else:
    
        kerker_dielectric.append((kerker_a * LL_G[j]**2)/(LL_G[j]**2 + kerker_b**2))


plt.ylim(0,1)
plt.xlim(0,10)

LL = plt.plot(LL_G,inverse_LL_dielectric,color="blue", label = "Levine-Louie")
kerker = plt.plot(LL_G,kerker_dielectric,color="red", label = "Kerker")

#first_legend = plt.legend([LL, kerker], loc=1)
plt.legend(loc=4)
plt.xlabel('$|G|$ $(\AA^{-1})$')
plt.ylabel('Inverse dielectric $(\epsilon{^-1})$')
plt.title('Levine-Louie (d=0.08 qf=$0.98 ( 0.2G^2 + 2.25 )$) vs Kerker (a=0.8, b=1.5) Inverse Dielectric', fontsize=12)
plt.savefig('k_default_LL_0.08_match_dielectric.pdf')
#plt.show()

