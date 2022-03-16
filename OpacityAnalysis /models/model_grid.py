
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import asin,pi, sqrt
import matplotlib.patches as patches
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

clight    = 2.99792458e10 # en cm
#Msun      = 1.989e33

## PARAMETERS ##
vmax      = 0.9*clight
vmin_dyn  = 0.08*clight
vmax_wind = 0.08*clight
vmin      = 0.025*clight
t         = 3600*24 # time equal to one day
nb_cell   = 100 # numbers of cells in one direction
R         = t*vmax     #le rayon de la grille
step      = 2*R/(nb_cell)
phi       = np.pi/4   # ouverture angulaire de 45° % a l'axe de la fusion
alpha1    = 3
alpha2    = 6



x = np.zeros(nb_cell)
for a in range(nb_cell):
    x[a] = R - a*step

y = np.zeros(nb_cell)
for a in range(nb_cell):
    y[a] = R - a*step

z = np.zeros(nb_cell)
for a in range(nb_cell):
    z[a] = R - a*step



def density(r,alpha):
    return  A*r**(-alpha)     # on prend A = 1

X = []
Y = []
Z = []
Ye_tot = []
dens_tot = []
for l in x:
    for m in y:
        for n in z:
            r = sqrt((l-step/2)**2 + (m-step/2)**2 + (n-step/2)**2)
            if r>R or r<(vmin*t):               #outside the sphere
                A = 1.
                dens  = density(r,alpha1)
                gamma = 0.4
            else:
                A = 10**50
                if (vmin*t)<r<(vmin_dyn*t):  #in the wind outflow
                    dens  = density(r,alpha1)
                    gamma = 0.25
                else:                         #in the dynamical outflow
                    dens = density(r,alpha2)
                    theta = asin((n-step/2)/r)
                    if abs(theta) < phi:     # red ejecta
                        gamma = 0.2
                    else:
                        gamma = 0.4
            Ye_tot.append(gamma)
            dens_tot.append(dens)
            X.append(l-step/2)
            Y.append(m-step/2)
            Z.append(n-step/2)

Y_grid, Z_grid = np.meshgrid(np.sort(np.unique(Y)),np.sort(np.unique((Z))))

size1 = np.shape(Y_grid)[0]
size2 = np.shape(Z_grid)[0]

Gamma  = np.full(np.shape(Y_grid), 0.5)
Rho    = np.full(np.shape(Y_grid), 0.5)
for index1, yl, zl in zip(range(size1), Y_grid, Z_grid):
    for index2, valY, valZ in zip(range(size2), yl, zl):
        whereY       = Y == valY
        whereZ       = Z == valZ
        whereY_and_Z = np.logical_and(whereY,whereZ)
        position_list = np.array(np.where(whereY_and_Z))[0]
        index = [(nb_cell/2)-1]
        position = np.array(np.where(whereY_and_Z))[index]
        valGamma = Ye_tot[position]
        valRho   = dens_tot[position]
        Gamma[index1,index2] = valGamma
        Rho[index1,index2]   = valRho


v_y = Y_grid/(t*clight)
v_z = Z_grid *1/t *1/clight

v_ejecta_min = vmin *1/clight
v_wind_max   = vmax_wind * 1/clight
v_ejecta_max = vmax/clight

##### PLOT 1 : MAPPING the DENSITY ######
ax = plt.subplot(111)
plt.xticks(np.arange(-0.8,0.9,0.2))
plt.yticks(np.arange(-0.8,0.9,0.2))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(axis='both', labelsize = 10)

#a1 = patches.Ellipse((0,0), v_ejecta_min*2,v_ejecta_min*2, edgecolor="bisque", facecolor='None',linewidth=1,alpha=1)#,label="Lanthanide-free")
#plt.gca().add_patch(a1)
a2 = patches.Ellipse((0,0), v_wind_max*2,v_wind_max*2, edgecolor="yellow", facecolor='yellow', fill = True, linewidth=3,alpha=1,label="wind outflow")
plt.gca().add_patch(a2)

#plt.fill_between(a,b)
a3 = patches.Ellipse((0,0), v_ejecta_max*2,v_ejecta_max*2, edgecolor="powderblue", facecolor='None',linestyle = '--', linewidth=2,alpha=7 ,label=" border line of the \n dynamical outflow")
plt.gca().add_patch(a3)

dt = ax.contourf(v_y, v_z, np.log10(Rho), 100, cmap='plasma')# vmin = -46, vmax = -45)#vmin=np.log10(min(dens_tot)), vmax=np.log10(max(dens_tot)))

col = plt.colorbar(dt)
col.set_label(r'$\rho$ (not normalized)', size=16)
ax.set_xlabel(r'$v_{y} \,( in \, \rm{c})$', size=16)
ax.set_ylabel(r'$v_{z} \,( in \, \rm{c})$', size=16)
col.ax.tick_params(labelsize=12)
ax.set_title('Density map at 1 day')
plt.legend(loc=1, fontsize= 8,framealpha = 1)
plt.savefig("mapping_density_test.pdf", bbox_inches='tight')




##### PLOT 2 : MAPPING the DENSITY ######
ax2 = plt.subplot(111)
#ax.yaxis.set_ticks_position('both')
#ax.xaxis.set_ticks_position('both')
#ax.tick_params(which='both', direction='in', labelsize=16)
plt.xticks(np.arange(-0.8,0.9,0.2))
ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
ax2.tick_params(axis='both', labelsize = 10)

a12 = patches.Ellipse((0,0), v_ejecta_min*2,v_ejecta_min*2, edgecolor="bisque", facecolor='None',linewidth=1,alpha=1)#,label="Lanthanide-free")
plt.gca().add_patch(a2)
a22 = patches.Ellipse((0,0), v_wind_max*2,v_wind_max*2, edgecolor="bisque", facecolor='None',linewidth=10,alpha=1,label="Wind outflow")
plt.gca().add_patch(a22)

#plt.fill_between(a,b)
a32 = patches.Ellipse((0,0), v_ejecta_max*2,v_ejecta_max*2, edgecolor="bisque", facecolor='None',linewidth=1,alpha=1 ,label="Dynamical outflow")
plt.gca().add_patch(a32)

dt = ax.contourf(v_y, v_z, Gamma, 100, cmap='plasma')# vmin = -46, vmax = -45)#vmin=np.log10(min(dens_tot)), vmax=np.log10(max(dens_tot)))

col = plt.colorbar(dt)
col.set_label(r'$\rho$', size=16)
ax.set_xlabel(r'$v_{y} \,( in \, \rm{c})$', size=16)
ax.set_ylabel(r'$v_{z} \,( in \, \rm{c})$', size=16)
col.ax.tick_params(labelsize=12)
ax.set_title('Gamma map at 1 day')
plt.savefig("mapping_density_test.pdf", bbox_inches='tight')
