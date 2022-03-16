# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10  19:24 2020
@author: Lina ISSA
"""
import numpy as np
import matplotlib.pyplot as plt
from math import asin,pi, sqrt
from plotUtilities import *
import matplotlib.patches as patches
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib.colors import (Normalize,ListedColormap, BoundaryNorm)

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
phi       = pi/4   # ouverture angulaire de 45° % a l'axe de la fusion
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
Mej_dyn  = []
Mej_wind = []
with open("grid_input_opang45.txt", "w+", encoding = 'utf-8') as model:
    model.write(f"{nb_cell} {R:e} {t} \n")
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
                    #Mej_wind.append(dens*step**3)
                    else:                         #in the dynamical outflow
                        dens = density(r,alpha2)
                        theta = asin((n-step/2)/r)
                        #Mej_dyn.append((dens*step**3))
                        if abs(theta) < phi:     # red ejecta
                            gamma = 0.2
                        else:
                            gamma = 0.4
        Ye_tot.append(gamma)
            dens_tot.append(dens)
            X.append(l-step/2)
            Y.append(m-step/2)
            Z.append(n-step/2)
            model.write(f'{(l-step/2):e} {(m-step/2):e} {(n-step/2):e} {dens:e} {gamma}\n')

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
        if np.shape(np.where(whereY_and_Z)[0])[0]!=0:
            #for i in range(np.shape(np.where(whereY_and_Z)[0])[0]):
            position = np.array(np.where(whereY_and_Z))[0][0]
            valGamma = Ye_tot[position]
            valRho   = dens_tot[position]
            Gamma[index1,index2] = valGamma
            Rho[index1,index2]   = valRho


v_y = Y_grid/(t*clight)
v_z = Z_grid *1/t *1/clight

v_ejecta_min = vmin *1/clight
v_wind_max   = vmax_wind * 1/clight
v_ejecta_max = 0.8


x = plt.subplot(111)
#ax.yaxis.set_ticks_position('both')
#ax.xaxis.set_ticks_position('both')
#ax.tick_params(which='both', direction='in', labelsize=16)
plt.xticks(np.arange(-0.8,0.9,0.2))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(axis='both', labelsize = 10)

a1 = patches.Ellipse((0,0), v_ejecta_min*2,v_ejecta_min*2, edgecolor="bisque", facecolor='None',linewidth=1,alpha=1)#,label="Lanthanide-free")
plt.gca().add_patch(a1)
a2 = patches.Ellipse((0,0), v_wind_max*2,v_wind_max*2, edgecolor="bisque", facecolor='None',linewidth=10,alpha=1,label="Wind outflow")
plt.gca().add_patch(a2)

#plt.fill_between(a,b)
a3 = patches.Ellipse((0,0), v_ejecta_max*2,v_ejecta_max*2, edgecolor="bisque", facecolor='None',linewidth=1,alpha=1 ,label="Dynamical outflow")
plt.gca().add_patch(a3)

dt = ax.contourf(v_y, v_z, np.log10(Rho), 50, cmap='plasma')# vmin = -46, vmax = -45)#vmin=np.log10(min(dens_tot)), vmax=np.log10(max(dens_tot)))

col = plt.colorbar(dt)
col.set_label(r'$\rho$', size=16)
ax.set_xlabel(r'$v_{y} \,( in \, \rm{c})$', size=16)
ax.set_ylabel(r'$v_{z} \,( in \, \rm{c})$', size=16)
col.ax.tick_params(labelsize=12)
ax.set_tilte('Density mapp at t= 1 day')
#cont = plt.contour(Y_grid, Z_grid, Rho) #,levels=[0.3,0.5,1], colors='black', linestyle="--")
#cont.clabel(fontsize=10)
plt.savefig("mapping_density_test.pdf", bbox_inches='tight')
