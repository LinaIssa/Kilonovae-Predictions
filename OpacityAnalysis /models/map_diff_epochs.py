# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import asin,pi, sqrt
import matplotlib.patches as patches
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator, LogLocator)

clight    = 2.99792458e10 # en cm
#Msun      = 1.989e33

## PARAMETERS ##
vmax      = 0.9*clight
vmin_dyn  = 0.08*clight
vmax_wind = 0.08*clight
vmin      = 0.025*clight
day       = 3600*24 # time equal to one day
nb_cell   = 50 # numbers of cells in one direction
phi       = pi/4   # ouverture angulaire de 45° % a l'axe de la fusion
alpha1    = 3
alpha2    = 6



def density(r,alpha):
    return  A*r**(-alpha)     # on prend A = 1

EPOCH = [1, 4, 6, 8, 10, 12]
for (i, epoch)  in enumerate(EPOCH):
    t    = epoch * day
    R    = t*vmax     #le rayon de la grille
    step = 2*R/(nb_cell)
    x = np.zeros(nb_cell)
    for a in range(nb_cell):
        x[a] = R - a*step
    
    y = np.zeros(nb_cell)
    for a in range(nb_cell):
        y[a] = R - a*step
    
    z = np.zeros(nb_cell)
    for a in range(nb_cell):
        z[a] = R - a*step

    X = []
    Y = []
    Z = []
    Ye_tot = []
    dens_tot = []
    dataR = []
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
                        theta = asin((n-step/2)/r)
                        dens = density(r,alpha2) * np.sin(np.pi/2-theta)**2
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
                dataR.append(r)

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
            position_list = np.array(np.where(whereY_and_Z))[0] # liste des nb_cell positions où le couple (Y,Z) est dégénré pour chaque valeur de X
            index = int((nb_cell/2)-1) # choix de index guidé par la valeur souhaitée à X
            position = position_list[index]
            #print( X[position])
            valGamma = Ye_tot[position]
            valRho   = dens_tot[position]
            Gamma[index1,index2] = valGamma
            Rho[index1,index2]   = valRho

    v_y = Y_grid/(t*clight)
    v_z = Z_grid *1/t *1/clight
    v_ejecta_min = vmin *1/clight
    v_wind_max   = vmin_dyn * 1/clight
    v_ejecta_max = vmax /clight
    
    numplot      = 231 + i
    
    ax = plt.subplot(numplot)
    
    plt.xticks(np.arange(-0.8,0.9,0.4))
    plt.yticks(np.arange(-0.8,0.9,0.4))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.tick_params(axis='both', labelsize = 7)
    a2 = patches.Ellipse((0,0), v_wind_max*2,v_wind_max*2, edgecolor="yellow", facecolor='yellow', fill = True,linewidth=3,alpha=1,label="wind outflow")
    plt.gca().add_patch(a2)
    a3 = patches.Ellipse((0,0), v_ejecta_max*2,v_ejecta_max*2, edgecolor="powderblue", facecolor='None', linestyle = '--', linewidth=1,label=" border line of the \n dynamical ejecta")
    plt.gca().add_patch(a3)
    
    dt = ax.contourf(v_y, v_z,np.log10(Rho), 100, cmap='plasma', vmin = -40, vmax=0)
    
    
    #cont = plt.contour(v_y, v_z, np.log10(Rho), levels = [-42, -40], colors='black')
    #cont.clabel(fontsize=7)
    
    if i == 0 or i == 3:
        ax.set_ylabel(r'$v_{z} \,( in \, \rm{c})$', size=15)
    col = plt.colorbar(dt) #,cax = plt.axes([-40, -35, -30, -25, -20, -15, -10,  -5,   0]))
    col.ax.tick_params(labelsize=10)
    if i == len(EPOCH) -1:
        col.set_label(r'$\log(\rho)$ (not normalized)', size=15)#, x = 3, y = 1.5)

    ax.text(0.35, 0.95, 'Day {}'.format(epoch), bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=8)


ax.xaxis.set_label_coords(-1,-0.2)
ax.set_xlabel(r'$v_{y} \,( in \, \rm{c})$', size=15)
plt.subplots_adjust(left = -0.1, right = 1.2, bottom=0,top=1, wspace = 0.3, hspace= 0.2)

plt.legend( loc= 4, bbox_to_anchor=(1.2,-0.35), fontsize= 8,framealpha = 1)#,  shadow = True)
plt.savefig("mapping_density_test.pdf", bbox_inches='tight')
