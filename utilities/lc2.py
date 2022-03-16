# -*- coding: utf-8 -*-
"""
Created on Dec the 27 2019
@author: Lina ISSA
"""
import numpy as np
import sncosmo as sn
import matplotlib.pyplot as plt
import matplotlib
import glob
from plotUtility import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib.colors import (ListedColormap, BoundaryNorm)
from matplotlib import cm
import sys
from palettable.cartocolors.diverging import ArmyRose_7
from palettable.scientific.sequential import Acton_11
from pathlib import Path
#-------------------------------------------------------------------------------------------------------------
#                   Extract Data
#--------------------------------------------------------------------------------------------------------------
Nph   = input('number of photons ?')
M_ej  = input('mass of the ejecta?')
Phi   = input('opening angle?')
Temp  = input('temperature ?')
kappa = input('opacity?')

print(f'filename is: nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec.txt')

file = Path('outputs_file/') / f'nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec.txt'
#file   = 'nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_spec_ref.txt'
#-------------------------------------------------------
nb_obs  = np.genfromtxt(fname=file, max_rows=1, dtype='int')
nb_bins = np.genfromtxt(fname=file, max_rows=1, skip_header=1, dtype='int')
Nstep   = np.genfromtxt(fname=file, max_rows=1, skip_header=2, dtype='int', usecols=(0))
ti      = np.genfromtxt(fname=file,max_rows= 1, skip_header=2, dtype='float', usecols=(1))
tf      = np.genfromtxt(fname=file,max_rows= 1, skip_header=2, dtype='float', usecols=(2))
wave    = np.genfromtxt(fname=file,max_rows=nb_bins, skip_header=3, dtype='float', usecols= (0))
a = int(nb_obs)
#print(nb_obs, '\n', nb_bins, '\n', wave)
#-------------------------------------------------------

step   = (tf-ti)/Nstep
ph     = np.arange(ti+step/2, tf + step/2, step)
index  = np.arange(Nstep)
cost   = np.zeros((nb_obs, Nstep))
if nb_obs!=1:
    for i in range(nb_obs):
        cost[i,:]=1./(nb_obs-1)*i
else:
    cost[0,:]=1
#print(step, '\n', ph, '\n', index)
#-------------------------------------------------------------------------------------------------------------
#                   Spectra
#--------------------------------------------------------------------------------------------------------------
intensity =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity[l,i] = list(np.genfromtxt(fname=file, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

#--------------------------------------------------------------------------------------------------------------
#                   Light Curves
#--------------------------------------------------------------------------------------------------------------
filter = 'sdss::u'
#source = np.zeros((nb_obs,Nstep))
m      = np.zeros((nb_obs, Nstep), dtype=list)

for l in range(nb_obs):
    fluxes_l = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l[i] = (intensity[l,i])
# print(f'observers number {l}', fluxes_l)
    source_l = sn.TimeSeriesSource(ph,wave,fluxes_l)
    m[l,:]   = list(source_l.bandmag(filter,"ab",ph))      # génère une ndarray

#---------------------------------------------------------------------------------------------------------------
#               PLOT
#--------------------------------------------------------------------------------------------------------------


plt.rcParams['figure.figsize'] = [11, 11] #[largeur, hauteur]
f = plt.figure()
plt.subplots_adjust(wspace=0.05, hspace=0.05) #wspace = écartement horizontal, #vspace=vertical ?

X = [list(ph)]*a
dataX = list(np.copy(X))
dataX.append(X)

Y = [m[l] for l in range(nb_obs)]
dataY = list(np.copy(Y))
dataY.append(np.ravel(m))

dataX = [ph, ph, [ph,ph]]
dataY = [m[0], m[1], np.ravel(m)]
dataY= [m[l] for l in range(nb_obs)]


def selectRGB(cmap,N):
    nb = cmap.N//(N-1)
    ls =[]
    for i in range(0,cmap.N, nb-1):
        ls.append(cmap.colors[i])
    return ls

map = matplotlib.cm.get_cmap('plasma')
values = selectRGB(map, 11)
cmap = ListedColormap(values)
norm = BoundaryNorm(np.arange(0,1,0.1),cmap.N)

Ax, LW= asManyPlots(111, dataX, dataY, xlabel='time since merger (in days)', ylabel='$m_{g}$',
                    markerSize=np.arange(0,nb_obs+1)*0,
                    color = values +[cost],
                    linestyle=['-' for l in range(nb_obs)]+ [None], linewidth= 1, textsize=18,
                    cmap=cmap, norm=norm,
                    plotFlag=[True for l in range(nb_obs)]+[False], colorbarLabel = r'$cos\theta_{\rm{obs}}$',
                    colorbarLabelSize = 18,
                    showColorbar=True)



plt.gca().invert_yaxis()            # reverse y axis
plt.xticks(np.arange(0,15,2))
Ax.xaxis.set_minor_locator(MultipleLocator(1))
Ax.yaxis.set_minor_locator(MultipleLocator(1))
f.show()
output =f'lc_nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_{nb_obs}obs.pdf'
plt.savefig(output, bbox_inches='tight')

#cl = ['red','blue', 'cyan', 'green', 'purple', 'orange', 'olive', 'black', 'pink','brown', 'yellow']
#cl = ['blue', 'green']
#cl = Acton_11
#mamap = get_mpl_colormap(ArmyRose_7)
#cmap = ListedColormap(cl)
#tu peux essayer cmap = ListedColormap(plasma.colors, N=2)
#bounds = np.arange(0,1, 0.1)

#norm   = BoundaryNorm(bounds,cmap.N)


#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




#plt.rcParams['figure.figsize'] = [11, 11] #[largeur, hauteur]
#plt.subplots_adjust(wspace=0.05, hspace=0.05) #wspace = écartement horizontal, #vspace=vertical ?

#time = np.zeros((nb_obs,Nstep), dtype=list)
#for i in range(nb_obs):
#    time[i,:] = ph

#time = []
#for l in range(nb_obs):
#    for t in list(ph):
#        time.append(t)


#magn = []
    #for p in m:
    #   for q in p:
#    magn.append(q)

#cost_Z = []
    #for a in cost:
    #    for b in a:
#    cost_Z.append(b)

#X,Y = np.meshgrid(np.unique(time), np.unique(magn))
#Z   = np.full(np.shape(X), 0.5)
#print(np.shape(X))

#size1 = np.shape(X)[0]
#size2 = np.shape(Y)[0]

#for index1, x, y in zip(range(size1), X, Y):
#    for index2, valx, valy in zip(range(size1), x, y):
#        whereX = time == valx
#whereY = magn == valy
        #        whereX_and_Y = np.logical_and(whereX,whereY)
        #        if np.shape(np.array(np.where(whereX_and_Y))[0])[0]!=0:
        #    position = np.array(np.where(whereX_and_Y))[0][0]
        #    valZ = cost_Z[position]
#    Z[index1, index2] = valZ

#ax = plt.subplot(111)
#ax.yaxis.set_ticks_position('both')
#ax.xaxis.set_ticks_position('both')
#ax.tick_params(which='both', direction='in', labelsize=26)
#dt = ax.contourf(X, Y, Z, cmap='plasma')
#col = plt.colorbar(dt)
#col.set_label(r'$|\Delta \log P| \, (\rm{dyn/cm^{2}})$', size=26)
#ax.set_xlabel(r'$\log \rho \, \rm{(g/cm^{3})}$', size=26)
#ax.set_ylabel(r'$\log T \, \rm{(K)}$', size=26)
#col.ax.tick_params(labelsize=22)
#cont = ax.contour(X,Y, Z)
#cont.clabel(fontsize=10)
#plt.gca().invert_yaxis()            # reverse y axis

#plt.show()


