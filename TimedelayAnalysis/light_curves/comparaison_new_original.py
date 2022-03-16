import numpy as np
import sncosmo as sn
import matplotlib.pyplot as plt
import matplotlib
import glob
from plotUtilities import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib.colors import (ListedColormap, BoundaryNorm)
from matplotlib import cm
import sys
from palettable.cartocolors.diverging import ArmyRose_7
from palettable.scientific.sequential import Acton_11
import os
from pathlib import Path
#-------------------------------------------------------------------------------------------------------------
#                   Extract Data
#--------------------------------------------------------------------------------------------------------------
Nph    = input('number of photons ?')
M_ej   = input('mass of the ejecta?')
Phi    = input('opening angle?')
Temp   = input('temperature ?')
kappa  = input('opacity?')

print(f'filename is: nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec.txt')
data_folder = Path('outputs_file/')
file        = data_folder / f'nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec_ref.txt'
file_delay  = data_folder / f'nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec_delay.txt'

print(file, file_delay)
#file ='nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_spec_ref.txt'
#file_delay = 'nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_spec_delay.txt'
#-------------------------------------------------------
nb_obs  = np.genfromtxt(fname=file, max_rows=1, dtype='int')
nb_bins = np.genfromtxt(fname=file, max_rows=1, skip_header=1, dtype='int')
Nstep   = np.genfromtxt(fname=file, max_rows=1, skip_header=2, dtype='int', usecols=(0))
ti      = np.genfromtxt(fname=file,max_rows= 1, skip_header=2, dtype='float', usecols=(1))
tf      = np.genfromtxt(fname=file,max_rows= 1, skip_header=2, dtype='float', usecols=(2))
wave    = np.genfromtxt(fname=file,max_rows=nb_bins, skip_header=3, dtype='float', usecols= (0))
a = int(nb_obs)
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

#-------------------------------------------------------------------------------------------------------------
#                   Spectra
#--------------------------------------------------------------------------------------------------------------
intensity =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity[l,i] = list(np.genfromtxt(fname=file, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

intensity_delay =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity_delay[l,i] = list(np.genfromtxt(fname=file_delay, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

#--------------------------------------------------------------------------------------------------------------
#                   Light Curves
#--------------------------------------------------------------------------------------------------------------

filter = ['sdss::u','sdss::g','sdss::i', 'swope2::h']
m = []
m_delay = []
for (f, band) in enumerate(filter):

    m_f       = np.zeros((nb_obs, Nstep), dtype=list)
    m_delay_f = np.zeros((nb_obs, Nstep), dtype=list)

    for l in range(nb_obs):
        fluxes_l       = np.zeros((Nstep, nb_bins),dtype=list)
        fluxes_delay_l = np.zeros((Nstep, nb_bins),dtype=list)
        for i in index:
            fluxes_l[i] = (intensity[l,i])
            source_l = sn.TimeSeriesSource(ph,wave,fluxes_l)
            fluxes_delay_l[i] = (intensity_delay[l,i])
            source_delay_l = sn.TimeSeriesSource(ph,wave,fluxes_delay_l)
        m_f[l,:]       = list(source_l.bandmag(band,"ab",ph))
        m_delay_f[l,:] =  list(source_delay_l.bandmag(band,"ab",ph))
    m.append(m_f)
    m_delay.append(m_delay_f)


#---------------------------------------------------------------------------------------------------------------
#               PLOT
#--------------------------------------------------------------------------------------------------------------

fig = plt.figure()
plt.subplots_adjust(wspace=0.35, hspace=0.4) #wspace = écartement horizontal, #vspace=vertical ?

for l in range(nb_obs):
    if l == 0:
        obser = 'edge on'
    else :
        obser = 'face on'
    
    for (band,f, i) in zip(filter, [0,1,2,3], range(l*4,l*4+4,1)):
        dataY       = list(m[f][l])
        dataY_delay = list(m_delay[f][l])
        numPlot = 241 + i
        if i == 0 or i == 4:
            Ax, LW = asManyPlots(numPlot, [ph, ph], [dataY, dataY_delay], ylabel = 'Magnitude',
                                 markerSize = [0,0], color=['crimson', 'sandybrown'], linestyle =['-','-'], linewidth=1, textsize= 8, legendTextSize=6,
                                 plotFlag=[True, True], showLegend=False, tickSize = 7)
        else:
            Ax, LW = asManyPlots(numPlot, [ph, ph], [dataY, dataY_delay],
                                 markerSize = [0,0], color=['crimson', 'sandybrown'], linestyle =['-','-'], linewidth=1, textsize= 8, legendTextSize=7,
                                 plotFlag=[True, True], showLegend=False, tickSize = 7)

        if band == 'sdss::u':
            Ax.text(1.5, 30.5 ,f'filter: {band}', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=6)
        if band == 'sdss::g':
            Ax.text(4.5, 20 ,f'filter: {band}', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=6)
        if i == 2:
            Ax.text(4.5, 18.5,f'filter: {band}', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=6)
        if i == 3:
            Ax.text(3, 19.3,f'filter: {band}', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=6)
        if i == 6:
            Ax.text(5, 18,f'filter: {band}', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=6)
        if i == 7:
            Ax.text(4, 17.8,f'filter: {band}', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=6)

        plt.gca().invert_yaxis()
        Ax.xaxis.set_minor_locator(MultipleLocator(1))
        Ax.yaxis.set_minor_locator(MultipleLocator(1))
        if i == 0:
            Ax.legend(['original code', 'time delay'], fontsize = 5)
            Ax.set_title(f'Lights curves drawn {obser}', x = 2.5, y =1.15, fontsize = 10)
        if i == 4:
            Ax.set_xlabel('time since merger (in days)', fontsize = 8)
            Ax.xaxis.set_label_coords(2.5,-0.15)

fig.suptitle(f'Lights curves drawn {obser}', x =0.5, y =0.5, fontsize = 10)
fig.show()
output = f'discrepancies_in_different_band{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}.pdf'

plt.savefig(output, bbox_inches='tight')
