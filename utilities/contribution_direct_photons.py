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
#-------------------------------------------------------------------------------------
#                                   Extract Data
#------------------------------------------------------------------------------------

Nph   = input('number of photons ?')
M_ej  = input('mass of the ejecta?')
Phi   = input('opening angle?')
Temp  = input('temperature ?')
kappa = input('opacity?')
f1    = input(' which filter #1 ?')
f2    = input('which filter #2 ?')

print(f'filename is: nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec.txt')
data_folder = Path('outputs_file/')
file        = data_folder / f'nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_spec_delay.txt'
file_direct = data_folder / f'nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_direct_photons.txt'
file_other  = data_folder / f'nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}_other_photons.txt'

file = 'nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_spec_delay.txt'

#file_direct ='nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_direct_photons.txt'
#file_other ='nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_other_photons.txt'

#f1 = 'u'
#f2 = 'z'
#-----------------------------------------------------------------------------------------------------------
nb_obs  = np.genfromtxt(fname=file, max_rows=1, dtype='int')
nb_bins = np.genfromtxt(fname=file, max_rows=1, skip_header=1, dtype='int')
Nstep   = np.genfromtxt(fname=file, max_rows=1, skip_header=2, dtype='int', usecols=(0))
ti      = np.genfromtxt(fname=file,max_rows= 1, skip_header=2, dtype='float', usecols=(1))
tf      = np.genfromtxt(fname=file,max_rows= 1, skip_header=2, dtype='float', usecols=(2))
wave    = np.genfromtxt(fname=file,max_rows=nb_bins, skip_header=3, dtype='float', usecols= (0))
a = int(nb_obs)

#-----------------------------------------------------------------------------------------------------------

step   = (tf-ti)/Nstep
ph     = np.arange(ti+step/2, tf + step/2, step)
index  = np.arange(Nstep)
cost   = np.zeros((nb_obs, Nstep))
if nb_obs!=1:
    for i in range(nb_obs):
        cost[i,:]=1./(nb_obs-1)*i
else:
    cost[0,:]=1


#--------------------------------------------------------------------------------------------------------------
#                   Spectra
#--------------------------------------------------------------------------------------------------------------
intensity =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity[l,i] = list(np.genfromtxt(fname=file, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

intensity_direct =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity_direct[l,i] = list(np.genfromtxt(fname=file_direct, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

intensity_other = np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity_other[l,i] = list(np.genfromtxt(fname=file_other, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mp

#--------------------------------------------------------------------------------------------------------------
#                   Light Curves
#--------------------------------------------------------------------------------------------------------------

#   ************************
#   * Premier filtre choisi*
#   ************************

filter_1 = f'sdss::{f1}'
m_1      = np.zeros((nb_obs, Nstep), dtype=list)

for l in range(nb_obs):
    fluxes_l = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l[i] = (intensity[l,i])

    source_l = sn.TimeSeriesSource(ph,wave,fluxes_l)
    m_1[l,:]   = list(source_l.bandmag(filter_1,"ab",ph))

m_direct_1      = np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    fluxes_l_direct = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l_direct[i] = (intensity_direct[l,i])

    source_l_direct = sn.TimeSeriesSource(ph,wave,fluxes_l_direct)
    m_direct_1[l,:]   = list(source_l_direct.bandmag(filter_1,"ab",ph))


m_other_1      = np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    fluxes_l_other = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l_other[i] = (intensity_other[l,i])
    # print(f'observers number {l}', fluxes_l)
    source_l_other   = sn.TimeSeriesSource(ph,wave,fluxes_l_other)
    m_other_1[l,:]   = list(source_l_other.bandmag(filter_1,"ab",ph))

#   *************************
#   * Deuxième filtre choisi*
#   *************************
filter_2 = f'sdss::{f2}'
m_2      = np.zeros((nb_obs, Nstep), dtype=list)

for l in range(nb_obs):
    fluxes_l = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l[i] = (intensity[l,i])
   
    source_l = sn.TimeSeriesSource(ph,wave,fluxes_l)
    m_2[l,:]   = list(source_l.bandmag(filter_2,"ab",ph))

m_direct_2      = np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    fluxes_l_direct = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l_direct[i] = (intensity_direct[l,i])
 
    source_l_direct = sn.TimeSeriesSource(ph,wave,fluxes_l_direct)
    m_direct_2[l,:] = list(source_l_direct.bandmag(filter_2,"ab",ph))


m_other_2      = np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    fluxes_l_other = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_l_other[i] = (intensity_other[l,i])

    source_l_other = sn.TimeSeriesSource(ph,wave,fluxes_l_other)
    m_other_2[l,:]   = list(source_l_other.bandmag(filter_2,"ab",ph))      # génère une ndarray


#---------------------------------------------------------------------------------------------------------------
#               PLOT
#--------------------------------------------------------------------------------------------------------------


fig = plt.figure()
plt.subplots_adjust(wspace=0.25, hspace=0.2) #wspace = écartement horizontal, #vspace=vertical ?

for l in range(nb_obs):
    
    ratio_blue = 10**(-m_direct_1[l]/(2.5))/10**(-m_other_1[l]/(2.5))
    ratio_red  = 10**(m_direct_2[l]/(-2.5))/10**(m_other_2[l]/(-2.5))
    
    
    y = np.full(np.shape(ratio_red),1)
    numPlots = 211 + l
    
    Ax, Lw = asManyPlots(numPlots, [ph,ph, ph], [ratio_blue, ratio_red, y], color=['teal', 'indianred', 'sandybrown'], markerSize=[0,0,0], linestyle =['-','-','dashed'], plotFlag =[True, True, True],
                         showLegend = True,
                         label =[f'in sdss:: {f1}', f'in sdss:: {f2}', 'ratio equal to 1'], ylabel =r'$F_{\rm{direct}}/F_{\rm{nondirect}}$',
                         linewidth = 1,
                         legendTextSize = 8, textsize=8, tickSize = 10 )
    Ax.yaxis.set_minor_locator(MultipleLocator(1))
    Ax.xaxis.set_minor_locator(MultipleLocator(1))
    if l == 0:
        Ax.text(1, 30,'edge on', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=7)

    else:
        Ax.text(1, 7,'face on', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=7)
Ax.set_xlabel('Time since merger (in days)')
fig.show()
output1 = f'ratio_fluxs_nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}.pdf'

plt.savefig(output1 , bbox_inches='tight')

#--------------------------------------------------------------------------

#fig, ax = plt.subplots(2,2, figsize=(10,4))
f = plt.figure()
plt.subplots_adjust(wspace=0.25, hspace=0.2) #wspace = écartement horizontal, #vspace=vertical ?

#plt.subplots(221)

dataY_1 = np.zeros((nb_obs, 3), dtype=list)
dataY_2 = np.zeros((nb_obs, 3), dtype=list)

for l in range(nb_obs):
    dataY_1[l][0] = list(m_1[l])
    dataY_1[l][1] = list(m_direct_1[l])
    dataY_1[l][2] = list(m_other_1[l])

    dataY_2[l][0] = list(m_2[l])
    dataY_2[l][1] = list(m_direct_2[l])
    dataY_2[l][2] = list(m_other_2[l])


Ax1, LW = asManyPlots(221, [ph, ph, ph], [dataY_1[0][0], dataY_1[0][1], dataY_1[0][2]], ylabel=f'$m_{f1}$',
                         markerSize = [0,0,0], color=['crimson', 'coral', 'darkseagreen'], linestyle =['-','-', '-'], linewidth=1, textsize= 8, legendTextSize=6,
                      plotFlag=[True, True, True], showLegend=[True], label=[ 'total fluxes',  'direct photons','non direct photons'], tickSize = 10)

Ax1.text(1.5, 33,f'filter: {f1} \n edge on ', bbox=dict(facecolor='skyblue', alpha=0.2, boxstyle='round'), size=7)
plt.gca().invert_yaxis()

    
Ax2, Lw = asManyPlots(222, [ph, ph, ph], [dataY_2[0][0], dataY_2[0][1], dataY_2[0][2]], ylabel=f'$m_{f2}$',
            markerSize = [0,0,0], color=['crimson', 'coral', 'darkseagreen'], linestyle =['-','-', '-'], textsize= 8, linewidth=1, legendTextSize=6,
            plotFlag=[True, True, True], showLegend=False,
            tickSize= 10)

Ax2.text(2, 24,f'filter : {f2} \n edge on', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=7)
plt.gca().invert_yaxis()
Ax3, Lw = asManyPlots(223, [ph, ph, ph], [dataY_1[1][0], dataY_1[1][1], dataY_1[1][2]], ylabel=f'$m_{f1}$',
                                markerSize = [0,0,0], color=['crimson', 'coral', 'darkseagreen'], linestyle =['-','-', '-'], linewidth=1, textsize= 8, legendTextSize=6,
                                plotFlag=[True, True, True], showLegend=False, tickSize=10)

Ax3.text(1.5, 33,f'filter: {f1} \n face on ', bbox=dict(facecolor='skyblue', alpha=0.2, boxstyle='round'), size=7)
plt.gca().invert_yaxis()
Ax4, Lw = asManyPlots(224, [ph, ph, ph], [dataY_2[1][0], dataY_2[1][1], dataY_2[1][2]], xlabel='Time since merger (in days)', ylabel=f'$m_{f2}$',
                      markerSize = [0,0,0], color=['crimson', 'coral', 'darkseagreen'], linestyle =['-','-', '-'], linewidth=1, textsize= 10, legendTextSize=6,
                      plotFlag=[True, True, True], showLegend=False, tickSize = 10)
Ax4.xaxis.set_label_coords(-0.25,-0.15)
Ax4.text(2, 24,f'filter : {f2} \n face on', bbox=dict(facecolor='salmon', alpha=0.2, boxstyle='round'), size=7)

plt.gca().invert_yaxis()            # reverse y axis
Ax1.xaxis.set_minor_locator(MultipleLocator(1))
Ax1.yaxis.set_minor_locator(MultipleLocator(1))
Ax2.xaxis.set_minor_locator(MultipleLocator(1))
Ax2.yaxis.set_minor_locator(MultipleLocator(1))
Ax3.xaxis.set_minor_locator(MultipleLocator(1))
Ax3.yaxis.set_minor_locator(MultipleLocator(1))
Ax4.xaxis.set_minor_locator(MultipleLocator(1))
Ax4.yaxis.set_minor_locator(MultipleLocator(1))


f. show()
output2 = f'contribution_of_photons_nph{Nph}_mej{M_ej}_phi{Phi}_T{Temp}_aopac{kappa}.pdf'
plt.savefig(output2, bbox_inches='tight')



#-------------------------------------------------------
