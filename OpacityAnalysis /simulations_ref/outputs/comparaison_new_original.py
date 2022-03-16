
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
Nph      = input('number of photons ?')
MejDyn  = input('mass of the dynamical ejecta?')
MejWind = input('mass of the wind ejecta?')
Phi     = input('opening angle?')

print(f'filename is: nph{Nph}_mejdyn{MejDyn}_mejwind{MejWind}_phi{Phi}_spec.txt')
data_fold = Path('outputs_file/')
file_ref  = data_fold / f'nph{Nph}_mejdyn{MejDyn}_mejwind{MejWind}_phi{Phi}_spec_sin_ref.txt'
#file_sin  = data_fold / f'nph{Nph}_mejdyn{MejDyn}_mejwind{MejWind}_phi{Phi}_spec_sin.txt'
file_sin   = data_fold /f'nph1.0e+04_mejdyn0.010_mejwind0.030_phi45_spec.txt'
data_time = []
data_magn = []
for f in ['g', 'r', 'i', 'h']:
    file_data_f = data_fold / f'data_{f}.txt'
    time_data_f = np.genfromtxt(fname = file_data_f, skip_header = 1, dtype = 'float', usecols = (0))
    magn_data_f = np.genfromtxt(fname = file_data_f, skip_header = 1, dtype ='float', usecols = (1))
    data_time.append(time_data_f)
    data_magn.append(magn_data_f)

#file_uni  = data_funier / 'nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_spec_ref.txt'
#file_sin  = data_funier /'nph1.0e+04_mej0.04_phi30_T5.0e+03_aopac1.0_spec_sin_2.txt'
#-------------------------------------------------------
nb_obs  = np.genfromtxt(fname=file_ref, max_rows=1, dtype='int')
#nb_obs = 1
nb_bins = np.genfromtxt(fname=file_ref, max_rows=1, skip_header=1, dtype='int')
Nstep   = np.genfromtxt(fname=file_ref, max_rows=1, skip_header=2, dtype='int', usecols=(0))
ti      = np.genfromtxt(fname=file_ref,max_rows= 1, skip_header=2, dtype='float', usecols=(1))
tf      = np.genfromtxt(fname=file_ref,max_rows= 1, skip_header=2, dtype='float', usecols=(2))
wave    = np.genfromtxt(fname=file_ref,max_rows=nb_bins, skip_header=3, dtype='float', usecols= (0))
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
intensity_ref =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity_ref[l,i] = list(np.genfromtxt(fname=file_ref, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

intensity_sin =np.zeros((nb_obs, Nstep), dtype=list)
for l in range(nb_obs):
    for (t,i) in zip(ph,index):
        intensity_sin[l,i] = list(np.genfromtxt(fname=file_sin, skip_header=3 + l*nb_bins, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100))    # factor of scale at 40 Mpc

#--------------------------------------------------------------------------------------------------------------
#                   Light Curves
#--------------------------------------------------------------------------------------------------------------

filter = ['sdss::g','sdss::r', 'sdss::i','swope2::h']
#filter = ['sdss::g']
m_ref = []
m_sin = []
for (f, band) in enumerate(filter):

    m_ref_f = np.zeros((nb_obs, Nstep), dtype=list)
    m_sin_f = np.zeros((nb_obs, Nstep), dtype=list)

    for l in range(nb_obs):
        fluxes_ref_l = np.zeros((Nstep, nb_bins),dtype=list)
        fluxes_sin_l = np.zeros((Nstep, nb_bins),dtype=list)
        for i in index:
            fluxes_ref_l[i] = (intensity_ref[l,i])
            source_ref_l = sn.TimeSeriesSource(ph,wave,fluxes_ref_l)
            fluxes_sin_l[i] = (intensity_sin[l,i])
            source_sin_l = sn.TimeSeriesSource(ph,wave,fluxes_sin_l)
        m_ref_f[l,:]       = list(source_ref_l.bandmag(band,"ab",ph))
        m_sin_f[l,:] =  list(source_sin_l.bandmag(band,"ab",ph))
    m_ref.append(m_ref_f)
    m_sin.append(m_sin_f)


#---------------------------------------------------------------------------------------------------------------
#               PLOT
#--------------------------------------------------------------------------------------------------------------

fig = plt.figure()
plt.subplots_adjust(wspace=0.15, hspace=0.2) #wspace = écartement horizontal, #vspace=vertical ?



for l in range(nb_obs):
    if l == 0 and nb_obs!=1:
        obser = 'edge on'
    elif l>1 and nb_obs!=1 :
        obser = 'face on'
    if nb_obs ==1:
        obser = 'face on'
l = 0
obser = 'face on'
for (band,f, i) in zip(filter, [0,1,2,3], range(l*4,l*4+4,1)):
    dataY_ref = list(m_ref[f][1])
    dataY_sin = list(m_sin[f][1])
    numPlot = 221 + i
    Ax, LW = asManyPlots(numPlot, [ph, ph, data_time[f]],[dataY_ref, dataY_sin, data_magn[f]],
                         markerSize = [0,0,8], marker ='o', color=['lightsteelblue', 'navy','slateblue'], linestyle =['-','-', 'None'],
                         fillstyle = ['full', 'full', 'none'],unfilledFlag = ['False', 'False', 'True'],
                         linewidth=1.5, textsize= 8, legendTextSize=6,
                         plotFlag=[True, True, True], showLegend=False, tickSize = 7,
                         xlim = [0,15], ylim = [16,24])
    Ax.text(4.5, 17 ,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=6)
    plt.gca().invert_yaxis()
#Ax.xaxis.set_major_locator(MultipleLocator(2))
    plt.xticks(np.arange(0,15,2))
    Ax.xaxis.set_minor_locator(MultipleLocator(1))
    Ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    if i==0:
        Ax.set_ylabel('Magnitude', fontsize =8)
    if i==2:
        Ax.set_ylabel('Magnitude', fontsize =8)
Ax.legend([r' Opacity model from Bulla +2019 ', 'New opacity', 'data from GW170817'], fontsize = 6, loc=8)
Ax.set_xlabel('Time since merger (in days)', fontsize=10)
Ax.xaxis.set_label_coords(0,-0.1)
Ax.set_title(f'Light curves drawn {obser}', x = 0, y =2.3, fontsize = 10)


fig.show()
output = f'lc_oldPossis_vs_NewPossis_1e6_same_geo.pdf'
plt.savefig(output, bbox_inches='tight')
sys.exit(0)
    
for (band,f, i) in zip(filter, [0,1,2,3], range(l*4,l*4+4,1)):
    dataY_uni = list(m_uni[f][l])
    dataY_sin = list(m_sin[f][l])
    numPlot = 141 + i
    if i == 0 or i == 4:
        Ax, LW = asManyPlots(numPlot, [ph, ph, time_data],[dataY_uni, dataY_sin, magn_data] , ylabel = 'Magnitude', xlabel = 'time since merger (in days)',
                                markerSize = [0,0,8], marker ='o', color=['rebeccapurple', 'royalblue','thistle'], linestyle =['-','-', 'None'],
                                fillstyle = ['full', 'full', 'none'],unfilledFlag = ['False', 'False', 'True'],
                                lisinidth=2, textsize= 8, legendTextSize=6, title = f'Lights curves drawn {obser}', titlesize = 12,
                                plotFlag=[True, True, True], showLegend=False, tickSize = 7,
                                xlim = [0,7], ylim = [16,24])
    else :
        Ax, LW = asManyPlots(numPlot, [ph, ph, time_data],[dataY_uni, dataY_sin, magn_data] ,
                            markerSize = [0,0,8], marker ='o', color=['rebeccapurple', 'royalblue','plum'], linestyle =['-','-', 'None'],
                                fillstyle = ['full', 'full', 'none'],unfilledFlag = ['False', 'False', 'True'],
                                lisinidth=2, textsize= 8, legendTextSize=6,
                                plotFlag=[True, True, True], showLegend=False, tickSize = 7,
                                xlim = [0,7], ylim = [16,24])

        #Ax, LW = asManyPlots2(numplot, [ph, ph, time_data],[dataY_uni, dataY_sin, magn_data], dataProperties = {'color' :  ['sandybrown','crimson','blue'], 'markerSize': [8, 8, 8],  'unfillMarker': [False, False, True], 'type' : ['plot', 'plot', 'plot']},
        #               generalProperties = {'linestyle': ['-', '-', 'None'], 'textsize': 8}, axesProperties = {'yLabelTextSize' : 6})
        #else:
        #    Ax, LW = asManyPlots(numPlot, [ph, ph], [dataY_uni, dataY_sin],
        #                         markerSize = [0,0], color=['crimson', 'sandybrown'], linestyle =['-','-'], lisinidth=1, textsize= 8, legendTextSize=7,
        #                         plotFlag=[True, True], showLegend=False, tickSize = 7)

        if band == 'sdss::u':
            Ax.text(1.5, 30.5 ,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=6)
        if band == 'sdss::g':
            Ax.text(4.5, 18 ,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=8)
        if i == 2:
            Ax.text(4.5, 18.5,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=6)
        if i == 3:
            Ax.text(3, 19.3,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=6)
        if i == 6:
            Ax.text(5, 18,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=6)
        if i == 7:
            Ax.text(4, 17.8,f'filter: {band}', bbox=dict(facecolor='thistle', alpha=0.2, boxstyle='round'), size=6)

        plt.gca().invert_yaxis()
        Ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        Ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        if i == 0:
            Ax.legend([r'with $\epsilon_{\rm{thermal}} = 0.5$', r' with $\epsilon_{\rm{thermal}}$ from Barnes et al. 2016', 'data from GW170817'], fontsize = 8)
            Ax.set_title(f'Lights curves drawn {obser}', x = 2.5, y =1.15, fontsize = 10)
            Ax.set_xlabel('time since merger (in days)', fontsize = 8)
            Ax.xaxis.set_label_coords(2.5,-0.15)


#fig.suptitle(f'Lights curves drawn {obser}', x =0.5, y =0.5, fontsize = 10)
fig.show()
#output = f'diff_band_effects_of_the_sin_eps.pdf'
output = f'sin_eps_{M_ej}.pdf'
#plt.savefig(output, bbox_inches='tight')
