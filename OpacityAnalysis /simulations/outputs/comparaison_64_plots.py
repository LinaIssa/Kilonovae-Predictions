import numpy as np
import sncosmo as sn
import matplotlib.pyplot as plt
import matplotlib
import glob
from plotUtilities import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import matplotlib.colors as colors
import matplotlib.cm as cmx
import sys
from palettable.cartocolors.diverging import ArmyRose_7, Tropic_6
from palettable.scientific.sequential import Batlow_3, LaPaz_20
import os
from pathlib import Path


#-------------------------------------------------------------------------------------------------------------
#                   Extract Data
#--------------------------------------------------------------------------------------------------------------

Nph = f'{1.0e+05:.1e}'
MejDyn  = f'{0.010:.3f}'
MejWind = f'{0.030:.3f}'
Phi = 45

kappaLF = 100
#KappaLF = [1, 10, 100, 1000]
GammaLF = [-1.0,-0.7,-0.4,-0.1]
KappaLR = [1, 10, 100, 1000]
GammaLR = [-1.0,-0.7,-0.4,-0.1]

Nb_para  = 4
Nb_obser = 1
Nb_sim   = 64


data_fold1 = Path('outputs_file/')
file_ref  = data_fold1 / f'nph{Nph}_mejdyn{MejDyn}_mejwind{MejWind}_phi{Phi}_spec_ref.txt'

data_fold2 = Path('outputs_nph5/')

data_time = []
data_magn = []
for f in ['g', 'r', 'i', 'h']:
    file_data_f = data_fold1 / f'data_{f}.txt'
    time_data_f = np.genfromtxt(fname = file_data_f, skip_header = 1, dtype = 'float', usecols = (0))
    magn_data_f = np.genfromtxt(fname = file_data_f, skip_header = 1, dtype ='float', usecols = (1))
    data_time.append(time_data_f)
    data_magn.append(magn_data_f)


#-------------------------------------------------------
nb_bins = np.genfromtxt(fname=file_ref, max_rows=1, skip_header=1, dtype='int')
Nstep   = np.genfromtxt(fname=file_ref, max_rows=1, skip_header=2, dtype='int', usecols=(0))
ti      = np.genfromtxt(fname=file_ref,max_rows= 1, skip_header=2, dtype='float', usecols=(1))
tf      = np.genfromtxt(fname=file_ref,max_rows= 1, skip_header=2, dtype='float', usecols=(2))
wave    = np.genfromtxt(fname=file_ref,max_rows=nb_bins, skip_header=3, dtype='float', usecols= (0))


#-------------------------------------------------------
step   = (tf-ti)/Nstep
ph     = np.arange(ti+step/2, tf + step/2, step)
index  = np.arange(Nstep)


#-------------------------------------------------------------------------------------------------------------
#                   Spectra
#--------------------------------------------------------------------------------------------------------------
intensity_ref = [] #np.zeros((1, Nstep), dtype=list)
for i in index:
    intensity_ref.append(list(np.genfromtxt(fname=file_ref, skip_header=3, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100)))    # factor of scale at 40 Mpc


#--------------------------------------------------------------------------------------------------------------
#                   Light Curves
#--------------------------------------------------------------------------------------------------------------

filter = ['sdss::g','sdss::r', 'sdss::i','swope2::h']
#filter = ['sdss::g']
m_ref = []
for (f, band) in enumerate(filter):
    fluxes_ref = np.zeros((Nstep, nb_bins),dtype=list)
    for i in index:
        fluxes_ref[i,:] = intensity_ref[i] # list des intensités I_lambda pour un temps donné
    source_ref = sn.TimeSeriesSource(ph,wave,fluxes_ref)
    m_ref_f      = list(source_ref.bandmag(band,"ab",ph))
    m_ref.append(m_ref_f)



SimulationsData = []
    #for kappaLF in KappaLF:
for gammaLF in GammaLF:
    for kappaLR in KappaLR:
        for gammaLR in GammaLR:
            file_sim   = data_fold2 /f'nph1.0e+05_kappaLF{kappaLF}_gammaLF{gammaLF}_kappaLR{kappaLR}_gammaLR{gammaLR}_spec.txt'
            intensity_sim = []
            for i in index:
                intensity_sim.append(list(np.genfromtxt(fname=file_sim, skip_header=3, max_rows=nb_bins, dtype='float', usecols= (1+ i*3))/(16*100)))    # factor of scale at 40 Mpc
            m_sim = []
            for (f, band) in enumerate(filter):
                fluxes_sim = np.zeros((Nstep, nb_bins),dtype=list)
                for i in index:
                    fluxes_sim[i,:] = intensity_sim[i] # list des intensités I_lambda pour un temps donné
                source_sim = sn.TimeSeriesSource(ph,wave,fluxes_sim)
                m_sim_f      = list(source_sim.bandmag(band,"ab",ph))
                m_sim.append(m_sim_f)
            SimulationsData.append(m_sim)


#---------------------------------------------------------------------------------------------------------------
#               PLOT
#--------------------------------------------------------------------------------------------------------------

fig = plt.figure()
plt.subplots_adjust(wspace=0.15, hspace=0.2) #wspace = écartement horizontal, #vspace=vertical ?

values = range(len(SimulationsData))
#cm = plt.get_cmap('plasma')
#cm = LaPaz_20.mpl_colormap
#cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
#print scalarMap.get_clim()

for (band,f, i) in zip(filter, [0,1,2,3], range(0,4,1)):
    dataY_ref = m_ref[f]
    numPlot = 221 + i
    for s in range(len(SimulationsData)):
        dataY_sim = SimulationsData[s][f]
        #colorVal = scalarMap.to_rgba(values[s])
        Ax, LW = asManyPlots(numPlot, [ph, ph, data_time[f]],[dataY_sim, dataY_ref,data_magn[f]], markerSize = [0,0,5], marker ='o', color=['lavender', 'navy','navy'], linestyle =['-','-', 'None'], fillstyle = ['full', 'full', 'none'],unfilledFlag = ['False', 'False', 'True'], linewidth=0.6, textsize= 8, legendTextSize=6, plotFlag=[True, True, True], showLegend=False, tickSize = 7,xlim = [0,15], ylim = [16,24])
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
Ax.legend([r'simulations runs for $\kappa_{\rm{LF}} = $' + f'{kappaLF} ' + r'$\rm{cm^2/g}$', 'opacity model from Bulla +2019', 'data from GW170817'], fontsize = 6, loc=8)
Ax.set_xlabel('Time since merger (in days)', fontsize=10)
Ax.xaxis.set_label_coords(0,-0.1)
Ax.set_title(r'Light curves drawn for $\theta_{\rm{obs}} = 17 \degree$', x = 0, y =2.3, fontsize = 10)


fig.show()
output = f'lc_nph{Nph}_kappaLF{kappaLF}.pdf'
plt.savefig(output, bbox_inches='tight')

#
