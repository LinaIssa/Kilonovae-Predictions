# -*- coding: utf-8 -*-
"""
Created on Dec the 27 2019
@author: Lina Issa
"""

import numpy as np
import sncosmo as sn
import matplotlib.pyplot as plt
import glob
import itertools

#********************************************



filenames = sorted(glob.glob('flux_spectra_*.txt'))
# ask nb_bins
nb_bins = 50
fluxes = []
wavelength = np.loadtxt(fname = 'flux_spectra_1day.txt', skiprows=1)[:,2]
ph = np.array([1,2,4,7])
filter = 'sdss::g'
for f in filenames:
    data = np.loadtxt(fname=f, skiprows=1)
    intensity = list(data[:,1])
    fluxes.append(intensity)    # liste de spectres pour chaque temps 
source = sn.TimeSeriesSource(ph,wavelength,np.array(fluxes))
m = source.bandmag(filter,"ab",ph)
print(list(m))

