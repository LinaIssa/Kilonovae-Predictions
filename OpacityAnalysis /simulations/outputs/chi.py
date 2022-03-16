import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import stats
import pandas as pd
from pylab import *
import time
import os
import sncosmo
import matplotlib.cm as cm
import pandas as pd
import matplotlib.gridspec as gridspec
import numpy.ma as ma
import extinction

plt.rcdefaults()
plt.rc('lines', linewidth=1.5)
plt.rc('axes', labelsize=19)
plt.rc('xtick', labelsize=19)
plt.rc('ytick', labelsize=19)
plt.rc('text', usetex=True)
plt.rc('font',family='serif')
plt.rc('font', serif='Times')
plt.rc("legend", fontsize=10.5)

home = "/Users/linaissa/documents/ARPE/my_work/opacities_I/simulations/outputs/"

directory = "outputs_nph5/"

tshift = 0

#ep_max = 70

kappaLF = [10]#, 10, 100, 1000, 10000]#, 100000]
kappaLR = [10]#, 10, 100, 1000] #10000, 100000]
gammaLF = [-0.1]#, -0.7, -0.4,-0.1, 0.2]
gammaLR = [-0.4]#, -0.7, -0.4,-0.1,]

H0 = 73
CLIGHT = 2.99792458e5
dMpc = 40 
z = dMpc / CLIGHT * H0

### FILTERS
filters = ["sdss::u","sdss::g","sdss::r","sdss::i","sdss::z", "swope2::y","swope2::J","swope2::h"]#,"cspk"]
#filters = ['sdss::g','sdss::r', 'sdss::i','swope2::h']

### DATA
datadir = home+"outputs_file/"

mjd = []
phase = []

file = open(datadir+"AT2017gfo_phot_compiled_sjs.dat","r")
for i, line in enumerate(file):
    if not line.startswith("#"):
    	mjd.append( float(line.split()[0]) )
    	phase.append( float(line.split()[1]) )

data_all = []
dataerr_all = []
ph_data_all = []

for ifilt in range(0,len(filters)):

	##### DATA
	data = []
	dataerr = []
	file = open(datadir+"AT2017gfo_phot_compiled_sjs.dat","r")

	for kk, line in enumerate(file):
	    if not line.startswith("#"):
	    	if ifilt<6:
	    		ind = 2+2*ifilt
	    	else:
	    		ind = 4+2*ifilt

	    	if (str(line.split()[ind])=="NaN" or float(line.split()[ind+1])==9999):
	    		data.append( 99.0 )
	    		dataerr.append( 0.0 )	    
	    	else:
	    		data.append( float(line.split()[ind]) )
	    		dataerr.append( float(line.split()[1+ind]) )
	
	ph_data = np.array(phase)[np.where(np.array(data)!=99.0)[0]]
	dataerr = np.array(dataerr)[np.where(np.array(data)!=99.0)[0]]
	data = np.array(data)[np.where(np.array(data)!=99.0)[0]]

	data_all.append(data)
	dataerr_all.append(dataerr)
	ph_data_all.append(ph_data)


##### FIT

## Get relevant parameters

filename = "nph1.0e+05_kappaLF"+"{:.0f}".format(kappaLF[0])+"_gammaLF"+"{:.1f}".format(gammaLF[0])+"_kappaLR"+"{:.0f}".format(kappaLR[0]) +"_gammaLR"+"{:.1f}".format(gammaLR[0])
a = open(home+directory+filename+"_spec.txt")

l=0

for ii in a:
	if l==0:
		Nobs=int(ii)
	elif l==1:
		Nwave=int(ii)
	elif l==2:
		Ntime=int(ii.split()[0])
		ti = double(ii.split()[1])
		tf = double(ii.split()[2])
	l+=1
step_time = (tf-ti) / double(Ntime)

chi2 = [ [ [ [0 for x in range(len(kappaLF))] for y in range(len(kappaLR))] for w in range(len(gammaLF))] for z in range(len(gammaLR))]
#chi2_red = [ [ [[0 for x in range(Nobs)] for y in range(len(phi))] for w in range(len(mej2))] for k in range(len(mej1))]

## Run over different models

for igammaLR in range(0,len(gammaLR)):
    print(gammaLR[igammaLR])
    for igammaLF in range(0,len(gammaLF)):
        for ikappaLR in range(0,len(kappaLR)):
            for ikappaLF in range(0,len(kappaLF)):

			#if (phi[iphi]==0 or phi[iphi]==90):
			#	Nobs=1
                filename = "nph1.0e+05_kappaLF"+"{:.0F}".format(kappaLF[ikappaLF])+"_gammaLF"+"{:.1f}".format(gammaLF[igammaLF])+"_kappaLR"+"{:.0f}".format(kappaLR[ikappaLR]) +"_gammaLR"+"{:.1f}".format(gammaLR[igammaLR])


                a = genfromtxt(home+directory+filename+"_spec.txt",skip_header=3)

                ndata = 0

                for ifilt in range(0,len(filters)):

                    w = a[0:Nwave,0]
                    # Redshift SEDs
                    w *= (1+z)

                    ph = []
                    fl = np.ones((Ntime,len(w)))
					
                    for itime in range(0,Ntime):

                        time = ti + step_time * (itime + 0.5) #+ ph_data[0] + tshift
                        I = a[0:Nwave,1+itime*3] * (1./dMpc)**2

                        ph.append(time)
                        fl[itime] = I

                    source = sncosmo.TimeSeriesSource(ph,w,fl)
                    m = source.bandmag(filters[ifilt],"ab",ph)
                    fmodel  = 10**(-m/2.5)

                    # Compare model to data
                    phdata = ph_data_all[ifilt]
                    mdata = data_all[ifilt]
                    mdataerr = dataerr_all[ifilt]
                    fdata = 10**(-mdata/2.5)
                    fdataerr = log(10) * fdata * mdataerr

                    for idata in range(0,len(phdata)):
						#if (phdata[idata] < ep_max):

                        diff = fabs(np.array(ph) - phdata[idata])
                        ind = np.where(diff==min(diff))[0][0]

                        if fmodel[ind]!=fmodel[ind]:
                            fmodel[ind] = 1e100

                        chi2[igammaLR][igammaLF][ikappaLR][ikappaLF] += ( fmodel[ind] - fdata[idata])**2 / fdataerr[idata]**2
                        ndata += 1
#print(ndata, chi2,m[ind], mdata[idata], mdataerr[idata],(m[ind] - mdata[idata])**2/mdataerr[idata]**2 )
#print(fmodel[ind], fdata[idata], fdataerr[idata])
		
				#if Nobs>1:
				#	chi2_red[imej1][imej2][iphi][iobs] = chi2[imej1][imej2][iphi][iobs] / (ndata - 3)
				#else:
				#	for j in range(0,len(chi2_red[imej1][imej2][iphi])):
				#		chi2_red[imej1][imej2][iphi][j] = chi2[imej1][imej2][iphi][iobs] / (ndata - 3)


######

#P3D = exp(- 0.5 * np.array(chi2_red))

P3D = exp(- 0.5 * np.array(chi2))

P3D /= sum(P3D)
P_kappaLF = P3D.sum(axis=0).sum(axis=0).sum(axis=0)
P_kappaLR = P3D.sum(axis=3).sum(axis=1).sum(axis=0)
P_gammaLF = P3D.sum(axis=3).sum(axis=2).sum(axis=0)
P_gammaLR = P3D.sum(axis=3).sum(axis=2).sum(axis=1)

####### PLOT

fig=plt.figure(figsize=(16,10))

# KappaLF

ax2 = fig.add_subplot(2,2,1)
#ax2.set_xlim(0,0.03)
ax2.set_xlabel("$\kappa_{\rm{lf}$")
ax2.set_ylabel("P",labelpad=10)

majorLocator = MultipleLocator(0.005)
ax2.xaxis.set_major_locator(majorLocator)

ax2.scatter(kappaLF,P_kappaLF)

# KappaLR

ax3 = fig.add_subplot(2,2,2)
#ax3.set_xlim(0,0.15)
ax3.set_xlabel("$\kappa_{\rm{lr}$")
ax3.set_ylabel("P",labelpad=10)

majorLocator = MultipleLocator(0.01)
ax3.xaxis.set_major_locator(majorLocator)

ax3.scatter(kappaLR,P_kappaLR)



ax4 = fig.add_subplot(2,2,3)
#ax4.set_xlim(15,75)
ax4.set_xlabel("gamma_{\rm{lf}")
ax4.set_ylabel("P",labelpad=10)

majorLocator = MultipleLocator(15)
ax4.xaxis.set_major_locator(majorLocator)

ax4.scatter(gammaLF,P_gammaLF)

subplots_adjust(hspace=0.25)
subplots_adjust(wspace=0.35)


ax5 = fig.add_subplot(2,2,4)
#ax5.set_xlim(15,75)
ax5.set_xlabel("gamma_{\rm{lr}")
ax5.set_ylabel("P",labelpad=10)

majorLocator = MultipleLocator(15)
ax5.xaxis.set_major_locator(majorLocator)

ax5.scatter(gammaLR,P_gammaLR)

subplots_adjust(hspace=0.25)
subplots_adjust(wspace=0.35)


print(np.max(P3D))
print(np.where(P3D==np.max(P3D)))
print(kappaLF[np.where(P3D==np.max(P3D))[3][0]])
print(kappaLR[np.where(P3D==np.max(P3D))[2][0]])
print(gammaLF[np.where(P3D==np.max(P3D))[1][0]])
print(gammaLR[np.where(P3D==np.max(P3D))[0][0]])


fig.show()

#### SAVE

#filename = "viewangle_mej1+mej2+phi"
#
#out = open(filename+".txt","w")
#for i in range(0,len(P_angle)):
#	out.write("%f %f \n"%(mu[i],P_angle[i]))
#
#out.close()

#fig.savefig(filename+".pdf")
