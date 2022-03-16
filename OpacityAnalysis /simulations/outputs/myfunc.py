###### MATTIA BULLA FUNCTIONS ########

import numpy as np
import scipy as sp
from pylab import *
import sncosmo

def HMS2deg(ra='', dec=''):
  RA, DEC, rs, ds = '', '', 1, 1
  if dec:
    D, M, S = [float(i) for i in dec.split(":")]
    if str(D)[0] == '-':
      ds, D = -1, abs(D)
    deg = D + (M/60) + (S/3600)
    DEC = '{0}'.format(deg*ds)
  
  if ra:
    H, M, S = [float(i) for i in ra.split(":")]
    if str(H)[0] == '-':
      rs, H = -1, abs(H)
    deg = (H*15) + (M/4) + (S/240)
    RA = '{0}'.format(deg*rs)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

def ebv_weight_avg(phase,ebv,sebv,phase_max):

    # Average epochs earlier than phase_max
    if (phase[0]<phase_max):
        ebv = ebv[np.where(phase<phase_max)]
        sebv = sebv[np.where(phase<phase_max)]
    # Average all epochs
    else:
        ebv = ebv
        sebv = sebv
    
    w = 1. / sebv**2

    ebv_avg = np.average(ebv,weights=w)

    return ebv_avg

def ebv_weight_avg_last(phase,ebv,sebv,phase_max):

    # Average epochs earlier than phase_max
    if (phase[-1]>phase_max):
        ebv = ebv[np.where(phase>phase_max)]
        sebv = sebv[np.where(phase>phase_max)]
    # Average all epochs
    else:
        ebv = ebv
        sebv = sebv
    
    w = 1. / sebv**2

    ebv_avg = np.average(ebv,weights=w)

    return ebv_avg

##########################################################
##################### POLARIZATION #######################
##########################################################

def chi_func_0_plus180(qt,ut):
    
    chit=[]

    for j in range(0,len(qt)):

        # This finds the polarisation angle chi, angle between the electric field and the north (l)
        # Chi goes from the north in the counterclockwise direction (see eq 6 of Bulla+2015) 
        if (qt[j]>0 and ut[j]<0):
            pol_ang = 1./2. * ( 2*np.pi - np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) ) 
            chit.append(pol_ang*180/math.pi)
        elif (qt[j]<0 and ut[j]<0):
            pol_ang = 1./2. * ( np.pi + np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) ) 
            chit.append(pol_ang*180/math.pi)
        elif (qt[j]<0 and ut[j]>0):
            pol_ang = 1./2. * ( np.pi - np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) ) 
            chit.append(pol_ang*180/math.pi)
        elif (qt[j]>0 and ut[j]>0):
            pol_ang = 1./2. * np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) 
            chit.append(pol_ang*180/math.pi)
        elif (qt[j]==0):
            pol_ang = 0.25 * math.pi
            if (ut[j]<0):
                pol_ang = 0.75 * math.pi
            chit.append(pol_ang*180/math.pi)
        elif (ut[j]==0):
            pol_ang = 0
            if (qt[j]<0):
                pol_ang = 0.5 * math.pi
            chit.append(pol_ang*180/math.pi)
        # This is for nan values
        else:
            chit.append('nan')

    return chit


def chi_func_min90_plus90(qt,ut):
    
    chit=[]

    for j in range(0,len(qt)):

        # This finds the polarisation angle chi, angle between the electric field and the north (l)
        # Chi goes from the north in the counterclockwise direction (see eq 6 of Bulla+2015) 
        if (qt[j]>0 and ut[j]<0):
            pol_ang = 1./2. * ( 2*np.pi - np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) ) 
            if (pol_ang > math.pi/2.):
                chit.append(pol_ang*180/math.pi-180)    
            else:
                chit.append(pol_ang*180/math.pi)
        elif (qt[j]<0 and ut[j]<0):
            pol_ang = 1./2. * ( np.pi + np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) ) 
            if (pol_ang > math.pi/2.):
                chit.append(pol_ang*180/math.pi-180)    
            else:
                chit.append(pol_ang*180/math.pi)
        elif (qt[j]<0 and ut[j]>0):
            pol_ang = 1./2. * ( np.pi - np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) ) 
            if (pol_ang > math.pi/2.):
                chit.append(pol_ang*180/math.pi-180)    
            else:
                chit.append(pol_ang*180/math.pi)
        elif (qt[j]>0 and ut[j]>0):
            pol_ang = 1./2. * np.arctan(np.absolute(ut[j])/np.absolute(qt[j])) 
            if (pol_ang > math.pi/2.):
                chit.append(pol_ang*180/math.pi-180)    
            else:
                chit.append(pol_ang*180/math.pi)
        elif (qt[j]==0):
            pol_ang = 0.25 * math.pi
            if (ut[j]<0):
                pol_ang = 0.75 * math.pi
            if (pol_ang > math.pi/2.):
                chit.append(pol_ang*180/math.pi-180)    
            else:
                chit.append(pol_ang*180/math.pi)
        elif (ut[j]==0):
            pol_ang = 0
            if (qt[j]<0):
                pol_ang = 0.5 * math.pi
            if (pol_ang > math.pi/2.):
                chit.append(pol_ang*180/math.pi-180)    
            else:
                chit.append(pol_ang*180/math.pi)
        # This is for nan values
        else:
            chit.append('nan')

    return chit


##########################################################
##################### LC TEMPLATES #######################
##########################################################

def read_hsiao(tfit_min,tfit_max,deg,stretch):

    ind_tBmax_templ = 20
    tBmax_templ = 20

    nwave = 5
    indB = 1

    tobs = [0] * nwave
    lc = [0] * nwave

    a=genfromtxt("/Users/mattia/Documents/Work/Work_modelling/toymodels/dust3d/work/obs/templates/hsiao-20070811_lc.dat")

    for i in range(0,nwave):
        
        tobs[i] = a[:,0] 
        lc[i] = -2.5 * log10(a[:,i+1])

        #### CALCULATE DELTAM15 TEMPLATE
        if (i==indB):
        
            w=list(zip(tobs[indB]-tobs[indB][ind_tBmax_templ],lc[indB]))
            zz=list(filter(lambda r: r[0]>tfit_min and r[0]<tfit_max,w))
            x,y=list(zip(*zz))

            p = np.polyfit(x,y,deg=deg, rcond=None, full=False, w=None, cov=False)
            pn = np.poly1d(p)
            xx = linspace(tfit_min,tfit_max,1e4)

            maxb=min(pn(xx))

            w=list(zip(xx,pn(xx)))
            z1=list(filter(lambda r: r[1]==maxb,w))

            deltam15_original = pn(z1[0][0]+15)-z1[0][1]                

            #ax_lc.plot(tobs[indB],lc[indB],color="black",linewidth=2,alpha=0.5)
            #ax_lc.plot(xx,pn(xx),color="black",linewidth=2)


        ### STRETCH LIGHT CURVE

        tobs[i] = tobs[i] * stretch
        tobs[i] = tobs[i] - tobs[i][0]

    
    tBmax_templ = tBmax_templ * stretch

    return tobs,lc,deltam15_original,ind_tBmax_templ,tBmax_templ


def dm15(t,mag,tfit_min,tfit_max,deg):

    w=list(zip(t,mag))
    zz=list(filter(lambda r: r[0]>tfit_min and r[0]<tfit_max,w))
    x,y=list(zip(*zz))

    p = np.polyfit(x,y,deg=deg, rcond=None, full=False, w=None, cov=False)
    pn = np.poly1d(p)
    xx = linspace(tfit_min,tfit_max,1e4)

    maxb=min(pn(xx))

    w=list(zip(xx,pn(xx)))
    z1=list(filter(lambda r: r[1]==maxb,w))

    tBmax = z1[0][0]
    deltam15 = pn(z1[0][0]+15)-z1[0][1] 

    return deltam15, tBmax , pn           



##########################################################
####################### PLOTTING #########################
##########################################################


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def rebin(vector,bins_ratio):

    old_size = int(len(vector))
    new_size =  int(old_size / bins_ratio)

    new_vector = []

    if (old_size%new_size!=0):
        
        ### Fill vector with 0 at the end
        for i in range(0,new_size-old_size%new_size):       
            
            if (i==0):
                new_vector=np.append(vector,1e50)
            else:
                new_vector=np.append(new_vector,1e50)

        x=new_vector.reshape(new_size,int(old_size/new_size)+1)
        final_vector=x.mean(1)

        final_vector = final_vector[np.where(final_vector!=1e50)][0:-1]
    
    else:
        x=vector.reshape(new_size,int(old_size/new_size))
        final_vector=x.mean(1)

    return final_vector


def fill_between_steps(x, y1, y2=0, h_align='mid', ax=None, **kwargs):
    ''' Fills a hole in matplotlib: fill_between for step plots.

    Parameters :
    ------------

    x : array-like
        Array/vector of index values. These are assumed to be equally-spaced.
        If not, the result will probably look weird...
    y1 : array-like
        Array/vector of values to be filled under.
    y2 : array-Like
        Array/vector or bottom values for filled area. Default is 0.

    **kwargs will be passed to the matplotlib fill_between() function.

    '''
    # If no Axes opject given, grab the current one:

    x = np.array(x)

    if ax is None:
        ax = plt.gca()
    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = sp.repeat((x[1:] - x[:-1]), 2)
    xstep = sp.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    # Now: add one step at end of row.
    xx = sp.append(xx, xx.max() + xstep[-1])

    # Make it possible to chenge step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)
    if type(y2) == sp.ndarray:
        y2 = y2.repeat(2)

    # now to the plotting part:
    ax.fill_between(xx[:-1], y1[:-1], y2=y2[:-1], **kwargs)

    return ax



############################
####### SELECTION ##########
############################

def select_amin_amax(array_all,array_sel,amin,amax):

    array_new = [0] * len(array_all)
    
    for i in range(0,len(array_all)):
        array_new[i] = array_all[i][np.logical_and(array_sel>amin,array_sel<amax)]

    return array_new

def select_amin(array_all,array_sel,amin):

    array_new = [0] * len(array_all)
    
    for i in range(0,len(array_all)):
        array_new[i] = array_all[i][np.where(array_sel>amin)]

    return array_new

def select_equal(array_all,array_sel,value):

    array_new = [0] * len(array_all)
    
    for i in range(0,len(array_all)):
        array_new[i] = array_all[i][np.where(array_sel==value)]

    return array_new


################################################################
#**************************************************************#
#################### ZTF SUPERNOVA ANALYSIS ####################
#**************************************************************#
################################################################



################################################################
####### IDENTIFY SNe ON THE MARSHALL THAT ARE SAME SN ##########
################################################################

def duplicate(SN,ra,dec,dtheta):

    from itertools import combinations
    
    SN_comb = []
    ra_comb = []
    dec_comb = []

    for i in list(combinations(SN,2)):
        SN_comb.append(i)
    for i in list(combinations(ra,2)):
        ra_comb.append(i)
    for i in list(combinations(dec,2)):
        dec_comb.append(i)

    t = np.pi/180.
    
    theta = []

    for i in range(0,len(SN_comb)):
        costheta = np.sin(dec_comb[i][0]*t)*np.sin(dec_comb[i][1]*t)+np.cos(dec_comb[i][0]*t)*np.cos(dec_comb[i][1]*t)*np.cos(ra_comb[i][0]*t-ra_comb[i][1]*t)

        # Theta in arcsec
        theta.append( np.arccos(costheta) / t * 3600 )
    
    SN_dupl = np.array(SN_comb)[np.where(np.array(theta)<dtheta)]
    theta_dupl = np.array(theta)[np.where(np.array(theta)<dtheta)]

    return SN_dupl, theta_dupl


################################################################
##### Table of data in sncosmo format (merging duplicates) #####
################################################################

def table_sncosmo_dupl(data,SN_dupl,nobj,idupl):

    import marshaltools as mt
    from astropy.table import vstack

    ### No duplicates
    if (nobj==1):  

        # Light curves
        classification = data.classification
        z = data.redshift
        ebv_mw = data.mwebv

        lc_data = data.table_sncosmo
       
    ### Duplicates
    elif nobj==2:  

        lcs = mt.ProgramList("Cosmology")
        
        lc1 = lcs.get_lightcurve(SN_dupl[idupl][0])
        lc1_data = lc1.table_sncosmo

        lc2 = lcs.get_lightcurve(SN_dupl[idupl][1])
        lc2_data = lc2.table_sncosmo

        classification = lc1.classification
        z1 = lc1.redshift
        z2 = lc2.redshift
        if z1 == None:
            z = z2
        elif z2 == None:
            z = z1
        else:
            z = 0.5 * (float(z1) + float(z2))

        lc_data = vstack([lc1_data,lc2_data])
        lc_data.meta['z'] = z
        lc_data.meta['mwebv'] = lc1.mwebv   
        
        lc_data.sort('mjd')

        idupl+=1
    else:
        lc_data = None
        classification = None
    
    return lc_data, classification, idupl


################################################################
############### Print messages and saved parameters ############
################################################################

def print_and_save(sn,isn,nSN,data=None,fit=None,par1=None,par2=None,save=None,message=None,filename=None):

    if save == 1:
        out = open(filename+".txt","a+")

    # SN included in the analysis, good fit
    if message == 0:

        print("%s (%d out of %d)"%(sn,isn+1,nSN))

        if save == 1:

            out.write("%s "%sn)
            
            for i in range(0,len(fit[1].parameters)):
                out.write("%f "%fit[1].parameters[i])
            
            out.write("%f %f "%(par1,par2))
            
            lc_datag = data[np.where(data['band']=="p48g")]
            lc_datar = data[np.where(data['band']=="p48r")]
            
            if len(lc_datag)!=0 and len(lc_datar)!=0:
                out.write("%f "%(lc_datag['mjd'][0] - fit[1].parameters[1]))
                out.write("%f "%(lc_datar['mjd'][0] - fit[1].parameters[1]))
            elif len(lc_datag)==0:
                out.write("1e6 ")
                out.write("%f "%(lc_datar['mjd'][0] - fit[1].parameters[1]))
            elif len(lc_datar)==0:
                out.write("%f "%(lc_datag['mjd'][0] - fit[1].parameters[1]))
                out.write("1e6 ")
            else:
                out.write("1e6 ")
                out.write("1e6 ")  

            out.write("%d "%len(lc_datag) )
            out.write("%d "%len(lc_datar) )

            out.write("%f %d \n"%(fit[0].chisq,fit[0].ndof))

            fig = sncosmo.plot_lc(data, model=fit[1],errors=fit[0].errors,xfigsize=15,tighten_ylim=True)
            fig.savefig("fits/"+sn+".pdf",bbox_inches='tight')
                              
            
    elif message == -1:

        print("%s not considered, duplicate "%sn)
        
        if save == 1:
            out.write("#%s not considered, duplicate \n"%sn) 


    # Not a normal Ia
    elif message == 1:

        print("%s not considered, not a SN Ia norm (%s)"%(sn,par1))
        
        if save == 1:
            out.write("#%s not considered, not a SN Ia norm (%s) \n"%(sn,par1))     

    # Missing redshift
    elif message == 2:

        if par2 == 1:
            print("%s not considered, missing redshift (spectra available) "%sn)
            if save == 1:
                out.write("#%s not considered, missing redshift (spectra available) \n"%sn) 
        else:
            print("%s not considered, missing redshift "%sn)
            if save == 1:
                out.write("#%s not considered, missing redshift \n"%sn)     

    # Bad fit
    elif message == 3:

        print("%s not considered, bad fit"%sn)

        if save == 1:
            out.write("#%s not considered, bad fit \n"%sn)

    # No data points
    elif message == 4:

        print("%s not considered, no data points \n"%sn)

        if save == 1:
            out.write("#%s not considered, no data points \n"%sn)


    if save == 1:
        out.close()


################################################################
############# Line argument parameters for hdiag.py ############
################################################################

def hdiag_params():

    x1cut = raw_input("Would you like a cut on x1 [y/n] ? ")

    if x1cut == "y":
        x1_min = np.float(input("Enter min(x1): "))
        x1_max = np.float(input("Enter max(x1): "))
    else:
        x1_min = -11
        x1_max = 11

    chicut = raw_input("Would you like a cut on the reduced chisquare from SALT2 fitting [y/n] ? ")

    if chicut == "y":
        chisqred_max = np.float(input("Enter max(chi_sq_red): "))
    else:
        chisqred_max = 1e5

    phasecut = raw_input("Would you like a cut on the maximum phase for the earliest g/r data point [y/n] ? ")

    if phasecut == "y":
        phasemax_g = np.float(input("Enter maximum phase in g band: "))
        phasemax_r = np.float(input("Enter maximum phase in r band: "))
    else:
        phasemax_g = 1e5
        phasemax_r = 1e5

    ncut = raw_input("Would you like a cut on the minimum number of g/r data points [y/n] ? ")

    if ncut == "y":
        nmin_g = np.int(input("Enter minimum number of data points in g band: "))
        nmin_r = np.int(input("Enter minimum number of data points in r band: "))
    else:
        nmin_g = 0
        nmin_r = 0

    return  x1_min, x1_max, chisqred_max, phasemax_g, phasemax_r, nmin_g, nmin_r


def _get_bandmag(band, magsys, t=0, rest_frame=True, **kwargs):
    """
    Returns mag at max for the model, band and magsys
    Arguments:
    model  -- sncosmo model, e.g. SALT2
    band   -- sncosmo band object or string, e.g. 'bessellb'
    magsys -- magnitude system, e.g. 'ab'
    Keyword arguments:
    t -- time relative to t0 (observer-frame), at which to evaluate
    rest_frame -- default: True, overrides the redshifts
    """
    
    model = sncosmo.Model(source='salt2')
    if rest_frame:
        kwargs['z'] = 0

    model.set(**kwargs)
    return model.bandmag(band,magsys,kwargs['t0'] + t)

def _get_bandmag_gradient(band, magsys, param, sig, fixed_param, 
                          t=0, rest_frame=True):
    """
    Return gradient of _get_bandmag as function of param
    param, sig must be dictionaries of means and uncertainties
    Best use odicts to make sure that the order of the components is correct
    """

    model = sncosmo.Model(source='salt2')
    out = []
    
    if rest_frame:
        if 'z' in param.keys():
            param['z'] = 0
        if 'z' in fixed_param.keys():
            fixed_param['z'] = 0

    model.set(**fixed_param)
    for key,val in param.items():
        model.set(**param)
        h = sig[key] / 100.
        
        model.set(**{key: val - h})
        m0 = model.bandmag(band, magsys, param['t0'] + t)

        model.set(**{key: val + h})
        m1 = model.bandmag(band, magsys, param['t0'] + t)
        
        out.append((m1 - m0) / (2. * h))

    return np.array(out)

################################################################
######################## SNID TYPING ###########################
################################################################

def snid_select(file,rlap,sn_templ,clas_templ,z,dz,phase,lap_min,rlap_min):
    """
    Select all the SNID templates with lap>lap_min rlap > rlap_min
    """
    
    stop = 0
    
    with open(file) as fp:  
        for cnt, line in enumerate(fp):

            if line.startswith('#---'):
                stop = 1

            if (cnt>69 and stop==0): 

                laptmp = float(line.split()[3])
                rlaptmp = float(line.split()[4])
                quality = str(line.split()[9])
                ztmp = float(line.split()[5])

                if (ztmp > 0 and ztmp < 0.15 and laptmp>=lap_min and rlaptmp>=rlap_min and quality=="good"):

                    rlap.append(rlaptmp)
                    sn_templ.append( str(line.split()[1]) ) 
                    clas_templ.append( str(line.split()[2]) ) 
                    z.append( ztmp )
                    dz.append( float(line.split()[6]) )
                    phase.append( float(line.split()[7]) )



def snid_type(rlap,sn_templ,clas_templ,z,dz,ntempl):
    """
    Select best SNID type
    """
   
    import pandas as pd
    import collections

    ### Order by decreasing rlap
    rlap,sn_templ,clas_templ,z,dz = zip(*sorted(zip(rlap,sn_templ,clas_templ,z,dz),reverse=True))
    
    ### Remove duplicate templates
    d = {'col1': sn_templ, 'col2': clas_templ, 'col3': rlap, 'col4': z, 'col5': dz}
    df = pd.DataFrame(data=d)
    dnew = df.drop_duplicates(cols='col1')
    sn_templ = dnew['col1']
    clas_templ = dnew['col2']
    rlap = dnew['col3']
    z = dnew['col4']
    dz = dnew['col5']

    ### Classify
    types = collections.Counter(clas_templ)
    rlap_mean = []
    clas_type = []
    z_type = []
    dz_type = []

    # How many template to use
    nmax_first = types.most_common(len(types))[0][1]

    if (nmax_first > ntempl):
        nchosen = ntempl
    else:
        nchosen = nmax_first

    for itype in range(0,len(types)):
        # type
        clastmp = types.most_common(len(types))[itype][0]
        # rlap values for given type
        rlap_type = np.array(rlap)[np.where(np.array(clas_templ) == clastmp)]
        ztmp = np.array(z)[np.where(np.array(clas_templ) == clastmp)]
        dztmp = np.array(dz)[np.where(np.array(clas_templ) == clastmp)]
        # Save clas
        clas_type.append( clastmp )
        # Take only the first nchosen template for each 
        #if len(rlap_type)>=nchosen:
        z_type.append( ztmp[:nchosen] )
        dz_type.append( dztmp[:nchosen] )
        rlap_mean.append( mean(rlap_type[:nchosen]) ) # mean rlap value
    
    # Save best type and relevant parameters
    clas = np.array(clas_type)[np.where(np.array(rlap_mean)==max(np.array(rlap_mean)))][0]
    rlap_clas = np.array(rlap_mean)[np.where(np.array(rlap_mean)==max(np.array(rlap_mean)))][0]
    z_clas = np.array(z_type)[np.where(np.array(rlap_mean)==max(np.array(rlap_mean)))][0]
    dz_clas = np.array(dz_type)[np.where(np.array(rlap_mean)==max(np.array(rlap_mean)))][0]

    # Classification probability (take the first nchosen entries)
    first_ntempl = np.array(rlap[:nchosen])
    prob = sum( first_ntempl[np.where(clas_templ[:nchosen]==clas)] ) / sum(first_ntempl)

    # Weighted mean redshift
    w = 1/np.array(dz_clas)**2
    zmean = average(np.array(z_clas),weights=w)
    dzmean = 1. / sqrt( sum(w) )

    return clas,rlap_clas,prob,zmean,dzmean,nchosen


def snid_prob_type(rlap,sn_templ,clas_templ,ntempl,which_type,z,phase):
    """
    Select a given (not necessarily best, see snid_type) type and output its probability
    """
   
    import pandas as pd

    ### Order by decreasing rlap
    rlap2,sn_templ2,clas_templ2,z2,phase2 = zip(*sorted(zip(rlap,sn_templ,clas_templ,z,phase),reverse=True))
    
    ### Remove duplicate templates
    d = {'col1': sn_templ2, 'col2': clas_templ2, 'col3': rlap2, 'col4': z2 ,'col5': phase2}
    df = pd.DataFrame(data=d)
    dnew = df.drop_duplicates(cols='col1')
    clas_templ2 = dnew['col2']
    rlap2 = dnew['col3']
    z2 = dnew['col4']
    phase2 = dnew['col5']

    # Classification probability (take the first ntempl entries)
    clas = np.array(clas_templ2[:ntempl])
    rlap_all = np.array(rlap2[:ntempl])
    rlap_type = []
    z_type = []
    phase_type = []
    for i in range(0,len(clas)):
        if (which_type in clas[i]):
            rlap_type.append(rlap_all[i])
            z_type.append(np.array(z2[:ntempl])[i])
            phase_type.append(np.array(phase2[:ntempl])[:ntempl][i])

    rlap_clas = mean(rlap_type)
    prob = sum( rlap_type ) / sum(rlap_all)

    z_type = mean(z_type)
    phase_type = mean(phase_type)

    return rlap_clas,z_type,phase_type,prob,ntempl




