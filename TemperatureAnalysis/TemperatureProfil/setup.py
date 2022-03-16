import numpy as np
import scipy as sp
import myfunc as mf
import matplotlib.pyplot as plt
from pylab import *

Npar = 4
Nobs = 11

Mej = np.array([0.04,0.02,0.06,0.08,0.10])
opang_lr = np.array([30])
T0 = np.array([5000])
alpha_opac = np.array([1.0])

Nsim = len(opang_lr) * len(Mej) * len(T0) * len(alpha_opac)

a = open("setup.txt","w")
a.write("%d %d %d \n"%(Npar,Nobs,Nsim))

for j in range(0,len(Mej)):
	for i in range(0,len(opang_lr)):
		for k in range(0,len(T0)):
			for l in range(0,len(alpha_opac)):
				a.write("%g %g %g %g \n"%(Mej[j],opang_lr[i],T0[k],alpha_opac[l]))

a.close()