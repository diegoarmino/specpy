
import numpy as np
import math as mth
import cmath as cmth
import scipy.constants as cnst
from scipy.fftpack import fft, ifft,  fftfreq, fftshift,dct,idct
import matplotlib.pyplot as plt

from ramanlib import uvvis
##############################################################################
#     APPLICATIONS: hexatriene
#     Petrenko & Neese 2007
##############################################################################
omega_nm = np.array([354.,444.,934.,1192.,1285.,1295.,1403.,1581.,1635.])
delta_nm = np.array([0.55,0.23,0.23,0.82,0.485,0.02,0.085,0.38,1.32])
gamma    = 160.
omega0=39805
fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
states_list = np.identity(9)
states_list.tolist()
xmax=48000
xmin=38000
ymax=300
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
freq,spec = uvvis(delta_nm,omega_nm,omega0,gamma)
ax.plot(freq, spec)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION INTENSITY')
plt.show()
