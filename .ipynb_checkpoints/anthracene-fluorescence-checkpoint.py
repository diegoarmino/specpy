
import numpy as np
import math as mth
import cmath as cmth
import scipy.constants as cnst
from scipy.fftpack import fft, ifft,  fftfreq, fftshift,dct,idct
import matplotlib.pyplot as plt
import sys
import os.path
sys.path.append(os.path.abspath('.'))
sys.path.insert(0, "/home/diegoa/dev/specpy")
from ramanlib import fluorescence,uvvis,read_table


##############################################################################
#     APPLICATIONS: hexatriene
#     Petrenko & Neese 2007
##############################################################################

delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/anth-gaus/b3lyp-631/disp.dat')

gamma    = 160.
omega0=26695.
states_list = np.identity(len(omega_nm))
states_list.tolist()
ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uvvis = uvvis(delta_nm,omega_nm,omega0,gamma)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
xmin=25000
xmax=33000
ymax=450
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
#ax.plot(ffreq, fluo)
ax.plot(ufreq, uvvis)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
#ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.savefig('/home/diegoa/dev/specpy/anthracene-uvvis-fluo.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)

plt.show()
