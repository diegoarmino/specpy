
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
from ramanlib import fluorescence,uvvis,read_table,uvvis_multi_electronic_state

##############################################################################
#  BP86 TZVP PARAMETRS
##############################################################################
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/b3lyp-631/displacements.dat')
gamma    = [800.,600.,120.]
#omega0     = [5.76,7.52,7.99] # in units of eV
omega0     = [5.76,7.13,8.0] ##,6.17,6.45] # in units of eV
omega0     = np.asarray(omega0)
omega0     = omega0*8065.540107 # convert to cm-1


#omega0=10**7/218.74
states_list = np.identity(len(omega_nm))
states_list.tolist()
#ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uuvvis =  uvvis_multi_electronic_state(delta_nm,omega_nm,omega0,gamma)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
xmax=10.
xmin=4.
ymax=400
ymin=0
plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
#ax.plot(ffreq, fluo)
ax.plot(ufreq/8065.540107 , uuvvis)
ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)
#ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax.set_xlabel(r'$\omega_I$ (eV)')
#ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.show()



fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
xmax=100
xmin=300
ymax=400
ymin=0
#plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
#ax.plot(10**7/ffreq, fluo)
ax.plot(10**7/ufreq, uuvvis)
#ax.set_xlim(xmin,xmax)
ax.set_xlim(xmax,xmin)
#ax.set_ylim(ymin,ymax)
#ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
#ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.show()
