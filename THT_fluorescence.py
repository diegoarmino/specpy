
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
omega_nm = np.array([354.,444.,934.,1192.,1285.,1295.,1403.,1581.,1635.])
delta_nm = np.array([0.55,0.23,0.23,0.82,0.485,0.02,0.085,0.38,1.32])
gamma    = 160.
omega0=39805
states_list = np.identity(9)
states_list.tolist()
ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uuvvis = uvvis(delta_nm,omega_nm,omega0,gamma)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
xmax=52000
xmin=28000
ymax=100
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
ax.plot(ffreq, fluo)
ax.plot(ufreq, uuvvis)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.show()

##############################################################################
#  BP86 TZVP PARAMETRS
##############################################################################
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/disp.dat')
gamma    = 160.
omega0=39805
states_list = np.identity(len(omega_nm))
states_list.tolist()
ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uuvvis = uvvis(delta_nm,omega_nm,omega0,gamma)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,8))
xmax=52000
xmin=28000
ymax=400
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
ax.plot(ffreq, fluo)
ax.plot(ufreq, uuvvis)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.show()
