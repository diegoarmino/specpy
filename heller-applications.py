import numpy as np
import math as mth
import cmath as cmth
import scipy.constants as cnst
from scipy.fftpack import fft, ifft,  fftfreq, fftshift,dct,idct
import matplotlib.pyplot as plt

import ramanlib

##############################################################################
#     APPLICATIONS: Fictitious system in Heller and Tannor 1982
##############################################################################

omega_nm = np.array([920.0,1550.0])
delta_nm = np.array([0.9,0.556])
gamma    = 375.0
#freq,spec = uvvis(delta_nm,omega_nm,18500,gamma,state)
fig1,ax = plt.subplots(ncols=3,nrows=1,figsize=(16,4))
states_list = [[0,1],[1,0]]
xmax=25000
xmin=12000
ymax=1.2e12
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
cnt=0
for st in states_list:
    freq,spec = raman(delta_nm,omega_nm,18500,gamma,st,0,18000.)
    ax[cnt].plot(freq, spec)
    ax[cnt].set_xlim(xmin,xmax)
    ax[cnt].set_ylim(ymin,ymax)
    ax[cnt].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
    ax[cnt].set_aspect((xmax-xmin)/(ymax-ymin),'box')
    cnt+=1
freq,spec = uvvis(delta_nm,omega_nm,18500,gamma)
ax[cnt].plot(freq, spec)
ymin=0
ymax=260
ax[cnt].set_xlim(xmin,xmax)
ax[cnt].set_ylim(ymin,ymax)
ax[cnt].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax[cnt].set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax[0].set_ylabel(r'RAMAN INTENSITY $I^{(0,0) \rightarrow (0,1)}$')
ax[1].set_ylabel(r'RAMAN INTENSITY $I^{(0,0) \rightarrow (1,0)}$')
ax[2].set_ylabel(r'ABSORPTION INTENSITY')
plt.show()


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


##############################################################################
#     RAMAN SHIFT SPECTRA
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


##############################################################################
#     RAMAN SHIFT SPECTRA
##############################################################################

omega_nm = np.array([354.,444.,934.,1192.,1285.,1295.,1403.,1581.,1635.])
delta_nm = np.array([0.55,0.23,0.23,0.82,0.485,0.02,0.085,0.38,1.32])
#ex_lambdas = np.array([35651., 38805., 39651., 39809., 39952.,40323.,41034.,41425.,42123.,42644.])
ex_lambdas = np.array([35651.])
gamma    = 160.
omega0=39805
fig1,ax = plt.subplots(ncols=1,nrows=len(ex_lambda),figsize=(8,8))
states_list = np.identity(9)
states_list.tolist()
xmax=48000
xmin=38000
ymax=300
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

ex_range = [xmin,xmax]
sc_range = [0,4000]
gamma_scat = 2.

frec,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas)
for i in range(len(ex_lambda)):
    ax.plot(freq, spec[i])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
    ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')
plt.show()
