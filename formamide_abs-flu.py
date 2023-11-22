
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
from ramanlib import fluorescence,uvvis,read_table,uvvis_multi_electronic_state,read_table2

##############################################################################
#  CAM-B3LYP/6-31G*
##############################################################################
system = 'formamide/'
method = 'cam-b3lyp-tz-solv/'
method = 'cam-b3lyp-631/'
method = 'pbe-631/'
method = 'pbe-631-solv/'
method = 'pbe-dzvp/'
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/'+system+method+'dispgrad.dat')

ex_lambdas = np.array([192.,200.,209.,218.,229.,252.,325.,407.])

trans_dip  = read_table2('/home/diegoa/dev/specpy/'+system+method+'trans-dip.dat')

omega0  = read_table2('/home/diegoa/dev/specpy/'+system+method+'omega0.dat')
omega0  = np.multiply(omega0,8065.540107) # convert to cm-1
omega0  = omega0.tolist()

gamma    = [400.,400.,400.,400.,400.,100.]

#omega0=10**7/218.74
states_list = np.identity(len(omega_nm))
states_list = states_list.tolist()
#ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uuvvis =  uvvis_multi_electronic_state(delta_nm,omega_nm,omega0,gamma,trans_dip)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,6))
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
plt.savefig('formamide_uvvis_ev.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
plt.show()



fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,5))
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

ax.set_ylabel(r'ABSORPTION')

plt.savefig('formamide_uvvis_nm.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()





##############################################################################
#  CAM-B3LYP/aug-ccpvtz
##############################################################################

delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements.dat')

gamma      = [800.]*2
gamma      = [800.]
#omega0     = [5.76,7.52,7.99] # in units of eV
#omega0     = [5.76,7.83] ##,6.17,6.45] # in units of eV
#omega0     = np.asarray(omega0)
#omega0     = omega0*8065.540107 # convert to cm-1

#trans_dip=np.asarray([0.0076,0.4453/5.,0.0139])
#trans_dip=[0.0013,0.0258,0.0080,0.0995,0.0007,0.1071,0.3443,0.0192,0.0002]
trans_dip=[0.0013,0.3443]
trans_dip=[0.0013]
trans_dip=[0.0013*3.]

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,7.8286]
omega0     = [5.6578]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1

#omega0=10**7/218.74
states_list = np.identity(len(omega_nm))
states_list.tolist()
#ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uuvvis =  uvvis_multi_electronic_state(delta_nm,omega_nm,omega0,gamma,trans_dip)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,6))
xmax=9.
xmin=5.
ymax=400
ymin=0
plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
#ax.plot(ffreq, fluo)
#ax.plot(ufreq/8065.540107 , uuvvis)
ax.plot(ufreq , uuvvis)
ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)
#ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax.set_xlabel(r'$\omega_I$ (eV)')
#ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.savefig('formamide_uvvis_ev.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
plt.show()



gamma      = [800.]*2
gamma      = [800.]
#omega0     = [5.76,7.52,7.99] # in units of eV
#omega0     = [5.76,7.83] ##,6.17,6.45] # in units of eV
#omega0     = np.asarray(omega0)
#omega0     = omega0*8065.540107 # convert to cm-1

#trans_dip=np.asarray([0.0076,0.4453/5.,0.0139])
#trans_dip=[0.0013,0.0258,0.0080,0.0995,0.0007,0.1071,0.3443,0.0192,0.0002]
trans_dip=[0.0013,0.3443]
trans_dip=[0.3443]

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,7.8286]
omega0     = [7.8286]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1

#omega0=10**7/218.74
states_list = np.identity(len(omega_nm))
states_list.tolist()
#ffreq,fluo = fluorescence(delta_nm,omega_nm,omega0,gamma)
ufreq,uuvvis2 =  uvvis_multi_electronic_state(delta_nm,omega_nm,omega0,gamma,trans_dip)


fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,6))
xmax=9.
xmin=5.
ymax=2.
ymin=0
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
#ax.plot(ffreq, fluo)
ax.plot(ufreq/8065.540107 , uuvvis)
ax.plot(ufreq/8065.540107 , uuvvis2)
ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)
#ax.set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
ax.set_xlabel(r'$\omega_I$ (eV)')
#ax.set_aspect((xmax-xmin)/(ymax-ymin),'box')

ax.set_ylabel(r'ABSORPTION/FLUORESCENCE INTENSITY')
plt.savefig('formamide_uvvis_ev.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
plt.show()









fig1,ax = plt.subplots(ncols=1,nrows=1,figsize=(8,5))
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

ax.set_ylabel(r'ABSORPTION')

plt.savefig('formamide_uvvis_nm.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()
