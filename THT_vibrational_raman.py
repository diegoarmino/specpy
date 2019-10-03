
import numpy as np
import math as mth
import cmath as cmth
import scipy.constants as cnst
from scipy.fftpack import fft, ifft,  fftfreq, fftshift,dct,idct
import matplotlib.pyplot as plt

from ramanlib import raman_ex_disp


##############################################################################
#     RAMAN SHIFT SPECTRA
##############################################################################

omega_nm   = np.array([354.,444.,934.,1192.,1285.,1295.,1403.,1581.,1635.])
delta_nm   = np.array([0.55,0.23,0.23,0.82,0.485,0.02,0.085,0.38,1.32])
ex_lambdas = np.array([35651., 38805., 39651., 39809., 39952.,40323.,41034.,41425.,42123.,42644.])
ex_lambdas = np.flip(ex_lambdas)
gamma      = 160.
omega0     = 39805.

states_list = np.identity(9)
states_list.tolist()

exmax=48000.
exmin=38000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]

ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas)




##############################################################
#  PLOTTING
##############################################################


fig1,ax = plt.subplots(ncols=1,nrows=len(ex_lambdas),figsize=(16,16),sharex=True, gridspec_kw={'hspace': 0})
xmax=4000.
xmin=0.
ymin=0.
idx_v9 = (np.abs(freq - 1635.)).argmin()

omega_nm   = np.array([354.,444.,934.,1192.,1285.,1295.,1403.,1581.,1635.])
extra_st = []
extra_w  = []
for i in range(len(states_list)):
    st=[0]*len(omega_nm)
    for j in range(0,len(states_list[i]),2):
        st[states_list[i][j]-1]=states_list[i][j+1]
    extra_st.append(st)
extra_st = np.asarray(extra_st)

for k in range(len(states_list)):
    wtmp = sum(extra_st[k]*omega_nm)
    extra_w.append(wtmp)
omega_nm = np.append(omega_nm,extra_w,axis=0)

ticks =  np.array([354.,444.,934.,1192.,1295.,1403.,1635.,1989., 2384., 2827., 3270., 3624.])

for i in range(len(ex_lambdas)):
    ymax=np.amax(spec_array[i])
    ymax=spec_array[i][idx_v9]
    ymax=ymax+ymax
    ax[i].plot(freq, spec_array[i])
    ax[i].set_xlim(xmin,xmax)
    ax[i].set_ylim(ymin,ymax)
    ax[i].set_yticks([])
    ax[i].set_xticks(ticks)
    ax[i].set_xticklabels(ticks, rotation=75)
    ax[i].set_aspect((xmax-xmin)/((ymax-ymin)*8),'box')
    ax[i].grid(True, color='darkgrey', linestyle='-', linewidth=0.5)
    ax[i].text(500,2*ymax/3,str(ex_lambdas[i]), fontsize=18, horizontalalignment='right', verticalalignment='center')

ax[i].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('tht_raman.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()


xmax=exmax
xmin=exmin
#ymax=5000
ymin=0
fig1,ax = plt.subplots(ncols=3,nrows=3,figsize=(8,8),sharex=True, gridspec_kw={'hspace': 0})
cnt=1
for i in range(3):
    for j in range(3):
        ax[i,j].plot(ex_spec_list[0], ex_spec_list[cnt])
        ymax=np.amax(ex_spec_list[cnt])
        ymax=ymax+ymax/10
        ax[i,j].set_xlim(xmin,xmax)
        ax[i,j].set_ylim(ymin,ymax)
        ax[i,j].set_aspect((xmax-xmin)/(ymax-ymin),'box')
#        ax[i,j].set_ylabel(r' INTENSITY')
        cnt+=1
ax[i,j].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('tht_excitation.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)

plt.show()
