
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

from ramanlib import raman_ex_disp,read_table


##############################################################################
#     RAMAN SHIFT SPECTRA
##############################################################################
##############################################################################
#     ALL ELECTRONIC STATES
##############################################################################
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements.dat')

ex_lambdas = np.array([192.,200.,209.,218.,229.,252.,325.,407.])
ex_lambdas = 10**7/ex_lambdas
#ex_lambdas = np.flip(ex_lambdas,axis=None)
print(ex_lambdas)
# [52083 50000. 47846. 45871. 43668. 39682. 30769. 24570.]
#gamma      = [800.,600.,160.] # 6-31G*
#gamma      = [500.]*9
gamma      = [800.]*9

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,6.7266,6.7506,7.2737,7.4984,7.6410,7.8286,8.0887,8.0980]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
states_list = np.identity(len(omega_nm))
states_list.tolist()

exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,6.7266,6.7506,7.2737,7.4984,7.6410,7.8286,8.0887,8.0980]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
states_list = np.identity(len(omega_nm))
states_list.tolist()
print(omega0)
exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]
states_list = []
#trans_dip=[0.0076,0.4453/5.,0.0139] # cam-b3lyp/6-31G*
trans_dip=[0.0013,0.0258,0.0080,0.0995,0.0007,0.1071,0.3443,0.0192,0.0002]

ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)






##############################################################################
#    ONLY STATES 1 AND 7
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements.dat')

ex_lambdas = np.array([192.,200.,209.,218.,229.,10**7/41400.,252.,325.,407.])
ex_lambdas = 10**7/ex_lambdas
#ex_lambdas = np.flip(ex_lambdas,axis=None)
print(ex_lambdas)
#gamma      = [800.,600.,160.] # 6-31G*
#gamma      = [500.]*9
n=5.
m=1.

gamma      = [n*609./m,n*609.]
gamma      = [800.,800.]
#gamma      = [n*609.,n*609./5.]

states_list = np.identity(len(omega_nm))
states_list.tolist()

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,7.82860]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
omega0     = [39700.,51600.]
print(np.divide(10**7,omega0))

states_list = np.identity(len(omega_nm))
states_list.tolist()
print(omega0)
exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]
states_list = []
#trans_dip=[0.0076,0.4453/5.,0.0139] # cam-b3lyp/6-31G*
trans_dip=[0.3443/200.,0.3443]
trans_dip=[0.3443,0.3443]

ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)
ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,[omega0[0]],[gamma[0]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[0]])
ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,[omega0[1]],[gamma[1]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[1]])










##############################################################################
#    ONLY STATE 1
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements.dat')

ex_lambdas = np.array([192.,200.,209.,218.,229.,252.,325.,407.])
ex_lambdas = 10**7/ex_lambdas
#ex_lambdas = np.flip(ex_lambdas,axis=None)

#gamma      = [800.,600.,160.] # 6-31G*
#gamma      = [500.]*9
n=5.
m=10.

gamma      = [n*609./m]
#gamma      = [n*609.,n*609./5.]

states_list = np.identity(len(omega_nm))
states_list.tolist()

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
omega0     = [39700.]
print(np.divide(10**7,omega0))

states_list = np.identity(len(omega_nm))
states_list.tolist()
print(omega0)
exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]
states_list = []
#trans_dip=[0.0076,0.4453/5.,0.0139] # cam-b3lyp/6-31G*
trans_dip=[0.3443/200.]

ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)






##############################################################################
#    ONLY STATE 7
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements.dat')

ex_lambdas = np.array([192.,200.,209.,218.,229.,252.,325.,407.])
ex_lambdas = 10**7/ex_lambdas
#ex_lambdas = np.flip(ex_lambdas,axis=None)

#gamma      = [800.,600.,160.] # 6-31G*
#gamma      = [500.]*9
n=5.
m=10.

gamma      = [n*609.]
#gamma      = [n*609.,n*609./5.]

states_list = np.identity(len(omega_nm))
states_list.tolist()

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [7.82860]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
omega0     = [51600.]
print(np.divide(10**7,omega0))

states_list = np.identity(len(omega_nm))
states_list.tolist()
print(omega0)
exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]
states_list = []
#trans_dip=[0.0076,0.4453/5.,0.0139] # cam-b3lyp/6-31G*
trans_dip=[0.3443]

ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)



##############################################################
#  PLOTTING
##############################################################


fig1,ax = plt.subplots(ncols=1,nrows=len(ex_lambdas),figsize=(16,16),sharex=True, gridspec_kw={'hspace': 0})
xmax=4000.
xmin=0.
ymin=0.
idx_v9 = (np.abs(freq - 1283.)).argmin()

extra_st = []
extra_w  = []
if len(states_list) > 0:
    print('caca')
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

ticks =  np.array(omega_nm)

for i in range(len(ex_lambdas)):
    ymax=np.amax(spec_array[i])
    ymax=spec_array[i][idx_v9]
    ymax=ymax+ymax
    ax[i].plot(freq, spec_array[i])
    ax[i].set_xlim(xmin,xmax)
    ax[i].set_ylim(ymin,ymax)
    ax[i].set_yticks([])
    ax[i].set_xticks(ticks)
    ax[i].set_xticklabels(ticks.astype(int), rotation=75)
    ax[i].set_aspect((xmax-xmin)/((ymax-ymin)*5),'box')
    ax[i].grid(True, color='darkgrey', linestyle='-', linewidth=0.5)
    ax[i].text(500,2*ymax/3,str(int(10**7/ex_lambdas[i])), fontsize=18, horizontalalignment='right', verticalalignment='center')

ax[i].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('formamide_raman.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()












xmax=exmax
xmin=exmin
xmin=41000.
xmax=42000.
#ymax=5000
ymin=0
fig1,ax = plt.subplots(ncols=3,nrows=4,figsize=(8,8),sharex=True, gridspec_kw={'hspace': 0})
cnt=1
for i in range(4):
    for j in range(3):
        ax[i,j].plot(ex_spec_list[0], ex_spec_list[cnt])
    #    ymax=np.amax(ex_spec_list[cnt])
    #    ymax=ymax+ymax/10
        ax[i,j].set_xlim(xmin,xmax)
    #    ax[i,j].set_ylim(ymin,ymax)
        ax[i,j].set_yscale("log")
        #ax[i,j].set_aspect((xmax-xmin)/(ymax-ymin),'box')
#        ax[i,j].set_ylabel(r' INTENSITY')
        cnt+=1
ax[i,j].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('formamide_excitation.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)

plt.show()
