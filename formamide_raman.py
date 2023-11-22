
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

from ramanlib import raman_ex_disp,read_table,read_table2


##############################################################################
#     RAMAN SHIFT SPECTRA
##############################################################################

##############################################################################
#    ALL STATES SOLVENT
##############################################################################
system = 'formamide/'
system = 'formamide/dimer/'
method = 'cam-b3lyp-tz-solv/'
method = 'cam-b3lyp-631/'
method = 'pbe-631/'
method = 'pbe-dzvp/'
method = 'pbe-631-solv/'
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/'+system+method+'dispgrad.dat')
#omega0 = read_table2('/home/diegoa/dev/specpy/'+system+method+'elec-trans-ev.dat')
omega0 = read_table2('/home/diegoa/dev/specpy/'+system+method+'omega0.dat')
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
trans_dip = read_table2('/home/diegoa/dev/specpy/'+system+method+'trans-dip.dat')
ex_lambdas = np.array([192.,200.,209.,218.,229.,252.,325.,407.])
trans_dip = trans_dip[0]
omega0 = omega0[0]
gamma      = [400.]*10

#omega0     = [ 5.7717 ,  7.0116 , 7.0883 , 7.4425 , 7.7327 , 7.7938 , 8.0559 , 8.0613 ]
#trans_dip  = [0.0050 , 0.0151 , 0.0669 , 1.1231 , 0.7241 , 0.0029 , 0.6840 , 0.1480 ]

# for i in [1,2,4]:
#     delta_nm = np.delete(      delta_nm   , i  ,  0)
#     omega0 = np.delete(        omega0     , i  ,  0)
#     trans_dip = np.delete(     trans_dip  , i  ,  0)
#     gamma = np.delete(         gamma      , i  ,  0)


exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

states_list = []

re,im,ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)
re1,im1,ex_spec_list1,freq,spec_array1 = raman_ex_disp(np.array([delta_nm[0]]),omega_nm,[omega0[0]],[gamma[0]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[0]])
re2,im2,ex_spec_list2,freq,spec_array2 = raman_ex_disp(np.array(delta_nm[1:9]),omega_nm,omega0[1:9],gamma[1:9],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip[1:9])






# PLOTTING
# ------------------------------------------------------------------------------

xmin=45650.
xmax=45700.
xmin=20000.
xmax=60000.
#ymax=5000
ymin=0
fig1,ax = plt.subplots(ncols=3,nrows=30,figsize=(10,38),sharex=True, gridspec_kw={'hspace': 0})
cnt=1
for i in range(30):
    ax[i,0].plot(ex_spec_list[0], ex_spec_list[cnt])
    ax[i,1].plot(re[0], np.absolute(re[cnt]))
    ax[i,1].plot(re[0], np.absolute(re1[cnt]))
    ax[i,1].plot(re[0], np.absolute(re2[cnt]))

    ax[i,2].plot(im[0], np.absolute(im[cnt]))
    ax[i,2].plot(im[0], np.absolute(im1[cnt]))
    ax[i,2].plot(im[0], np.absolute(im2[cnt]))
    #    ymax=np.amax(ex_spec_list[cnt])
    #    ymax=ymax+ymax/10
    ax[i,0].set_xlim(xmin,xmax)
    ax[i,1].set_xlim(xmin,xmax)
    ax[i,2].set_xlim(xmin,xmax)

    idx1 = (np.abs(re[0] - xmin )).argmin()
    idx2 = (np.abs(re[0] - xmax )).argmin()
    max1=np.amax(np.absolute(re[cnt][idx2:idx1]))
    max2=np.amax(np.absolute(re1[cnt][idx2:idx1]))
    max3=np.amax(np.absolute(re2[cnt][idx2:idx1]))
    caca=[max1,max2,max3]
    ymax=np.amax(caca)
    #ymax=ymax+ymax/10
    ymin=0.
    ax[i,1].set_ylim(ymin,ymax)
    ax[i,1].axhline(0, color='black', lw=0.5)

    max1=np.amax(np.absolute(im[cnt][idx2:idx1]))
    max2=np.amax(np.absolute(im1[cnt][idx2:idx1]))
    max3=np.amax(np.absolute(im2[cnt][idx2:idx1]))
    caca=[max1,max2,max3]
    ymax=np.amax(caca)
    #ymax=ymax+ymax/10
    ymin=0.
    ax[i,2].set_ylim(ymin,ymax)
    ax[i,2].axhline(0, color='black', lw=0.5)
    ax[i,0].text(500,1*ymax/2,str(int(omega_nm[i])), fontsize=11, horizontalalignment='right', verticalalignment='center')


    #    ax[i,j].set_ylim(ymin,ymax)
    ax[i,0].set_yscale("log")
        #ax[i,j].set_aspect((xmax-xmin)/(ymax-ymin),'box')
#        ax[i,j].set_ylabel(r' INTENSITY')
    cnt+=1
ax[i,1].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('formamide_dimer_excitation.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)

plt.show()



print(ex_lambdas)




fig0,ax = plt.subplots(ncols=1,nrows=len(ex_lambdas),figsize=(25,15),sharex=True, gridspec_kw={'hspace': 0})
xmax=1750.
xmin=500.
ymin=0.
idx_v8 = (np.abs(freq - 1589.)).argmin()

extra_st = []
extra_w  = []
if len(states_list) > 0:
    for i in range(len(states_list)):
        st=[-1]*len(omega_nm)
        for j in range(0,len(states_list[i]),1):
            st[states_list[i][j]-2]=states_list[i][j+0]
        extra_st.append(st)
    extra_st = np.asarray(extra_st)

    for k in range(len(states_list)):
        wtmp = sum(extra_st[k]*omega_nm)
        extra_w.append(wtmp)
    omega_nm = np.append(omega_nm,extra_w,axis=0)

ticks =  np.array(omega_nm)
for i in range(len(ex_lambdas)):
    ymax=np.amax(spec_array[i])
    ymax=spec_array[i][idx_v8]
    ymax=ymax+ymax/10.
    ax[i].plot(freq, spec_array[i])
    ax[i].set_xlim(xmin,xmax)
    ax[i].set_ylim(ymin,ymax)
    ax[i].set_yticks([])
    ax[i].set_xticks(ticks)
    ax[i].set_xticklabels(ticks.astype(int), rotation=74)
    ax[i].set_aspect((xmax-xmin)/((ymax-ymin)*5),'box')
    ax[i].grid(True, color='darkgrey', linestyle='-', linewidth=1.)
    ax[i].text(500,1*ymax/2,str(int(10**7/ex_lambdas[i])), fontsize=17, horizontalalignment='right', verticalalignment='center')

ax[i].set_xlabel(r'$\omega_I$ (cm$^{-2}$)')
plt.savefig('formamide_raman.svg', dpi=299, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.0,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()







##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################






##############################################################################
#    ONLY STATES 1 AND 7
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements-not-grad.dat')

ex_lambdas = np.array([192.,200.,209.,218.,10**7/45685,229.,252.,325.,407.])
#print(ex_lambdas)
ex_lambdas = 10**7/ex_lambdas
#ex_lambdas = np.flip(ex_lambdas,axis=None)
#print(ex_lambdas)
#gamma      = [800.,600.,160.] # 6-31G*
#gamma      = [500.]*9
n=5.
m=1.

gamma      = [n*609./m,n*609.]
gamma      = [6000.,800.]
#gamma      = [n*609.,n*609./5.]

states_list = np.identity(len(omega_nm))
states_list.tolist()

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,7.82860]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
omega0     = [39700.,51600.]

states_list = np.identity(len(omega_nm))
states_list.tolist()
exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]
states_list = []
#trans_dip=[0.0076,0.4453/5.,0.0139] # cam-b3lyp/6-31G*
trans_dip=[0.0013*128.,0.3443]
re,im,ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)
re1,im1,ex_spec_list1,freq,spec_array1 = raman_ex_disp(np.array([delta_nm[0]]),omega_nm,[omega0[0]],[gamma[0]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[0]])
re2,im2,ex_spec_list2,freq,spec_array2 = raman_ex_disp(np.array([delta_nm[1]]),omega_nm,[omega0[1]],[gamma[1]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[1]])










##############################################################
#  PLOTTING
##############################################################



print(delta_nm)






xmin=45650.
xmax=45700.
xmin=20000.
xmax=55000.
#ymax=5000
ymin=0
fig1,ax = plt.subplots(ncols=3,nrows=12,figsize=(10,18),sharex=True, gridspec_kw={'hspace': 0})
cnt=1
for i in range(12):
    ax[i,0].plot(ex_spec_list[0], ex_spec_list[cnt])
    ax[i,1].plot(re[0], np.absolute(re[cnt]))
    ax[i,1].plot(re[0], np.absolute(re1[cnt]))
    ax[i,1].plot(re[0], np.absolute(re2[cnt]))

    ax[i,2].plot(im[0], np.absolute(im[cnt]))
    ax[i,2].plot(im[0], np.absolute(im1[cnt]))
    ax[i,2].plot(im[0], np.absolute(im2[cnt]))
    #    ymax=np.amax(ex_spec_list[cnt])
    #    ymax=ymax+ymax/10
    ax[i,0].set_xlim(xmin,xmax)
    ax[i,1].set_xlim(xmin,xmax)
    ax[i,2].set_xlim(xmin,xmax)

    idx1 = (np.abs(re[0] - xmin )).argmin()
    idx2 = (np.abs(re[0] - xmax )).argmin()
    max1=np.amax(np.absolute(re[cnt][idx2:idx1]))
    max2=np.amax(np.absolute(re1[cnt][idx2:idx1]))
    max3=np.amax(np.absolute(re2[cnt][idx2:idx1]))
    caca=[max1,max2,max3]
    ymax=np.amax(caca)
    #ymax=ymax+ymax/10
    ymin=0.
    ax[i,1].set_ylim(ymin,ymax)
    ax[i,1].axhline(0, color='black', lw=0.5)

    max1=np.amax(np.absolute(im[cnt][idx2:idx1]))
    max2=np.amax(np.absolute(im1[cnt][idx2:idx1]))
    max3=np.amax(np.absolute(im2[cnt][idx2:idx1]))
    caca=[max1,max2,max3]
    ymax=np.amax(caca)
    #ymax=ymax+ymax/10
    ymin=0.
    ax[i,2].set_ylim(ymin,ymax)
    ax[i,2].axhline(0, color='black', lw=0.5)


    #    ax[i,j].set_ylim(ymin,ymax)
    ax[i,0].set_yscale("log")
        #ax[i,j].set_aspect((xmax-xmin)/(ymax-ymin),'box')
#        ax[i,j].set_ylabel(r' INTENSITY')
    cnt+=1
ax[i,1].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('formamide_excitation.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)

plt.show()








fig0,ax = plt.subplots(ncols=1,nrows=len(ex_lambdas),figsize=(25,15),sharex=True, gridspec_kw={'hspace': 0})
xmax=3999.
xmin=0.
ymin=0.
idx_v8 = (np.abs(freq - 1818.)).argmin()

extra_st = []
extra_w  = []
if len(states_list) > 0:
    for i in range(len(states_list)):
        st=[-1]*len(omega_nm)
        for j in range(0,len(states_list[i]),1):
            st[states_list[i][j]-2]=states_list[i][j+0]
        extra_st.append(st)
    extra_st = np.asarray(extra_st)

    for k in range(len(states_list)):
        wtmp = sum(extra_st[k]*omega_nm)
        extra_w.append(wtmp)
    omega_nm = np.append(omega_nm,extra_w,axis=0)

ticks =  np.array(omega_nm)
for i in range(len(ex_lambdas)):
    ymax=np.amax(spec_array[i])
    ymax=spec_array[i][idx_v8]
    ymax=ymax+ymax/10.
    ax[i].plot(freq, spec_array[i])
    ax[i].set_xlim(xmin,xmax)
    ax[i].set_ylim(ymin,ymax)
    ax[i].set_yticks([])
    ax[i].set_xticks(ticks)
    ax[i].set_xticklabels(ticks.astype(int), rotation=74)
    ax[i].set_aspect((xmax-xmin)/((ymax-ymin)*5),'box')
    ax[i].grid(True, color='darkgrey', linestyle='-', linewidth=1.)
    ax[i].text(500,1*ymax/2,str(int(10**7/ex_lambdas[i])), fontsize=17, horizontalalignment='right', verticalalignment='center')

ax[i].set_xlabel(r'$\omega_I$ (cm$^{-2}$)')
plt.savefig('formamide_raman.svg', dpi=299, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.0,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()






##############################################################################
#    ONLY STATES 1 AND 7 INVERTED DISPLACEMENTS FOR STATES 1 AND 7
delta_nm, omega_nm = read_table('/home/diegoa/dev/specpy/formamide/cam-b3lyp-tz/displacements.dat')

ex_lambdas = np.array([192.,200.,209.,218.,10**7/42200.,229.,252.,325.,407.])
#print(ex_lambdas)
ex_lambdas = 10**7/ex_lambdas
#ex_lambdas = np.flip(ex_lambdas,axis=None)
#print(ex_lambdas)
#gamma      = [800.,600.,160.] # 6-31G*
#gamma      = [500.]*9
n=5.
m=1.

gamma      = [n*609./m,n*609.]
gamma      = [6000.,800.]
#gamma      = [n*609.,n*609./5.]

states_list = np.identity(len(omega_nm))
states_list.tolist()

#omega0     = [5.76,8.0,8.4] # in of eV 6-31G*
omega0     = [5.6578,7.82860]
omega0     = np.multiply(omega0,8065.540107) # convert to cm-1
omega0     = [39700.,51600.]

states_list = np.identity(len(omega_nm))
states_list.tolist()
exmax=90000.
exmin=20000.
ex_range = [exmin,exmax]
sc_range = [0.,4000.]
gamma_scat = 10.

#states_list = [[1,1,9,1],[4,2],[4,1,9,1],[9,2],[9,2,1,1]]
states_list = []
#trans_dip=[0.0076,0.4453/5.,0.0139] # cam-b3lyp/6-31G*
trans_dip=[0.0013*128.,0.3443]
re,im,ex_spec_list,freq,spec_array = raman_ex_disp(delta_nm,omega_nm,omega0,gamma,states_list,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip)
re1,im1,ex_spec_list1,freq,spec_array1 = raman_ex_disp(np.array([delta_nm[0]]),omega_nm,[omega0[0]],[gamma[0]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[0]])
re2,im2,ex_spec_list2,freq,spec_array2 = raman_ex_disp(np.array([delta_nm[1]]),omega_nm,[omega0[1]],[gamma[1]],states_list,gamma_scat,sc_range,ex_range,ex_lambdas,[trans_dip[1]])










##############################################################
#  PLOTTING
##############################################################










xmin=42000.
xmax=43000.
xmin=20000.
xmax=55000.
#ymax=5000
ymin=0
fig1,ax = plt.subplots(ncols=3,nrows=12,figsize=(10,18),sharex=True, gridspec_kw={'hspace': 0})
cnt=1
for i in range(12):
    ax[i,0].plot(ex_spec_list[0], ex_spec_list[cnt])
    ax[i,1].plot(re[0], np.absolute(re[cnt]))
    ax[i,1].plot(re[0], np.absolute(re1[cnt]))
    ax[i,1].plot(re[0], np.absolute(re2[cnt]))

    ax[i,2].plot(im[0], np.absolute(im[cnt]))
    ax[i,2].plot(im[0], np.absolute(im1[cnt]))
    ax[i,2].plot(im[0], np.absolute(im2[cnt]))
    #    ymax=np.amax(ex_spec_list[cnt])
    #    ymax=ymax+ymax/10
    ax[i,0].set_xlim(xmin,xmax)
    ax[i,1].set_xlim(xmin,xmax)
    ax[i,2].set_xlim(xmin,xmax)

    idx1 = (np.abs(re[0] - xmin )).argmin()
    idx2 = (np.abs(re[0] - xmax )).argmin()
    max1=np.amax(np.absolute(re[cnt][idx2:idx1]))
    max2=np.amax(np.absolute(re1[cnt][idx2:idx1]))
    max3=np.amax(np.absolute(re2[cnt][idx2:idx1]))
    caca=[max1,max2,max3]
    ymax=np.amax(caca)
    #ymax=ymax+ymax/10
    ymin=0.
    ax[i,1].set_ylim(ymin,ymax)
    ax[i,1].axhline(0, color='black', lw=0.5)

    max1=np.amax(np.absolute(im[cnt][idx2:idx1]))
    max2=np.amax(np.absolute(im1[cnt][idx2:idx1]))
    max3=np.amax(np.absolute(im2[cnt][idx2:idx1]))
    caca=[max1,max2,max3]
    ymax=np.amax(caca)
    #ymax=ymax+ymax/10
    ymin=0.
    ax[i,2].set_ylim(ymin,ymax)
    ax[i,2].axhline(0, color='black', lw=0.5)


    #    ax[i,j].set_ylim(ymin,ymax)
    ax[i,0].set_yscale("log")
        #ax[i,j].set_aspect((xmax-xmin)/(ymax-ymin),'box')
#        ax[i,j].set_ylabel(r' INTENSITY')
    cnt+=1
ax[i,1].set_xlabel(r'$\omega_I$ (cm$^{-1}$)')
plt.savefig('formamide_excitation.svg', dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        metadata=None)

plt.show()








fig0,ax = plt.subplots(ncols=1,nrows=len(ex_lambdas),figsize=(25,15),sharex=True, gridspec_kw={'hspace': 0})
xmax=3999.
xmin=0.
ymin=0.
idx_v8 = (np.abs(freq - 1283.)).argmin()

extra_st = []
extra_w  = []
if len(states_list) > 0:
    for i in range(len(states_list)):
        st=[0]*len(omega_nm)
        for j in range(0,len(states_list[i]),1):
            st[states_list[i][j]-2]=states_list[i][j+0]
        extra_st.append(st)
    extra_st = np.asarray(extra_st)

    for k in range(len(states_list)):
        wtmp = sum(extra_st[k]*omega_nm)
        extra_w.append(wtmp)
    omega_nm = np.append(omega_nm,extra_w,axis=0)

ticks =  np.array(omega_nm)
for i in range(len(ex_lambdas)):
    ymax=np.amax(spec_array[i])
    ymax=spec_array[i][idx_v8]
    ymax=ymax+ymax/10.
    ax[i].plot(freq, spec_array[i])
    ax[i].set_xlim(xmin,xmax)
    ax[i].set_ylim(ymin,ymax)
    ax[i].set_yticks([])
    ax[i].set_xticks(ticks)
    ax[i].set_xticklabels(ticks.astype(int), rotation=74)
    ax[i].set_aspect((xmax-xmin)/((ymax-ymin)*5),'box')
    ax[i].grid(True, color='darkgrey', linestyle='-', linewidth=1.)
    ax[i].text(500,1*ymax/2,str(int(10**7/ex_lambdas[i])), fontsize=17, horizontalalignment='right', verticalalignment='center')

ax[i].set_xlabel(r'$\omega_I$ (cm$^{-2}$)')
plt.savefig('formamide_raman.svg', dpi=299, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='svg',
        transparent=True, bbox_inches=None, pad_inches=0.0,
        metadata=None)
#    ax[i].set_ylabel(r'RAMAN INTENSITY')
plt.show()








#######################
