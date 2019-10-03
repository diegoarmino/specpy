import numpy as np
import math as mth
import cmath as cmth
import scipy.constants as cnst
from scipy.fftpack import fft, ifft,  fftfreq, fftshift,dct,idct
import matplotlib.pyplot as plt



def read_table(path):
    with open(path,mode='r') as f0:
        sn1_file = f0.read()
    f0.closed
    lsn1 = sn1_file.split('\n')
    sn1=[]
    for i in range(len(lsn1)-1):
        sn1.append(lsn1[i].split())
    tsn1=np.transpose(np.array(sn1))
    out=tsn1.astype(float)
    return out



def raman_full(tmax,omega_nm,omega_gm,delta_nm):
    c = cnst.c
    hbar = cnst.hbar
    pi = mth.pi
    fout = open("output.dat", "w")
    while (t < tmax):
        cosvar = cmth.cos(omega_nm * t)
        sinvar = cmth.sin(omega_nm * t)
        fratio = omega_nm/omega_gm
        i_fratio = 1/fratio

        alpha = -0.5j*(1j * cosvar - fratio * sinvar)/(1j*i_fratio*sinvar + cosvar)
        p_tnm = fratio * delta_nm * sinvar
        q_tnm = delta_nm * (1 - cosvar)
        gamma_tnm = 0.5j * cmth.log(1j * i_fratio * sinvar + cosvar) + 0.5 * p_tnm * (q_tnm - delta_nm)

        a0 = -alpha * q_tnm**2 - 1j * p_tnm * q_tnm + 1j * gamma_tnm
        a1 = 2 * alpha * q_tnm + 1j * p_tnm
        a2 = 0.5 - alpha

        Inm0 = cmth.exp( a0 + (  a1**2 / (4 * (1-a2))  ) )/cmth.sqrt(1-a2)
        Inm1 = -cmth.sqrt(2) * (a1 + a1 * a2 / (1-a2)) * Inm0

        print(t,"\t",gamma_tnm.real,"\t",gamma_tnm.imag, file=fout)

        t = t + dt
    fout.close()

def raman_simple(delta_nm,omega_nm_in,omega0_in,gamma_in,state):
    c = cnst.c
    ccmfs = c*100.0/1e15 # speed of light in cm/fs

    omega_nm = omega_nm_in * ccmfs # convert omega form cm-1 to fs-1
    omega0   = omega0_in   * ccmfs # convert to frequency fs-1
    gamma    = gamma_in    * ccmfs

    hbar = cnst.hbar
    pi = mth.pi
    pi2= 2*pi
    inv2pi= 1/pi2

#    dt = 1e-2/np.amax(omega_nm)
    dt = 1e-2/np.amax(omega0)
    period = pi2/np.amin(omega_nm)
    tmax = 10*period

    ovlap=[]

    t=0
    s=0.5*delta_nm**2
    fout = open("output.dat", "w")
    while (t < tmax):
        Inm=1.0
        arg0 = 1j * omega0 * t
        for i in range(len(omega_nm)):
            arg1 = 1j * omega_nm[i] * t
            arg2 = -s[i] * ( 1 - cmth.exp(-arg1)) - arg1 * 0.5 - arg0
            Inm0 = cmth.exp(arg2)

            if (state[i]==8):
                Inm1 = -cmth.sqrt(s[i]) * ( 1 - cmth.exp( -arg1 ) ) * Inm0
                Inm  = Inm * Inm1 * cmth.exp( -gamma * t ) * cmth.exp( arg1*0.5 )
            else:
                Inm  = Inm * Inm0 * cmth.exp( -gamma * t ) * cmth.exp( arg1*0.5 )

        Inm = Inm * cmth.exp( -arg0 )
        ovlap.append(Inm)
        print("%12.4f\t%12.4E"%(t*ccmfs, abs(Inm)),file=fout,sep='\t')
        t = t + dt
    ovlap = np.array(ovlap)
#    time  = np.array(time)
    fout.close()
#    time,correl_func = read_table("output.dat")
    #spec = dct(ovlap)#*cmth.sqrt(cmth.pi/2)
    spec = fft(ovlap)#*cmth.sqrt(cmth.pi/2)
    n=ovlap.size
    freq = fftfreq(n,d=(-dt*ccmfs))
#    freq = freq / ccmfs
    return freq,spec


def uvvis(delta_nm,omega_nm_in,omega0_in,gamma_in):
    c = cnst.c
    ccmfs = c*100.0/1e15 # speed of light in cm/fs

    omega_nm = omega_nm_in #* ccmfs # convert omega form cm-1 to fs-1
    omega0   = omega0_in   #* ccmfs # convert to frequency fs-1
    gamma    = gamma_in    #* ccmfs

    hbar = cnst.hbar
    pi = mth.pi
    pi2= 2*pi
    inv2pi= 1/pi2

#    dt = 1e-2/np.amax(omega_nm)
    dt = 5e-2/np.amax(omega_nm)
    period = pi2/np.amin(omega_nm)
    tmax = 20*period

    ovlap=[]

    t=0.0
    s=0.5*delta_nm**2
    fout = open("output.dat", "w")
    while (t < tmax):
        Inm=1.0
        arg=0.
        for i in range(len(omega_nm)):
            arg1 = 1.0j * omega_nm[i] * t
            arg2 = -s[i] * ( 1.0 - cmth.exp(-arg1))
            arg  = arg + arg2
        arg0 = 1.0j * omega0 * t
        arg = arg - arg0 - gamma * t
        Inm0 = cmth.exp(arg)

        ovlap.append(Inm0)
        print("%12.4f\t%12.4E"%(t*ccmfs, abs(Inm0)),file=fout,sep='\t')
        t = t + dt
    ovlap = np.array(ovlap)
#    time  = np.array(time)
    fout.close()
#    time,correl_func = read_table("output.dat")
    #spec = dct(ovlap)#*cmth.sqrt(cmth.pi/2)
    spec = fft(ovlap)
    n=ovlap.size
    freq = np.arange(-n/2,n/2)*2*np.pi/(np.amax(tmax))
    freq = -fftshift(freq)

    return freq,spec

def fluorescence(delta_nm,omega_nm_in,omega0_in,gamma_in):
    c = cnst.c
    ccmfs = c*100.0/1e15 # speed of light in cm/fs

    omega_nm = omega_nm_in #* ccmfs # convert omega form cm-1 to fs-1
    omega0   = omega0_in   #* ccmfs # convert to frequency fs-1
    gamma    = gamma_in    #* ccmfs

    hbar = cnst.hbar
    pi = mth.pi
    pi2= 2*pi
    inv2pi= 1/pi2

#    dt = 1e-2/np.amax(omega_nm)
    dt = 5e-2/np.amax(omega_nm)
    period = pi2/np.amin(omega_nm)
    tmax = 20*period

    ovlap=[]

    t=0.0
    s=0.5*delta_nm**2
    fout = open("output.dat", "w")
    while (t < tmax):
        Inm=1.0
        arg=0.
        for i in range(len(omega_nm)):
            arg1 = 1.0j * omega_nm[i] * t
            #arg2 = -s[i] * ( 1.0 - cmth.exp(-arg1)) # uvvis
            arg2 = -s[i] * ( 1.0 - cmth.exp(arg1))
            arg  = arg + arg2
        arg0 = 1.0j * omega0 * t
        arg = arg - arg0 - gamma * t
        Inm0 = cmth.exp(arg)

        ovlap.append(Inm0)
        print("%12.4f\t%12.4E"%(t*ccmfs, abs(Inm0)),file=fout,sep='\t')
        t = t + dt
    ovlap = np.array(ovlap)
#    time  = np.array(time)
    fout.close()
#    time,correl_func = read_table("output.dat")
    #spec = dct(ovlap)#*cmth.sqrt(cmth.pi/2)
    spec = fft(ovlap)
    n=ovlap.size
    freq = np.arange(-n/2,n/2)*2*np.pi/(np.amax(tmax))
    freq = -fftshift(freq)

    return freq,spec

def raman(delta_nm,omega_nm_in,omega0_in,gamma_in,state):
    c = cnst.c
    ccmfs = c*100.0/1e15 # speed of light in cm/fs

    omega_nm = omega_nm_in[0:len(delta_nm)] #* ccmfs # convert omega form cm-1 to fs-1
    omega0   = omega0_in   #* ccmfs # convert to frequency fs-1
    gamma    = gamma_in    #* ccmfs

    hbar = cnst.hbar
    pi = mth.pi
    pi2= 2*pi
    inv2pi= 1/pi2

#    dt = 1e-2/np.amax(omega_nm)
    dt = 5e-2/np.amax(omega_nm_in)
    period = pi2/np.amin(omega_nm)
    tmax = 50*period

    ovlap=[]

    t=0.0
    s=0.5*delta_nm**2
    fout = open("output.dat", "w")
    fac = [1.,1.,0.5,1./6.,1./24.,1./120.]
    while (t < tmax):
        Inm=1.0
        arg=0.0
        exterm = 1.0
        for i in range(len(state)):
            arg1 = 1.0j * omega_nm[i] * t
            arg2 = -s[i] * ( 1.0 - cmth.exp(-arg1))
            arg  = arg + arg2
            if state[i] != 0:
                exterm = exterm * fac[int(state[i])] * ( -cmth.sqrt(s[i]) * ( 1 - cmth.exp(- arg1) ) )**int(state[i])
        arg0 = 1.0j * omega0 * t
        arg = arg - arg0 - gamma * t

        Inm0 = cmth.exp(arg) * exterm

        ovlap.append(Inm0)
        print("%12.4f\t%12.4E"%(t*ccmfs, abs(Inm0)),file=fout,sep='\t')
        t = t + dt
    ovlap = np.array(ovlap)
    fout.close()
    ovlap=np.pad(ovlap,(ovlap.size,0),'constant',constant_values=0.) # DANGER DEBUG!
    n=ovlap.size
    freq = np.arange(-n/2,n/2)*2*np.pi/(2*tmax)
    freq = -fftshift(freq)
#    omega_tmp=omega_nm*np.array(state,dtype='float')
#    omega_s = omega0 - sum(omega_tmp)
    spec = abs(fft(ovlap)*np.sqrt(inv2pi) )**2
#    spec = spec * freq * omega_s  #!!DANGER this is the correct form
    #spec = spec * freq * omega_nm[0] * 0.5

    return freq,spec



def trim_spec(min,max,axis,spec):
    min_idx = (np.abs(axis - min)).argmin()
    max_idx = (np.abs(axis - max)).argmin()
    if max_idx < min_idx:
        tmp = min_idx
        min_idx = max_idx
        max_idx = tmp
    trim_ax = axis[min_idx-1:max_idx+1]
    trim_sp = spec[min_idx-1:max_idx+1]
    return trim_ax,trim_sp

def raman_ex_disp(delta_nm,omega_nm_in,omega0_in,gamma_in,states_extra,gamma_scat,sc_range,ex_range,ex_lambdas):
    state = np.identity(len(omega_nm_in),dtype=int)
    # ADDITIONAL STATES (NOT FUNDAMENTAL TRANSITIONS) TO COMPUTE
    extra_st = []
    extra_w  = []
    for i in range(len(states_extra)):
        st=[0]*len(omega_nm_in)
        for j in range(0,len(states_extra[i]),2):
            st[states_extra[i][j]-1]=states_extra[i][j+1]
        extra_st.append(st)
    extra_st = np.asarray(extra_st)
    state    = np.append(state,extra_st,axis=0)

    for k in range(len(states_extra)):
        wtmp = sum(extra_st[k]*omega_nm_in)
        extra_w.append(wtmp)
    omega_nm_in = np.append(omega_nm_in,extra_w,axis=0)

    sc_ax    = np.arange(sc_range[0],sc_range[1],0.2)
    lineshape     = []
    exlambda_idxs = []
    ex_spectrums  = []
    print('Computing Raman intensities for transition to state = ', state[0])
    ex_ax,ex_ls = raman(delta_nm,omega_nm_in,omega0_in,gamma_in,state[0])
    tex_ax,tex_ls = trim_spec(ex_range[1],ex_range[0],ex_ax,ex_ls)
    ex_ax=None
    ex_ls=None
    ex_spectrums.append(tex_ax)
    ex_spectrums.append(tex_ls)

    for lex in ex_lambdas:
        exlambda_idxs.append(np.abs(lex - tex_ax).argmin())
    lorentz = ((sc_ax - omega_nm_in[0])**2+gamma_scat**2)
    lorentz = gamma_scat/lorentz
    for idx in exlambda_idxs:
        lineshape.append( lorentz * tex_ls[ idx ] )
    lineshape = np.asarray(lineshape)

    for i in range(1,len(state)):
        print('Computing Raman intensities for transition to state = ', state[i])
        ex_ax,ex_ls = raman(delta_nm,omega_nm_in,omega0_in,gamma_in,state[i])
        tex_ax,tex_ls = trim_spec(ex_range[0],ex_range[1],ex_ax,ex_ls)
        ex_spectrums.append(tex_ls)
        ex_ax=None
        ex_ls=None
        lorentz = ((sc_ax - omega_nm_in[i])**2+gamma_scat**2)
        lorentz = gamma_scat/lorentz
        lineshape_tmp = []
        for idx in exlambda_idxs:
            lineshape_tmp.append( lorentz * tex_ls[ idx ] )
        lineshape_tmp = np.asarray(lineshape_tmp)
        lineshape = lineshape + lineshape_tmp
        lineshape_tmp = None

    return ex_spectrums,sc_ax,lineshape
