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
    colums=tsn1.astype(float)
    displacements = np.asarray(colums[0:-1])
    frequencies = np.asarray(colums[len(colums)-1])
    return ( displacements, frequencies )



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

def uvvis_multi_electronic_state(delta_nm,omega_nm_in,omega0_in,gamma_in,trans_dip):
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
    dt = 1e-2/np.amax(omega_nm)
    period = pi2/np.amin(omega_nm)
    tmax = 10.*period

#    ovlap=np.array([])
    ovlap=[]

    t=0.0
    s=0.5*delta_nm**2
    print('starting time series...')
    nstep=0
    while (t < tmax):
        Inm0 = 0.0
        for j in range(len(omega0_in)):
            arg  = 0.0
            for i in range(len(omega_nm)):
                arg1 = 1.0j * omega_nm[i] * t
                arg2 = -s[(j,i)] * ( 1.0 - cmth.exp(-arg1) )
                arg  = arg + arg2
            arg0 = 1.0j * omega0[j] * t
            arg  = arg - gamma[j] * t
            Inm0 = Inm0 + trans_dip[j] * cmth.exp( arg - arg0 )
        ovlap.append(Inm0)
        #ovlap = np.append(ovlap,[Inm0],axis=0)
        t = t + dt
        if (nstep%20000==0):
            print('Done with step ',nstep,'out of ',int(tmax/dt))
        nstep = nstep + 1
        if (nstep%(int(tmax/dt)-10)==0):
            print('LAST STEP')
    print('END OF TIME SERIES')
#    ovlap = np.array(ovlap)
    print('computing Fourier Transform')
    spec = fft(ovlap)
    #n=ovlap.size
    n=len(ovlap)
    freq = np.arange(-n/2,n/2)*2*np.pi/(np.amax(tmax))
    freq = -fftshift(freq)

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
    dt = 1e-2/np.amax(omega_nm)
    period = pi2/np.amin(omega_nm)
    tmax = 50*period

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
    dt = 1e-2/np.amax(omega_nm)
    period = pi2/np.amin(omega_nm)
    tmax = 50*period

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

def raman_amplitude_mes(delta_nm,omega_nm_in,omega0_in,gamma_in,state,trans_dip):
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
    dt = 1e-2/np.amax(omega_nm_in)
    period = pi2/np.amin(omega_nm)
    tmax = 50.*period

    ovlap=[]

    print(state)
    t=0.0
    s=np.multiply(0.5,np.square(delta_nm))
    fac = [1.,1.,0.5,1./6.,1./24.,1./120.]
    print('Starting time series...')
    nstep=0
    while (t < tmax):
        Inm0 = 0.0
        for j in range(len(omega0)):
            arg=0.0
            exterm = 1.0
            for i in range(len(state)):
                arg1 = 1.0j * omega_nm[i] * t
                arg2 = -s[(j,i)] * ( 1.0 - cmth.exp(-arg1))
                arg  = arg + arg2
                if state[i] != 0:
                    exterm = exterm * cmth.sqrt( fac[int(state[i])] ) * ( - delta_nm[(j,i)] * cmth.sqrt(0.5) * ( 1 - cmth.exp(- arg1) ) )**int(state[i])
            arg0 = 1.0j * omega0[j] * t
            arg = arg - arg0 - gamma[j] * t
            Inm0 = Inm0 + trans_dip[j] * cmth.exp(arg) * exterm ## DANGER!!
        ovlap.append(Inm0)
        t = t + dt
        if (nstep%20000==0):
            print('Done with step ',nstep,'out of ',int(tmax/dt))
        nstep = nstep + 1
        if (nstep%(int(tmax/dt)-10)==0):
            print('LAST STEP')
#    ovlap = np.array(ovlap)
    print('computing Fourier Transform')

    n=len(ovlap)
    #ovlap=([0.0]*n)+ovlap
    n=len(ovlap)
    #freq = np.arange(-n/2,n/2)*2*np.pi/(2*tmax)
    freq = np.arange(-n/2,n/2)*2*np.pi/(tmax)
    freq = -fftshift(freq)

    #spec = abs(fft(ovlap)*np.sqrt(inv2pi) )**2
    spec = fft(ovlap)*np.sqrt(inv2pi)
    plt.plot(freq,np.square(np.absolute(spec)))
    plt.plot(freq,np.absolute(spec))
    plt.plot(freq,np.real(spec))
    plt.plot(freq,np.imag(spec))
    plt.xlim(30000,80000)
    #plt.yscale('log')
    plt.show()
    return freq,spec

def raman_amplitude(delta_nm,omega_nm_in,omega0_in,gamma_in,state):
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

    #spec = abs(fft(ovlap)*np.sqrt(inv2pi) )**2
    spec = fft(ovlap)*np.sqrt(inv2pi)

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

def raman_ex_disp(delta_nm,omega_nm_in,omega0_in,gamma_in,states_extra,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip):
    state = np.identity(len(omega_nm_in),dtype=int)
    # ADDITIONAL STATES (NOT FUNDAMENTAL TRANSITIONS) TO COMPUTE
    if (len(states_extra)>0):
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

    lineshape     = []
    exlambda_idxs = []
    ex_spectrums  = []

    # Loop over the rest of vibrational transitions.
    for i in range(0,len(state)):
    #for i in range(1,len(state)): # DANGER!!!
        print('Computing Raman intensities for transition to state = ', state[i])
        # Compute square Raman amplitude (excitation spectrum) for the first
        # vibrational transition in the list.
        # Following loop sums Raman amplitudes (complex). At the end it is squared.
        ex_ax,ex_ls = raman_amplitude_mes(delta_nm,omega_nm_in,omega0_in,gamma_in,state[i],trans_dip)

        # Excitation spectrum as calculated by the Fourier transform is too big and
        # most of it has zero intensity, so we must trim it to only excitation
        # frequencies of interest.
        tex_ax, tex_ls = trim_spec(ex_range[1], ex_range[0], ex_ax, ex_ls)

        ex_ax = None
        ex_ls = None

        # Squaring absolute value of Raman amplitude.
        tex_ls = np.square(np.absolute(tex_ls))

        # Get the closest indices to the target excitation wavelengths.
        if i == 0:
            for lex in ex_lambdas:
                exlambda_idxs.append(np.abs(lex - tex_ax).argmin())
            # Create an x-axis for the dispersion Raman spectra.
            sc_ax = np.arange(sc_range[0],sc_range[1],0.2)

        # Append excitation spectra to list
        if i == 0: ex_spectrums.append(tex_ax)
        ex_spectrums.append(tex_ls)

        # Compute lorentzian functions for each vibrational transition over such x-axis.
        lorentz = np.square( np.add( sc_ax , -omega_nm_in[i] ) )
        lorentz = np.add( lorentz , gamma_scat**2 )
        lorentz = np.divide( 1.0 , lorentz )
        lorentz = np.multiply( gamma_scat, lorentz )

        # Multiply lorentzian by the square of the Raman amplitude and accumulate in
        # lineshape for each excitation frequency of interest.
        lineshape_tmp = []
        for idx in exlambda_idxs:
            lineshape_tmp.append( np.multiply( lorentz , tex_ls[ idx ] ) )

        # Accumulate with previous vibrational transitions.
        if (i == 0):
            lineshape = lineshape_tmp.copy()
        else:
            lineshape = np.add(lineshape , lineshape_tmp)
        lineshape_tmp = None

    print('FINNISHED SPECTRA CALCULATION!!!')
    return ex_spectrums,sc_ax,lineshape
    #
    # a=[1+1j,2+2j,3+3j,4+4j]
    # b=[5,6,7,8]
    # print(np.add(a,-1+1j))
    # print(np.multiply(2.0,np.divide(1.0,b)))

def raman_ex_disp_test(delta_nm,omega_nm_in,omega0_in,gamma_in,states_extra,gamma_scat,sc_range,ex_range,ex_lambdas,trans_dip):
    state = np.identity(len(omega_nm_in),dtype=int)
    # ADDITIONAL STATES (NOT FUNDAMENTAL TRANSITIONS) TO COMPUTE
    if (len(states_extra)>0):
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

    lineshape     = []
    exlambda_idxs = []
    ex_spectrums  = []

    # Loop over the rest of vibrational transitions.
    for i in range(0,len(state)):
    #for i in range(1,len(state)): # DANGER!!!
        print('Computing Raman intensities for transition to state = ', state[i])
        # Compute square Raman amplitude (excitation spectrum) for the first
        # vibrational transition in the list.
        # Following loop sums Raman amplitudes (complex). At the end it is squared.
        for j in range(lent(omega0_in)):
            ex_ax,ex_ls = raman_amplitude_mes(delta_nm,omega_nm_in,[omega0_in[j]],[gamma_in[j]],state[i],[trans_dip[j]])

            # Excitation spectrum as calculated by the Fourier transform is too big and
            # most of it has zero intensity, so we must trim it to only excitation
            # frequencies of interest.
            tex_ax, tex_ls = trim_spec(ex_range[1], ex_range[0], ex_ax, ex_ls)

            ex_ax = None
            ex_ls = None

            # Squaring absolute value of Raman amplitude.
        tex_ls = np.square(np.absolute(tex_ls))

        # Get the closest indices to the target excitation wavelengths.
        if i == 0:
            for lex in ex_lambdas:
                exlambda_idxs.append(np.abs(lex - tex_ax).argmin())
            # Create an x-axis for the dispersion Raman spectra.
            sc_ax = np.arange(sc_range[0],sc_range[1],0.2)

        # Append excitation spectra to list
        if i == 0: ex_spectrums.append(tex_ax)
        ex_spectrums.append(tex_ls)

        # Compute lorentzian functions for each vibrational transition over such x-axis.
        lorentz = np.square( np.add( sc_ax , -omega_nm_in[i] ) )
        lorentz = np.add( lorentz , gamma_scat**2 )
        lorentz = np.divide( 1.0 , lorentz )
        lorentz = np.multiply( gamma_scat, lorentz )

        # Multiply lorentzian by the square of the Raman amplitude and accumulate in
        # lineshape for each excitation frequency of interest.
        lineshape_tmp = []
        for idx in exlambda_idxs:
            lineshape_tmp.append( np.multiply( lorentz , tex_ls[ idx ] ) )

        # Accumulate with previous vibrational transitions.
        if (i == 0):
            lineshape = lineshape_tmp.copy()
        else:
            lineshape = np.add(lineshape , lineshape_tmp)
        lineshape_tmp = None

    print('FINNISHED SPECTRA CALCULATION!!!')
    return ex_spectrums,sc_ax,lineshape
