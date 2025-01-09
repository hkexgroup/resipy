#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:10:43 2025
Custom sequence generator 
@author: jimmy
"""
import os, time, warnings 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

def customprint(text,end='\n'):
    logfile = 'commandGen.log'
    if not os.path.exists(logfile):
        fh = open(logfile,'w')
        fh.close()
    fh = open(logfile,'a')
    fh.write(text)
    fh.write(end)
    fh.close()
    print(text,end=end)

def printgeninfo(array='DP-DP',nele=128,amin=1,amax=5,nmin=1,nmax=5,reverse=False):
    customprint('\n=========================')
    customprint('Creating custom sequence:')
    customprint('array type = %s'%array)
    customprint('Number of electrodes = %i'%nele)
    customprint('amin = %i, amax = %i'%(amin,amax))
    customprint('nmin = %i, nmax = %i'%(nmin,nmax))
    order = 'Forward'
    if reverse:
        order = 'Reverse'
    customprint('order = %s'%order)
    customprint('Generating ...')

def printseqinfo(createtime, ncmdline, nmeas, cycletime=0.8):
    customprint('Creation time = %f seconds'%createtime)
    customprint('CMD lines = %i'%ncmdline)
    customprint('Number of measurements = %i'%nmeas)
    # totaltime = (ncmdline * cycletime )/60 # ??? estimated duration in minutes 
    # cant figure out how sting estimates measurement time, using emperical relation
    # from sting software puts us at ~9.55 meas per minute
    totaltime = 0.105 * ncmdline
    customprint('Estimated duration = %f minutes'%totaltime)
    
class Generator():
    def __init__(self, elec):
        self.elec = elec
        self.seq = [] 
        
    def clear(self):
        self.seq = [] 
        
    def ddgen(self, amin=1, amax=5, nmin=1, nmax=8):
        
        nele = len(self.elec)
        labels = self.elec.label.values 
        
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmax,int)
        assert isinstance(nmax,int)
        assert isinstance(nele,int)
        
        e1 = 0 # start electrode index 
        en = nele # end electrde index 
        step = 1 # step
    
        nmeas = 0 # number of measurements 
        seq = [] # list to store single channel commands 
    
        # create lists of possible a and n values 
        da = amax - amin
        dn = nmax - nmin 
        ia = [amin + i for i in range(da+1)] # possible a values
        jn = [nmin + i for i in range(dn+1)] # possible n values
    
        t0 = time.time()
            
        for i in range(e1,en,step):
            c2 = i # current negative 
            for a in ia:
                c1 = c2 + a 
                if c1 >= nele or c1 < 1:
                    continue
                for n in jn:
                    p1 = c1 + (n*a)
                    p2 = p1 + a
                    if p1 >= nele or p1 < 1:
                        continue
                    elif p2 >= nele or p2 < 1:
                        continue
                    else:
                        seq.append([labels[c1],labels[c2],labels[p1],labels[p2]])
                        nmeas += 1
    
        t1 = time.time()
        dt = t1 - t0
        
        self.seq += seq 
                            
        return np.array(seq)
    
    def wsgen(self, amin=1, amax=5, nmin=1, nmax=8):
    
        nele = len(self.elec)
        labels = self.elec.label 
        
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmax,int)
        assert isinstance(nmax,int)
        assert isinstance(nele,int)
        
        e1 = 0 # start electrode index 
        en = nele # end electrde index 
        step = 1 # step
    
        nmeas = 0 # number of measurements 
        seq = [] # list to store single channel commands 
    
        # create lists of possible a and n values 
        da = amax - amin
        dn = nmax - nmin 
        ia = [amin + i for i in range(da+1)] # possible a values
        jn = [nmin + i for i in range(dn+1)] # possible n values
    
        t0 = time.time() # time the loop 
        
        for i in range(e1,en,step):
            c1 = i # current negative 
            for a in ia:
                for n in jn:
                    c2 = c1 + (2*n*a) + a 
                    p1 = c1 + (n*a)
                    p2 = p1 + a
                    if c2 >= nele:
                        continue
                        # skip onwards becuase last electrode is beyond end of line 
                    elif p1 < c1 or p2 > c2:
                        continue
                        # skip onwards as potential electrodes must be inside current electrodes 
                    else:
                        seq.append([labels[c1],labels[c2],labels[p1],labels[p2]])
                        nmeas += 1
    

        t1 = time.time()
        dt = t1 - t0
        
        self.seq += seq 
                            
        return np.array(seq)
    
    def wgen(self, amin=1, amax=5, nmin=None, nmax=None):
    
        nele = len(self.elec)
        labels = self.elec.label.values 
        
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nele,int)
        
        e1 = 0 # start electrode index 
        en = nele # end electrde index 
        step = 1 # step
    
        nmeas = 0 # number of measurements 
        seq = [] # list to store single channel commands 
    
        # create lists of possible a and n values 
        da = amax - amin
        ia = [amin + i for i in range(da+1)] # possible a values
    
        t0 = time.time() # time the loop 
        
        for i in range(e1,en,step):
            c1 = i # current negative 
            for a in ia:
                c2 = c1 + (2*a) + a 
                p1 = c1 + a
                p2 = p1 + a
                if c2 >= nele:
                    continue
                    # skip onwards becuase last electrode is beyond end of line 
                elif p1 < c1 or p2 > c2:
                    continue
                    # skip onwards as potential electrodes must be inside current electrodes 
                else:
                    seq.append([labels[c1],labels[c2],labels[p1],labels[p2]])
                    nmeas += 1
    
        t1 = time.time()
        dt = t1 - t0
        
        self.seq += seq 
                            
        return np.array(seq)
    
    
    def xdgen(self, amin=1, amax=5, nmin=1, nmax=8, discrepancy = 0.5):
        ux = np.unique(self.elec.x)
        ux = sorted(ux)
        
        # loop through x coordinates find the next borehole/line in x direction 
        seq = []
        for i, x0 in enumerate(ux):
            if (i + 1) >= len(ux):
                continue 
        
            bindx = x0 == self.elec.x 
            belec0 = self.elec[bindx].reset_index(drop=True)
            
            # find electrodes for the next borehole / line 
            bindx = ux[i+1] == self.elec.x
            belec1 = self.elec[bindx].reset_index(drop=True)
            
            # continue if not enough canidate x coordinates to create cross lines from 
            if len(belec0) < 2 or len(belec1) < 2:
                continue 
            
            # find elevation of first borehole electrodes 
            nelec0 = len(belec0)
            nelec1 = len(belec1)
            for axis in ['y','z']:
                if all(belec1[axis] == belec1[axis][0]):
                    continue 
                
                for k in range(nelec0):
                    if k >= nelec1:
                        continue 
                    # find index of electrode on the same plane (roughly)
                    z0 = belec0[axis][k]
                    zmin = z0 - discrepancy
                    zmax = z0 + discrepancy
                    bindx = (belec1[axis] > zmin) & (belec1[axis] < zmax) 
                    elec = belec1[bindx]
                    if len(elec) == 0:
                        continue 
                    i0 = k 
                    i1 = belec1[bindx].index[0]
                    
                    for a in range(amax):
                        # extend in the positive direction 
                        c1 = i0 + a 
                        c2 = i1 + a 
                        for j in range(amin,amax):
                            p1 = c1 + j 
                            if p1 < 0 or p1 >= nelec0:
                                continue 
                            p2 = c2 + j 
                            if p2 < 0 or p2 >= nelec1:
                                continue 
    
                            line = [belec0.label[c1], belec1.label[c2], belec0.label[p1], belec1.label[p2]]
                            if line not in seq: 
                                seq.append(line)
                        
                        # extend in the NEGATIVE direction 
                        c1 = i0 - a 
                        c2 = i1 - a 
                        for j in range(amin,amax):
                            p1 = c1 - j 
                            if p1 < 0 or p1 >= nelec0:
                                continue 
                            p2 = c2 - j 
                            if p2 < 0 or p2 >= nelec1:
                                continue 
    
                            line = [belec0.label[c1], belec1.label[c2], belec0.label[p1], belec1.label[p2]]
                            if line not in seq: 
                                seq.append(line)
                                
        self.seq += seq 
        
        return np.array(seq)
    
    
    def mggen(self, amin=1, amax=5, nmin=1, nmax=8, smin=1, smax=2):
        ''' Genetrate measurement matrix for multigradient array Torleif Dahlin.
        
        Parameters
        ----------
        elec_num : int
            Number of electrodes
        a : int
            Spacing between potential electrodes (in electrode spacing).
        n : int
            Multiplier for `a` to determine spacing from A to M.
        s : int
            Seperation factor for current electrodes, should be the intermediate
            numbers.
        '''
        nele = len(self.elec)
        labels = self.elec.label.values 
        
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmax,int)
        assert isinstance(nmax,int)
        assert isinstance(nele,int)
        
        e1 = 0 # start electrode index 
        en = nele # end electrde index 
        step = 1 # step
        
        seq = []
        
        # create lists of possible a and n values 
        da = amax - amin
        dn = nmax - nmin 
        ds = smax - smin 
        ia = [amin + i for i in range(da+1)] # possible a values
        jn = [nmin + i for i in range(dn+1)] # possible n values
        ks = [smin + i for i in range(ds+1)]
        
        for i in range(e1, en, step): 
            c1 = i
            for a in ia: 
                for n in jn: 
                    for s in ks: 
                        c2 = c1 + s
                        p1 = c1 + (n*a)
                        p2 = p1 + a 
                        c2 = p2 + (s-((n-1)*a))
                        
                        if p1 >= nele or p1 < 0:
                            continue 
                        if p2 >= nele or p2 < 0:
                            continue 
                        if c2 >= nele or c2 < 0: 
                            continue 
                        seq.append([labels[c1],labels[c2],labels[p1],labels[p2]])
                        
        self.seq += seq 
            
        return np.array(seq)
    
#%% convert labels to integers  
def seq2int(seq):
    """
    Convert sequence into integers. 

    Parameters
    ----------
    seq : nd array 
        N by 4 matrix of str

    Returns
    -------
    iseq : nd array 
        N by 4 matrix of ints 

    """
    # check for string numbers in the labels 
    hasString = False 
    if len(seq[0,0].split()) == 2:
        hasString = True 
    cacheString = {}
    iseq = np.zeros_like(seq, dtype=int)
    
    if hasString: 
        # search through and find unique electrode numbers associated with each string 
        for i in range(seq.shape[0]):
            for j in range(4):
                label = seq[i,j]
                string = int(label.split()[0])
                number = int(label.split()[1])
                if string in cacheString.keys():
                    cacheString[string].append(number)
                else:
                    cacheString[string] = [number]
        # create unique index for each electrode 
        lookup = {}
        count = 1 
        for string in sorted(cacheString.keys()):
            lookup[string] = {}
            for n in np.unique(cacheString[string]):
                lookup[string][n] = count 
                count += 1 
        # reassign sequence with unique values 
        for i in range(seq.shape[0]):
            for j in range(4):
                label = seq[i,j]
                string = int(label.split()[0])
                number = int(label.split()[1])
                iseq[i,j] = lookup[string][number] 

    else:
        for i in range(seq.shape[0]):
            for j in range(4):
                iseq[i,j] = int(seq[i,j])
    
    return iseq 
        
#%% add reciprocal pairs 
def reciprocalise(seq):
    """
    Add reciprocals to a single channel sequence 

    Parameters
    ----------
    seq : nd array 
        N by 4 matrix of ints 

    Returns
    -------
    nseq : nd array 
        sequence with reciprocal measurements appended. 

    """
    rseq = np.zeros_like(seq, dtype=int)
    
    # to make a reciprocal 
    # p1 >> c2 and p2 >> c1 or... 
    # n >> b and m >> a 
    rseq[:,0] = seq[:,3]
    rseq[:,1] = seq[:,2]
    rseq[:,2] = seq[:,1]
    rseq[:,3] = seq[:,0]
    
    nseq = np.vstack([seq,rseq])
    return nseq 
    
    
#%% multichannelise a sequence 
def multichannelise(elecdf, seq, maxcha=8): 
    """
    Compress single channel sequence into a more efficient sequence that can 
    leverage modern multichannel instruments. 

    Parameters
    ----------
    elec : pd.DataFrame
        Electrode datafame 
    seq : nd array 
        DESCRIPTION.
    maxcha : int, optional
        Maximum number of channels. The default is 8.

    Returns
    -------
    lines : list
        Lines that can be written to a multichannel sequence 

    """
    elec = elecdf[['x','y','z']].values 
    ### first determine if measurements are nested ###
    #find mid points of AB 
    array = seq-1 # array for indexing electrodes 
    ABmid = (elec[array[:,0]] + elec[array[:,1]]) / 2 # mid points of AB 
    MNmid = (elec[array[:,2]] + elec[array[:,3]]) / 2 # mid points of MN 
    
    ABrad = np.sqrt(np.sum((elec[array[:,0]] - ABmid)**2,axis=1)) # radius of AB circle 
    MNrad = np.sqrt(np.sum((elec[array[:,2]] - MNmid)**2,axis=1)) # radius of MN circle 
    
    A2mn = np.sqrt(np.sum((elec[array[:,0]] - MNmid)**2,axis=1)) # distance of A to mid point of MN 
    B2mn = np.sqrt(np.sum((elec[array[:,1]] - MNmid)**2,axis=1)) # distance of B to mid point of MN 
    N2ab = np.sqrt(np.sum((elec[array[:,2]] - ABmid)**2,axis=1)) # distance of N to mid point of AB 
    M2ab = np.sqrt(np.sum((elec[array[:,3]] - ABmid)**2,axis=1)) # distance of M to mid point of AB
    
    iABinMN = (A2mn < MNrad) & (B2mn < MNrad)
    iMNinAB = (N2ab < ABrad) & (M2ab < ABrad)
    inested = iABinMN | iMNinAB #if AB encompasses MN or MN encompasses AB 
    
    # note unested measurements can be multichannelised but may have 
    # some excess measurements with extreme geometric factors 
    
    # work out maximum channel digit 
    maxperline = maxcha + 2 
    cmax = np.max(seq)
    ndig = len(str(cmax))
    nmeas = array.shape[0]
    cpairidn = [0]*nmeas # nested pair ids 
    cpairidu = [0]*nmeas # unested pair ids 
    idtemplate= '{:0>%id}{:0>%id}'%(ndig, ndig)
    for i in range(nmeas):
        if inested[i]: # dont attempt to multichannelise nested measurements 
            cpairidn[i] = int(idtemplate.format(seq[i,0], seq[i,1]))
        else: 
            cpairidu[i] = int(idtemplate.format(seq[i,0], seq[i,1]))
        
    uniquepairsu = np.unique(cpairidu)
    uniquepairsn = np.unique(cpairidn)
    lines = []
    for uid in uniquepairsu: 
        if uid == 0:
            continue 
        index = uid == cpairidu # index all instances of unique c1 and c2 being used for nested measurements 
        subset = seq[index] # subset of sequence 
        subset = subset[np.argsort(subset[:,2])] # sort by p1 electrode 
        # if subset.shape[0]>3:
        lines.append([])
        c1, c2 = subset[0,0], subset[0,1]
        lines[-1].append(c1)
        lines[-1].append(c2) 
        for j in range(subset.shape[0]):
            p1 = subset[j,2]
            p2 = subset[j,3]
            if len(lines[-1]) > (maxperline):
                lines.append([])
                lines[-1].append(c1)
                lines[-1].append(c2) 
            if p1 not in lines[-1]: 
                lines[-1].append(p1)
            if len(lines[-1]) > (maxperline):
                lines.append([])
                lines[-1].append(c1)
                lines[-1].append(c2) 
                lines[-1].append(p1)
            if p2 not in lines[-1]: 
                lines[-1].append(p2)
    
    # add nested measurements 
    for uid in uniquepairsn: 
        if uid == 0:
            continue 
        index = uid == cpairidn # index all instances of unique c1 and c2 being used for nested measurements 
        subset = seq[index] # subset of sequence 
        lines.append([])
        c1, c2 = subset[0,0], subset[0,1]
        lines[-1].append(c1)
        lines[-1].append(c2) 
        for j in range(subset.shape[0]):
            p1 = subset[j,2]
            p2 = subset[j,3]
            if len(lines[-1]) > (maxperline-1):
                lines.append([])
                lines[-1].append(c1)
                lines[-1].append(c2)  
            if p1 in lines[-1] or p2 in lines[-1]:
                lines.append([])
                lines[-1].append(c1)
                lines[-1].append(c2)  
            
            lines[-1].append(p1)
            lines[-1].append(p2)
         
    # append 0s to each line
    for line in lines: 
        while len(line) < (maxperline+1):
            line.append(0)
            
    return lines 
    
#%% condition a sequence for IP reducing effects 
def checkCmdSep(lines, limit=None): 
    # limit is needed for efficiency 
    ncmdline = len(lines)
    if limit is None: 
        limit = ncmdline + 1  
    cmdsep = [ncmdline]*ncmdline
    for i in range(ncmdline):
        cline = lines[i]
        c1 = cline[0]
        c2 = cline[1]
        # loop through and grab next line 
        count = 0 
        for j in range(i, ncmdline):
            if j <= i:
                continue 
            nline = lines[j]
            pn = nline[2:]
            if c1 in pn or c2 in pn: 
                cmdsep[i] = (j - i)
                break 
            count += 1 
            # use counter to make the function more efficient in loops 
            # as we're only looking for minimum seperations 
            if count > limit:
                break 
    return cmdsep 

def conditionSequence(lines, plot=False):
    """
    Condition a multichannel sequence to avoid reusing potential electrodes 
    after they have just been used for injecting current. This avoids IP effects. 

    Parameters
    ----------
    lines : list 
        lines for mulitchannel measurements 

    Returns
    -------
    condLines : list 
        conditioned lines that minimise IP effects on measurements, improving
        reliability. 

    """
    # condition sequence to avoid using current electrodes
    # immiediately after measurement
    nswap = 0 # total number of swaps for all values of n, i, j etc... 
    ncmdline = len(lines)
    condLines = [lines[i] for i in range(ncmdline)]
    minsep= int(0.1*ncmdline)
    if minsep > 32:
        minsep = 32 
    condcmdsep = checkCmdSep(condLines , minsep)
    loop = 0 
    iswap = 9e9 
    iswapcache = -1 
    while min(condcmdsep) < minsep:
        loop += 1 
        print('Resorting iteration %i...'%loop)
        iswapcache = iswap 
        iswap = 0 # number of swaps for all values in loop 
        for i in range(ncmdline-1):
            for n in range(1,minsep):
                c1c = condLines[i][0] # current c1 electrode 
                c2c = condLines[i][1] # current c2 electrode
                j = i+n # shuffle index 
                if j >= ncmdline:
                    break
                reshuffle = True
                while reshuffle:
                    pn = condLines[j][2:] # next p electrodes
                    
                    # c1c or c2c cannot be in PN! 
                    if c1c in pn or c2c in pn:
                        # need to reshuffle order!
                        reshuffle = True
                        if j < (ncmdline-1):
                            j += 1
                        else:
                            j = 0 
                            break 
                    else:
                        reshuffle = False
                
                if j != (i+n): # consective line needs shuffling in this case
                    swapline = condLines[i+n]
                    condLines[i+n] = condLines[j]
                    condLines[j] = swapline
                    nswap += 1 
                    iswap += 1 
                    
        print('Total swaps = %i'%nswap) 
        print('Iteration swaps = %i'%iswap) 
        condcmdsep = checkCmdSep(condLines, minsep)

        if iswapcache <= iswap:# and loop>10: 
            print('No improvement in number of swaps over last 2 iterations, treating as sorted...')
            break 
        
    if plot: 
        cmdsep = checkCmdSep(lines)
        fig, (ax0,ax1) = plt.subplots(nrows=2)
        ax1.set_xlabel('Command Seperation')
        ax0.set_ylabel('Count')
        maxbin = np.max(np.unique(cmdsep)[:-1])
        ax0.hist(cmdsep, bins=np.arange(maxbin), density=False, label='Unconditioned', color='r')
        ax0.legend()
        ax0.set_title('Unconditioned versus Conditioned CMD Seperations')
        
        condcmdsep = checkCmdSep(condLines)
        ax1.hist(condcmdsep, bins=np.arange(maxbin), density=False, label='Conditioned')
        ax1.set_ylim(ax0.get_ylim())
        ax1.legend()
    
        usep, ucon = np.unique(condcmdsep,return_counts=True)
        for ax in (ax0, ax1):
            ax.set_ylim([0,max(ucon)])
        
    return condLines

#%% write to file 
def write2csv(lines, fname ):
    ncmdline = len(lines)
    fh = open(fname, 'w')
    nperline = len(lines[0])
    header = 'C1, C2, '
    i = 1 
    while len(header.split(',')) < nperline: 
        header += 'P%i, '%i
        i+=1 
    
    header += 'ChannelIDs\n'
    fh.write(header)
    
    for line in lines: 
        channels = '' 
        count = 0 
        for n in line: 
            fh.write('%i, '%n)
            if n > 0: 
                count += 1
        for i in range(count-3):
            channels += '%i'%(i+1)
        
        fh.write(channels)
        fh.write('\n')
        
    fh.close() 
    
#%% export sequence 
def exportSequence(fname, sequence, elec, integer = True, reciprocals = True, 
                   multichannel = True,  condition = True, maxcha = 8): 
    """
    Export single channel sequence (ie for forward modelling) into something 
    usable by a resistivity instrument 

    Parameters
    ----------
    fname : str
        Export file path 
    sequence : nd array 
        N by 4 array of electrode configurations 
    elec : pd.DataFrame
        Electrode data frame 
    integer : bool, optional
        Flag to convert sequence into integers before export. The default is True.
    reciprocals : bool, optional
        Flag to add reciprocal measurements. The default is True.
    multichannel : bool, optional
        Flag to convert measurements to a multichannel. The default is True.
    condition : bool, optional
        Flag to condition measurements for avoiding IP effects. The default is True.
    maxcha: int, optional 
        Maximum number of active channels of the resistivity instrument (normally
        8 for modern instruments). 

    """
    if not fname.endswith('.csv'):
        fname += '.csv'
        
    if integer: 
        sequence = seq2int(sequence)
        
    if reciprocalise: 
        sequence = reciprocalise(sequence)
        
    if multichannel: 
        sequence = multichannelise(elec, sequence, maxcha=maxcha)
    
    if condition:
        sequence = conditionSequence(sequence)
        
    if isinstance(sequence, list):
        write2csv(sequence, fname)
    else:
        tmp = {
            'a':sequence[:,0], 
            'b':sequence[:,1],
            'm':sequence[:,2],
            'n':sequence[:,3],
            }
        pd.DataFrame(tmp).to_csv(fname, index=False)
    
            
#%% test 
elecx = np.linspace(0,23,24)
elecy = np.zeros_like(elecx)
elecz = np.zeros_like(elecx)
elecl = ['%i'%(i+1) for i in range(len(elecx))]

elec = pd.DataFrame({
    'x':elecx, 
    'y':elecy, 
    'z':elecz, 
    'label':elecl, 
    })

dpdpParam = {
    'amin': 1,
    'amax': 8,
    'nmin': 1,
    'nmax': 8,
    }

# wenner schlumberger command gen
# parameters a and n etc 
wsParam = {
    'amin': 5,
    'amax': 20,
    'nmin': 5,
    'nmax': 20,
    }

gen = Generator(elec)

_ = gen.ddgen(**dpdpParam)
_ = gen.wsgen(**wsParam)

seq = np.array(gen.seq) 

iseq = seq2int(seq)
rseq = reciprocalise(iseq)     
mseq = multichannelise(elec, rseq)
cseq = conditionSequence(mseq, True)
# exportSequence('sequence.csv', seq, elec)

    
    
