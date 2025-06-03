#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:10:43 2025
Custom sequence generator 
@author: jimmy
"""
import time 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

## condition a sequence for IP reducing effects 
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
    
class Generator():
    def __init__(self, elec, seqIdx = None):
        self.elec = elec
        if seqIdx is None: 
            self.seqIdx = [np.arange(len(self.elec))]
        else:
            self.seqIdx = seqIdx 
        self.seq = np.array([]) 
        self.seqR = np.array([]) # reciprocal sequence 
        self.lines = [] # variable to store multichannel lines 
        self.linesR = [] # reciprocal lines 
        self.runtime = 0 
        
    def copy(self):
        g = Generator(self.elec, self.seqIdx)
        g.seq = self.seq.copy() 
        g.lines = [line for line in self.lines] 
        g.runtime = self.runtime 
        # handle reciprocal 
        g.seqR = self.seqR.copy() 
        g.linesR = [line for line in self.linesR] 
        return g 
        
    def clear(self):
        self.seq = np.array([]) 
        self.seqR = np.array([]) # reciprocal sequence 
        self.lines = [] # variable to store multichannel lines 
        self.linesR = [] # reciprocal lines 
        self.runtime = 0 
        
    def ddgen(self, amin=1, amax=5, nmin=1, nmax=8, si=0, e1=None, en=None):
        
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmax,int)
        assert isinstance(nmin,int)
        
        felec = self.elec.loc[self.seqIdx[si],:]
        nele = len(felec)
        labels = felec.label.values 
        
        if e1 is None: 
            e1 = 0 # start electrode index 
        if en is None: 
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
        self.runtime += dt 
        
        if self.seq.shape[0] > 0: 
            self.seq = np.vstack([self.seq, np.array(seq)])
        else:
            self.seq = np.array(seq)
                            
        return np.array(seq)
    
    def wsgen(self, amin=1, amax=5, nmin=1, nmax=8, si=0, e1=None, en=None):
    
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmin,int)
        assert isinstance(nmax,int)
        
        felec = self.elec.loc[self.seqIdx[si],:]
        nele = len(felec)
        labels = felec.label.values 
        
        if e1 is None: 
            e1 = 0 # start electrode index 
        if en is None: 
            en = nele # end electrde index 
        step = 1 
    
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
                    p1 = c1 + (n*a)
                    p2 = p1 + a
                    c2 = p2 + (n*a)
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
        self.runtime += dt 
        
        if self.seq.shape[0] > 0: 
            self.seq = np.vstack([self.seq, np.array(seq)])
        else:
            self.seq = np.array(seq)
                            
        return np.array(seq)
    
    def wgen(self, amin=1, amax=5, nmin=None, nmax=None, si=0, e1=None, en=None):
    
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        
        felec = self.elec.loc[self.seqIdx[si],:]
        nele = len(felec)
        labels = felec.label.values 
        
        if e1 is None: 
            e1 = 0 # start electrode index 
        if en is None: 
            en = nele # end electrde index 
        step = 1 
    
        nmeas = 0 # number of measurements 
        seq = [] # list to store single channel commands 
    
        # create lists of possible a and n values 
        da = amax - amin
        ia = [amin + i for i in range(da+1)] # possible a values
    
        t0 = time.time() # time the loop 
        
        for i in range(e1,en,step):
            c1 = i # current negative 
            for a in ia:
                p1 = c1 + a
                p2 = p1 + a
                c2 = p2 + a 
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
        
        self.runtime += dt 
        if self.seq.shape[0] > 0: 
            self.seq = np.vstack([self.seq, np.array(seq)])
        else:
            self.seq = np.array(seq)
                            
        return np.array(seq)
    
    
    def xdgen(self, amin=0, amax=8, nmin=1, nmax=3, discrepancy = 0.5, dump=print):
        
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmin,int)
        assert isinstance(nmax,int)
        
        t0 = time.time() 
        
        da = amax - amin
        dn = nmax - nmin 
        ia = [amin + i for i in range(da+1)] # possible a values
        jn = [nmin + i for i in range(dn+1)] # possible n values
        
        felec = self.elec # electrode data frame 
        
        ux = np.unique(felec.x)
        ux = sorted(ux)
        
        # loop through x coordinates find the next borehole/line in x direction 
        seq = []
        for i, x0 in enumerate(ux):
            if (i + 1) >= len(ux):
                continue 
        
            bindx = x0 == felec.x 
            belec0 = felec[bindx].reset_index(drop=True)
            
            # find electrodes for the next borehole / line 
            bindx = ux[i+1] == felec.x
            belec1 = felec[bindx].reset_index(drop=True)
            
            # continue if not enough canidate x coordinates to create cross lines from 
            if len(belec0) < 2 or len(belec1) < 2:
                continue 
            
            # find y coordinates / elevation of first 3d/borehole electrodes 
            nelec0 = len(belec0)
            nelec1 = len(belec1)
            for axis in ['y','z']:
                if all(belec1[axis] == belec1[axis][0]):
                    # skip if no change in axis 
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
                    c1 = i0 
                    
                    for a in ia:
                        # extend in the positive direction 
                        c2 = i1 + a 
                        for n in jn:
                            if c2 < 0 or c2 >= nelec1:
                                continue 
                            p1 = c1 + n 
                            if p1 < 0 or p1 >= nelec0:
                                continue 
                            p2 = c2 + n
                            if p2 < 0 or p2 >= nelec1:
                                continue 
                            line = [belec0.label[c1], belec1.label[c2], belec0.label[p1], belec1.label[p2]]
                            if line not in seq: 
                                seq.append(line)
                        
                        # extend in the NEGATIVE direction 
                        c2 = i1 - a 
                        for n in jn:
                            if c2 < 0 or c2 >= nelec1:
                                continue 
                            p1 = c1 - n 
                            if p1 < 0 or p1 >= nelec0:
                                continue 
                            p2 = c2 - n 
                            if p2 < 0 or p2 >= nelec1:
                                continue 
                            line = [belec0.label[c1], belec1.label[c2], belec0.label[p1], belec1.label[p2]]
                            if line not in seq: 
                                seq.append(line)
                                
        t1 = time.time()
        dt = t1 - t0
        if len(seq) == 0:
            t = 'Warning: Could not create any cross (equatorial) dipole measurements'
            dump(t)
            return np.array([])
        
        self.runtime += dt 
        if self.seq.shape[0] > 0: 
            self.seq = np.vstack([self.seq, np.array(seq)])
        else:
            self.seq = np.array(seq)
        
        return np.array(seq)
    
    
    def mggen(self, amin=1, amax=5, nmin=1, nmax=8, mmin=1, mmax=2, 
              si=0, e1=None, en=None):
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
        assert isinstance(amin,int)
        assert isinstance(amax,int)
        assert isinstance(nmax,int)
        assert isinstance(nmin,int)
        assert isinstance(mmin,int)
        assert isinstance(mmax,int)
        
        t0 = time.time() 
        felec = self.elec.loc[self.seqIdx[si],:]
        nele = len(felec)
        labels = felec.label.values 
        
        if e1 is None: 
            e1 = 0 # start electrode index 
        if en is None: 
            en = nele # end electrde index 
        step = 1 
        
        seq = []
        
        # create lists of possible a and n values 
        da = amax - amin
        dn = nmax - nmin 
        dm = mmax - mmin 
        ia = [amin + i for i in range(da+1)] # possible a values
        jn = [nmin + i for i in range(dn+1)] # possible n values
        km = [mmin + i for i in range(dm+1)]
        
        for i in range(e1, en, step): 
            c1 = i
            for a in ia: 
                for n in jn: 
                    for m in km: 
                        p1 = c1 + (n*a)
                        p2 = p1 + a 
                        # c2 = p2 + (s-((n-1)*a)) -->> formula in resipy image, but does not work! 
                        c2 = p2 + (m*a)
                        
                        if p1 >= nele or p1 < 0:
                            continue 
                        if p2 >= nele or p2 < 0:
                            continue 
                        if c2 >= nele or c2 < 0: 
                            continue 
                        seq.append([labels[c1],labels[c2],labels[p1],labels[p2]])
                        
                        
        t1 = time.time()
        dt = t1 - t0
        
        self.runtime += dt 
        if self.seq.shape[0] > 0: 
            self.seq = np.vstack([self.seq, np.array(seq)])
        else:
            self.seq = np.array(seq)
            
        return np.array(seq)
    
    def fromfile(self, fname):
        if not isinstance(fname, str): 
            raise ValueError('Argument "fname" should be of type "string".')
            
        t0 = time.time() 
        _seq = pd.read_csv(fname, header=0)
        
        if _seq.shape[1] != 4:
            raise ValueError('The file should be a CSV file with headers with exactly 4 columns '
                             '(a, b, m, n) with electrode numbers.')
        
        seq = np.asarray(_seq.values, dtype=str)
        
        # check for string numbers in the labels in the sequence 
        hasDataString = False 
        if len(seq[0,0].split()) == 2:
            hasDataString = True 
            
        # check electrodes have strings 
        hasElecString = False 
        if len(self.elec.label[0].split()) == 2: 
            hasElecString = True 
            
        # 4 scenarios 
        if hasDataString is False and hasElecString is False: 
            # do not need to worry about electrode strings at all in this case 
            # so do not edit the sequence and move on 
            pass 
        elif hasDataString is True and hasElecString is False: 
            # this means the user is importing a 3D sequence for 2D problem most likely --> error situation 
            raise Exception('Incoming data has string numbers but electrode data frame is stringless')
        elif hasDataString is True and hasElecString is True: 
            # the user has imported a 3D sequence and has 3D appropoate electrodes, no edits needed 
            pass 
        elif hasDataString is False and hasElecString is True: 
            # user has imported a 2D sequence but the electrodes are 3D, hence the sequence needs to be doctored 
            for i in range(seq.shape[0]):
                for j in range(4):
                    seq[i,j] = '1 %s'%seq[i,j]
                
        if self.seq.shape[0] > 0: 
            self.seq = np.vstack([self.seq, np.array(seq)])
        else:
            self.seq = np.array(seq)
            
        t1 = time.time()
        dt = t1 - t0
        self.runtime += dt 
        return seq 
        
    def generate(self, params = [('dpdp',1,8,1,8)], dump=print):
        """
        Loop through generation parameters anf generate measurement configurations 
        accordingly. 

        Parameters
        ----------
        params : list, optional
            DESCRIPTION. The default is [('dpdp',1,8,1,8)].

        Returns
        -------
        sequence: nd array
            N by 4 matrix, the measurement schedule for forward modelling 

        """
        self.clear()
        def logtext(p):
            # little helper function for logging progress 
            t='('
            for i in range(len(p)-1):
                t+= '%i, '%p[i]
            t+= '%i)'%p[-1]
            return t 

        for p in params: 
            config = p[0]
            if config in ['custom', 'custSeq']:
                dump('Reading custom sequence...')
                fname = p[1] 
                _  = self.fromfile(fname)
                # skip onwards if custom sequence 
                continue 
            if config in ['cross','xbh', 'xli', 'equat-dp']: 
                dump('Generating cross line/hole: %s'%logtext(p[1:5]))
                # special case where the sequence index can be ignored as all electrodes need to be considered for cross measurements 
                if len(p) == 6: 
                    _ = self.xdgen(p[1], p[2], p[3], p[4], p[5], dump=dump) 
                else: 
                    _ = self.xdgen(p[1], p[2], p[3], p[4], dump=dump)  
                continue 
            for i in range(len(self.seqIdx)):
                if config in ['dipole-dipole', 'dpdp']:
                    dump('Generating dipole-dipole: %s'%logtext(p[1:5]))
                    _ = self.ddgen(p[1], p[2], p[3], p[4], i) 
                elif config in ['wenner', 'w']:
                    dump('Generating wenner: %s'%logtext(p[1:3]))
                    _ = self.wgen(p[1], p[2], i)
                elif config in ['wenner-schlumberger', 'schlum', 'ws']: 
                    dump('Generating wenner-schlum: %s'%logtext(p[1:5]))
                    _ = self.wsgen(p[1], p[2], p[3], p[4], i) 
                elif config in ['multigradient', 'mg', 'multigrad']:
                    dump('Generating multi-grad: %s'%logtext(p[1:7]))
                    _ = self.mggen(p[1], p[2], p[3], p[4], p[5], p[6], i)
                else:
                    dump('Warning, config of type "%s" is not recognised'%config)
                    
        return self.seq.copy() 
    
    def checkString(self):
        # check for string numbers in the labels 
        hasString = False 
        if len(self.seq[0,0].split()) == 2:
            hasString = True 
        return hasString 
         
    ## convert labels to integers  
    def seq2int(self,ignorestring=False):
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
        if self.seq.shape[0] == 0:
            return 
        # check for string numbers in the labels 
        hasString = False 
        hasMoreThan1String = False 
        stringCache = [] 
        if len(self.seq[0,0].split()) == 2:
            hasString = True 
            plabels = self.seq.flatten()
            for label in plabels: 
                string = int(label.split()[0])
                if string not in stringCache: 
                    stringCache.append(string)
            if len(stringCache) > 1: 
                hasMoreThan1String = True 
                
        iseq = np.zeros_like(self.seq, dtype=int)
        template = '{:0>2d}{:0>3d}' 
        stringCache = []
        
        if hasString: 
            # search through and find unique electrode numbers associated with each string 
            for i in range(self.seq.shape[0]):
                for j in range(4):
                    label = self.seq[i,j]
                    string = int(label.split()[0])
                    number = int(label.split()[1])
                    
                    if ignorestring:
                        iseq[i,j] = number 
                    elif hasMoreThan1String:
                        iseq[i,j] = int(template.format(string, number))
                    else:
                        iseq[i,j] = number 
        else:
            for i in range(self.seq.shape[0]):
                for j in range(4):
                    iseq[i,j] = int(self.seq[i,j])
        
        self.seq = iseq 
        
        if self.seqR.shape[0] > 0:
            g = self.copy()
            g.seq = self.seqR.copy()
            g.seqR = np.array([])
            self.seqR = g.seq2int() 
            
        return iseq 
            
    ## add reciprocal pairs 
    def reciprocalise(self, split=False):
        """
        Add reciprocals to a single channel sequence 
    
        Parameters
        ----------
        split: bool
            If True, split the sequence into forward and recirprocal 
            measurements. Else append reciprocals to current measurement 
            sequence. 
    
        Returns
        -------
        nseq : nd array 
            reciprocal sequence if "split" is True, else newly appended sequence
    
        """
        rseq = np.zeros_like(self.seq, dtype=int)
        
        # to make a reciprocal 
        # p1 >> c2 and p2 >> c1 or... 
        # n >> b and m >> a 
        rseq[:,0] = self.seq[:,3]
        rseq[:,1] = self.seq[:,2]
        rseq[:,2] = self.seq[:,1]
        rseq[:,3] = self.seq[:,0]
        
        if split: 
            self.seqR = rseq 
            return rseq  
        else: 
            nseq = np.vstack([self.seq,rseq])
            self.seq = nseq
            return nseq 

        
    ## multichannelise a sequence 
    def multichannelise(self, maxcha=8, aggressive=False): 
        """
        Compress single channel sequence into a more efficient sequence that can 
        leverage modern multichannel instruments. 
    
        Parameters
        ----------
        maxcha : int, optional
            Maximum number of channels. The default is 8.
        agressive: bool, optional
            Don't preserve polarity of current electrodes to more agressively 
            pack the measurements. 
    
        Returns
        -------
        lines : list
            Lines that can be written to a multichannel sequence 
    
        """
        if maxcha < 1: 
            raise ValueError('Number of channels must exceed 1')
        elec = self.elec[['x','y','z']].values 
        ### first determine if measurements are nested ###
        #find mid points of AB 
        hasString = False 
        if isinstance(self.seq[0,0], str): 
            self.seq2int()
            
        array = self.seq-1 
        
        if np.max(array) >= elec.shape[0]:
            # need to account for electrode labels perhaps not being continuous 
            u = np.unique(self.seq.flatten())
            lookup = {}
            for i in range(u.shape[0]):
                lookup[u[i]] = i 
            for i in range(self.seq.shape[0]):
                for j in range(self.seq.shape[1]):
                    array[i,j] = lookup[self.seq[i,j]]
                    
            labels = self.elec.label.values.tolist()
            l1 = labels[0]
            template = '{:0>2d}{:0>3d}' 
            if len(l1.split()) == 2:
                for i, label in enumerate(labels):
                    string = int(label.split()[0])
                    number = int(label.split()[1])
                    labels[i] = int(template.format(string, number))
            else:
                labels = np.asarray(labels, dtype=int)
                    
            sortidx = np.argsort(labels)
            elec = elec[sortidx]
            
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
        cmax = np.max(self.seq)
        ndig = len(str(cmax))
        nmeas = array.shape[0]
        cpairidn = [0]*nmeas # nested pair ids 
        cpairidu = [0]*nmeas # unested pair ids 
        idtemplate= '{:0>%id}{:0>%id}'%(ndig, ndig)
        for i in range(nmeas):
            celec = [self.seq[i,0], self.seq[i,1]]
            if aggressive: 
                celec = sorted(celec)
            uid = int(idtemplate.format(celec[0], celec[1]))
            if inested[i]: # dont attempt to string nested measurements together 
                cpairidn[i] = uid 
            else: 
                cpairidu[i] = uid 
            
        # note, want to sort by current electrodes so that the 
        # potential electrodes trial the current #1 electrodes for unnested measurements 
        uniquepairsu = np.unique(cpairidu)
        uniquepairsn = np.unique(cpairidn)
        
        lines = []
        for uid in uniquepairsu: 
            if uid == 0:
                continue 
            index = uid == cpairidu # index all instances of unique c1 and c2 being used for nested measurements 
            subset = self.seq[index] # subset of sequence 
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
            subset = self.seq[index] # subset of sequence 
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
                
        self.lines = lines 
        
        # make reciprocal lines 
        if self.seqR.shape[0] > 0:
            g = self.copy()
            g.seq = self.seqR.copy()
            g.seqR = np.array([])
            g.lines = [] 
            g.linesR = [] 
            self.linesR = g.multichannelise(maxcha) 
            
        return lines 
    
    def condition(self, plot=False, dump=print):
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
        dump('Conditioning sequence to avoid IP effects')
        # condition sequence to avoid using current electrodes
        # immiediately after measurement
        nswap = 0 # total number of swaps for all values of n, i, j etc... 
        ncmdline = len(self.lines)
        condLines = [self.lines[i] for i in range(ncmdline)]
        minsep= int(0.1*ncmdline)
        if minsep > 32:
            minsep = 32 
        condcmdsep = checkCmdSep(condLines , minsep)
        loop = 0 
        iswap = 9e9 
        iswapcache = -1 
        while min(condcmdsep) < minsep:
            loop += 1 
            dump('Resorting iteration %i...'%loop)
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
                        
            dump('Total swaps = %i'%nswap) 
            dump('Iteration swaps = %i'%iswap) 
            condcmdsep = checkCmdSep(condLines, minsep)
    
            if iswapcache <= iswap:# and loop>10: 
                dump('No improvement in number of swaps over last 2 iterations, treating as sorted...')
                break 
            
        if plot: 
            cmdsep = checkCmdSep(self.lines)
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
        
        self.lines = condLines 
        
        if len(self.linesR) > 0:
            g = self.copy()
            g.seq = self.seqR.copy()
            g.lines = [line for line in self.linesR]
            g.seqR = np.array([])
            g.linesR = [] 
            self.linesR = g.condition() 
            
        return condLines
    
    ## write to file 
    def write2csv(self, fname, singlechannel = False):
        """
        Write commands to generic csv type file 

        Parameters
        ----------
        fname : str
            Name / path of output file 

        """
        
        #check for csv ext
        if not fname.lower().endswith('.csv'):
            fname = fname + '.csv'
        
        fnamef = fname 
        fnamer = None 
        if self.seqR.shape[0] > 0 or len(self.linesR) > 0: 
            if fnamef.endswith('.csv'):
                fnamef = fname.replace('.csv','_F.csv')
                fnamer = fname.replace('.csv','_R.csv')
            else:
                fnamef = fname.replace('.CSV','_F.CSV')
                fnamer = fname.replace('.CSV','_R.CSV')
        fnames = [fnamef, fnamer]
        
        if singlechannel:
            for i, seq in enumerate([self.seq, self.seqR]):
                fname = fnames[i]
                if fname is None:
                    continue 
                if seq.shape[0] == 0:
                    continue 
                tmp = {
                    'a':seq[:,0], 
                    'b':seq[:,1],
                    'm':seq[:,2],
                    'n':seq[:,3],
                    }
                pd.DataFrame(tmp).to_csv(fname, index=False)
            return 
        
        lines2write = [self.lines, self.linesR]
        for i, lines in enumerate(lines2write):
            if len(lines) == 0:
                continue 
            fname = fnames[i]
            if fname is None:
                continue 
            fh = open(fname, 'w')
            nperline = len(lines[0])
            header = 'C1, C2, '
            while len(header.split(',')) <= nperline: 
                header += 'P%i, '%i
                i+=1 
            
            fh.write(header)
            
            for line in lines: 
                for n in line: 
                    fh.write('%i, '%n)
                fh.write('\n')
                
            fh.close() 

    def write2prime(self, fname):
        """
        Write to PRIME/RESIMGR command file. Note, limited to a maximum 
        of about 900 lines.

        Parameters
        ----------
        fname : str
            Name / path of output file 

        """
        fnamef = fname 
        fnamer = None 
        fext = '.' + fname.split('.')[-1]
        if len(self.linesR) > 0: 
            fnamef = fname.replace(fext,'_F' + fext)
            fnamer = fname.replace(fext,'_R' + fext)
        fnames = [fnamef, fnamer]
        lines2write = [self.lines, self.linesR]
        # nperline = len(self.lines[0])
        
        for i, lines in enumerate(lines2write):
            if len(lines) == 0:
                continue 
            if fname is None:
                continue 
            cf = 1 # count number of files written 
            fname = fnames[i].replace(fext, '_p{:0>2d}'.format(cf)+fext)
            cl = 1 # count number of lines written 
            fh = open(fname, 'w')
            fh.write('### insert headers here ###\n\n')
            fh.write('# test 0\n\n')
            for j, line in enumerate(lines): 
                fh.write('test %i '%cl)
                for n in line: 
                    fh.write('%i '%n)
                fh.write('\n')
                cl+=1
                if cl>=900:
                    fh.close() # close the file
                    cf += 1
                    fname = fnames[i].replace(fext, '_p{:0>2d}'.format(cf)+fext)
                    fh = open(fname, 'w') # make a new file 
                    fh.write('# test 0\n\n')
                    cl = 1 # reset count to 1 
                    
            if not fh.closed: 
                fh.close() 
                
    def write2sting(self, fname):
        """
        Write to Sting command file. 
        Parameters
        ----------
        fname : str
            Name / path of output file 

        """
        fnamef = fname 
        fnamer = None 
        fext = '.' + fname.split('.')[-1]
        if len(self.linesR) > 0: 
            fnamef = fname.replace(fext,'_F' + fext)
            fnamer = fname.replace(fext,'_R' + fext)
        fnames = [fnamef, fnamer]
        lines2write = [self.lines, self.linesR]
        # nperline = len(self.lines[0])
        
        lines2write = [self.lines, self.linesR]
        for i, lines in enumerate(lines2write):
            if len(lines) == 0:
                continue 
            fname = fnames[i]
            if fname is None:
                continue 
            fh = open(fname, 'w')
            # write out header information 
            fh.write(";Automatically created command file from ResIPy\n") 
            fh.write(":header\n") 
            typeid = "R" # R for resistance (no other types available according to sting docs)
            progid = "CUSTOM%iF"%len(lines)
            if i == 1: 
                progid = "CUSTOM%iR"%len(lines)
            fh.write("progID=%s\n"%progid) 
            fh.write("type=%s\n"%typeid)
            fh.write("arraytype=0\n") 
            fh.write("Binf=0\n")
            fh.write("Ninf=0\n")
            fh.write("MUX=1\n") 
            fh.write("\n")
            
            # write out electrode geometry 
            fh.write(":geometry\n")
            nelec = len(self.elec)
            for i in range(nelec):
                fh.write("{:d},{:f},{:f},{:f}\n".format(
                    i+1, 
                    self.elec.x[i], 
                    self.elec.y[i], 
                    self.elec.z[i]))
            fh.write("\n")    
            
            # write out command sequence 
            nperline = len(lines[0])
            fh.write(":commands\n") 
            header = ';A,B,'
            while len(header.split(',')) <= nperline: 
                header += 'P%i,'%i
                i+=1 
            header += 'channels\n'
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
            
    def write2syscal(self, fname):
        """
        Write to Syscal type command file. 
        Parameters
        ----------
        fname : str
            Name / path of output file 

        """
        fnamef = fname 
        fnamer = None 
        fext = '.' + fname.split('.')[-1]
        if len(self.linesR) > 0: 
            fnamef = fname.replace(fext,'_F' + fext)
            fnamer = fname.replace(fext,'_R' + fext)
        fnames = [fnamef, fnamer]
        lines2write = [self.lines, self.linesR]
        # nperline = len(self.lines[0])
        
        lines2write = [self.lines, self.linesR]
        for i, lines in enumerate(lines2write):
            if len(lines) == 0:
                continue 
            fname = fnames[i]
            if fname is None:
                continue 
            fh = open(fname, 'w')
            # write out electrode geometry 
            fh.write("#\tX\tY\tZ\n")
            nelec = len(self.elec)
            for i in range(nelec):
                fh.write("{:d}\t{:f}\t{:f}\t{:f}\n".format(
                    i+1, 
                    self.elec.x[i], 
                    self.elec.y[i], 
                    self.elec.z[i]))
            fh.write("\n")    
            
            # write out command sequence 
            fh.write("#\tA\tB\tM\tN\n")
            c = 1 
            for line in lines: 
                a = line[0]
                b = line[1]
                for j in range(2, len(line)-1): 
                    m = line[j] 
                    n = line[j+1]
                    if m == 0 or n == 0: 
                        continue 
                    fh.write("{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(c, a, b, m, n))
                    c += 1 
            fh.close() 
            
    ## export sequence 
    def exportSequence(self, fname, ftype = 'generic', integer = True, ignorestring=False, 
                       reciprocals = True,  split=True, multichannel = True,  
                       condition = True, maxcha = 8, dump=print): 
        """
        Export single channel sequence (ie for forward modelling) into something 
        usable by a resistivity instrument 
    
        Parameters
        ----------
        fname : str
            Export file path 
        ftype : str
            N by 4 array of electrode configurations 
        elec : pd.DataFrame
            Electrode data frame 
        integer : bool, optional
            Flag to convert sequence into integers before export. The default is True.
        ignorestring: bool, optional
            Flag to ignore line number in the case of multiple electrode strings 
            (or lines) being present in the electrode labels. Default is False. 
        reciprocals : bool, optional
            Flag to add reciprocal measurements. The default is True.
        split: bool, optional 
            If reciprocals added, split them across two files 
        multichannel : bool, optional
            Flag to convert measurements to a multichannel. The default is True.
        condition : bool, optional
            Flag to condition measurements for avoiding IP effects. The default is True.
        maxcha: int, optional 
            Maximum number of active channels of the resistivity instrument (normally
            8 for modern instruments). 
    
        """

        if ftype.lower() == 'asis':
            fext = '.csv'
            integer = False 
            reciprocals = False 
            split = False 
            multichannel = False 
            condition = False 
        elif ftype.lower() == 'prime': 
            fext = '.txt'
            maxcha = 7 
        elif ftype.lower() == 'syscal':
            fext = '.txt'
        elif ftype.lower() == 'sting':
            #todo: expand! 
            fext = '.cmd' 
        else:
            fext = '.csv'

        if not fname.endswith(fext):
            fname += fext
            
        if integer or multichannel: 
            # have to convert to integer to get multichannel code working currently 
            if ignorestring:
                print('Debug: Ignoring string in PRIME sequence')
                dump('Debug: Ignoring string in PRIME sequence')
            self.seq2int(ignorestring)
            
        if reciprocals: 
            dump('Generating reciprocals...')
            self.reciprocalise(split)
            
        if multichannel or condition: 
            self.multichannelise(maxcha)
        else: 
            dump("Warning: only .csv format available for singlechannel measurements!")
            dump("writing to file...")
            self.write2csv(fname, True)
            return 
        
        if condition:
            self.condition(dump=dump)
            
        dump('Writing to file...')
        if ftype.lower() == 'prime':
            self.write2prime(fname)
        elif ftype.lower() == 'sting':
            self.write2sting(fname)
        elif ftype.lower() == 'syscal':
            self.write2syscal(fname)
        else:
            self.write2csv(fname)
        return 
    
