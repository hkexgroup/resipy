#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:21:54 2018

@author: ResIPy's core developers
"""
import sys, os, platform, time 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.ticker import MaxNLocator
import pandas as pd
from scipy.stats import norm, linregress
from scipy.stats import gaussian_kde
from scipy.linalg import lstsq
from functools import reduce

from resipy.parsers import (syscalParser, protocolParserLME, resInvParser,
                     primeParserTab, protocolParser,
                     stingParser, ericParser, lippmannParser, aresParser,
                     srvParser, bertParser, dasParser, electraParser)
from resipy.DCA import DCA
from resipy.interpolation import geometricMedian 
from resipy.saveData import to_csv

# show the deprecation warnings
import warnings
warnings.simplefilter('default', category=DeprecationWarning)

try:#import pyvista if avaiable
    import pyvista as pv
    try:
        from pyvistaqt import BackgroundPlotter # newer version
    except:
        from pyvista import BackgroundPlotter # older version
    pyvista_installed = True
except ModuleNotFoundError:
    pyvista_installed = False
    # warnings.warn('pyvista not installed, 3D meshing viewing options will be limited')
    
try:
    from resipy.cext.fastRecip import fastReciprocals
    fastrecip_installed = True 
except: 
    fastrecip_installed = False 
    
#replacement for numpy polyfit function which works on open blas 
def polyfit(x,y,deg=1):
    """Replacement function for numpy polyfit that works consistently (avoids 
    SVD convergence error prsent on some W10 computers) Nb: is not as robust
    as numpy polyfit function. 
    
    Parameters
    ----------
    x : array like
        x values of data
    y : TYPE
        y values of data
    deg : TYPE, optional
        DESCRIPTION. The default is 1.

    Raises
    ------
    ValueError
        If deg < 0 .
    TypeError
        array data is not as expected for fitting a polynomail.

    Returns
    -------
    coef : np.array 
        coefficents of polynomail 

    """
    x = np.asarray(x,dtype=float)
    y = np.asarray(y,dtype=float)
    
    if any(np.isnan(x)) or any(np.isinf(x)) or any(np.isnan(y)) or any(np.isinf(y)):
        # fit to just real numbers 
        keep_idx = ~np.isnan(x) & ~np.isinf(x) & ~np.isnan(y) & ~np.isinf(y)
        x = x[keep_idx]
        y = y[keep_idx]
        

    # check arguments.
    if deg < 0:
        raise ValueError("expected deg >= 0")
    if x.ndim != 1:
        raise TypeError("expected 1D vector for x")
    if x.size == 0:
        raise TypeError("expected non-empty vector for x")
    if y.ndim < 1 or y.ndim > 2:
        raise TypeError("expected 1D or 2D array for y")
    if x.shape[0] != y.shape[0]:
        raise TypeError("expected x and y to have same length")

    # set up least squares equation for powers of x, to be solved in the form Ax = y
    A = np.ones((len(x),deg+1))
    exps = np.arange(deg,0,-1) # exponents of polynomail 
    
    #construct A
    if deg == 0: # just fit y = mx
        A[:,0] = x
    else:
        count = 0
        for e in exps:
            A[:,count] = x**e
            count+=1

    #now solve 
    coef, resids, rank, s = lstsq(A, y, lapack_driver = 'gelss')
    
    if all(coef) == 0 or any(np.isnan(coef)) == True:
        #then try again 
        coef, resids, rank, s = lstsq(A, y, lapack_driver = 'gelss')
        
    return coef 

# def fixSequence(sequence, nelec):
#     """
#     Make sequence consecutive. 

#     Parameters
#     ----------
#     sequence : nd array 
#         N by 4 array which is the measurement sequence. 

#     Returns
#     -------
#     None.

#     """
#     uid = np.unique(sequence.flatten())
#     cid = np.arange(nelec)+1
#     # print('Number of elec in electrode df: %i, Number of elec in sequence: %i'%(len(uid),nelec))
#     newseq = sequence.copy()
#     for i in range(len(uid)):
#         if uid[i] != cid[i]:
#             replaceidx = sequence == uid[i]
#             newseq[replaceidx] = cid[i] 
#     return newseq 

def appendString(elec, df, string):
    """
    Append strings to the data frame and electrode labels in the case of 3D 
    surveys. 

    Parameters
    ----------
    elec : pandas.DataFrame
        Electrode data frame with columns, x,y,z
    df : pandas.DataFrame
        Resistance measurement dataframe with columns a,b,m,n
    """
    if 'label' not in elec.columns:
        elec['label'] = (np.arange(len(elec))+1).astype(str)
        
    def getLabelNumber(label):
        if len(label.split()) == 2:
            number = int(label.split()[-1])
        else:
            number = int(label)
        return number 
        
    for i,label in enumerate(elec.label):
        number = getLabelNumber(label)
        new_label = '{:d} {:d}'.format(string, number) 
        elec.loc[i,'label'] = new_label 
    
    
    ndata = len(df)
    for i in range(ndata):
        for e in ['a','b','m','n']:
            label = df[e].values[i]
            number = getLabelNumber(label)
            new_label = '{:d} {:d}'.format(string, number) 
            df.loc[i,e] = new_label 
            
    

class Survey(object):
    """Class that handles geophysical data and some basic functions. One 
    instance is created for each survey.
    
    Parameters
    ----------
    fname : str
        Name of the file where the data are.
    ftype : str
        Type of the data file. This setting is not read if a `parser` is
        not `None`.
    df : pandas.DataFrame, optional
        Pandas dataframe with 'a','b','m','n' columns as string and at least a column
        'resist' as float. Can be provided alternatively to fname.
    elec : pandas.DataFrame, optional
        Pandas dataframe with 'x','y','z','buried','remote' columns.
        'x','y','z' as float, 'buried' and 'remote' as bool. Should be provided
        alternatively to fname.
    name : str
        A personal name for the survey.
    spacing : float, optional
        This will be passed to the parser function to determine the
        electrode positions.
    parser : function, optional
        If provided, it should return tuple containing `elec` a 3 columns
        array containing electrodes position and `data` a pandas.DataFrame
        with the `a,b,m,n,resist` columns at least and `ip` if present.
    keepAll: bool, optional
        If `True` will keep all the measurements even the ones without
        reciprocal. Note that if none of the quadrupoles have reciprocal
        they will all be kept anyway.
    compRecip: bool, optional 
        Compute reciprocal errors, default is True. 
    string: int, optional
        Force the addition of string number in electrode labels and abmn 
        values. Useful when creating a survey that consists of multiple 
        files. 
    """
    def __init__(self, fname=None, ftype='', df=None, elec=None, name='',
                 spacing=None, parser=None, keepAll=True, debug=True,
                 compRecip=True, string=0):
        
        # set default attributes
        self.iBorehole = False # True is it's a borehole
        self.protocolIPFlag = False # True if the file to be parsed is a Protocol file with IP
        self.kFactor = 1 # factor to convert between chargeability (mV/mV) and phase (mrad)
        self.errorModel = None # function instantiated after fitting an error model with reciprocal errors
        self.phaseErrorModel = None # idem for IP
        self.iselect = None # use in filterManual()
        self.eselect = None # idem
        self.debug = debug # plotting all information message by default
        self.phiCbarmin = 0
        self.phiCbarMax = 25
        self.filt_typ = 'Raw'
        self.cbar = True
        self.filterDataIP = pd.DataFrame()

        # check arguments
        if name == '':
            if fname is not None:
                # tmp = fname.split(os.sep)[-1] #gets file name (without path)
                name = os.path.basename(os.path.splitext(fname)[0])
            else:
                name = 'Survey'
        self.name = name
        self.source = fname #save original file name 
        
        if spacing is not None:
            warnings.warn('The spacing argument is deprecated and will be removed in the next version.',
                      DeprecationWarning)
        
        # parsing data to form main dataframe and electrode array
        if fname is not None:
            avail_ftypes = ['Syscal','ProtocolDC','ResInv', 'BGS Prime', 'RESIMGR', 
                            'ProtocolIP','Sting', 'ABEM-Lund', 'Lippmann', 
                            'ARES', 'E4D', 'BERT', 'Electra', 'DAS-1']# add parser types here! 
            if parser is not None:
                elec, data = parser(fname)
            else:
                if ftype == 'Syscal':
                    elec, data = syscalParser(fname)
                    self.kFactor = 1.2
                elif ftype =='ProtocolDC':
                    elec, data = protocolParser(fname, ip=False)
                elif ftype == 'ResInv':
                    elec, data = resInvParser(fname)
                elif ftype == 'BGS Prime':
                    elec, data = primeParserTab(fname)
                elif ftype == 'RESIMGR': # same format as PRIME 
                    elec, data = primeParserTab(fname)
                elif ftype == 'ProtocolIP':
                    elec, data = protocolParser(fname, ip=True)
                    self.protocolIPFlag = True
                elif ftype == 'forwardProtocolDC':
                    elec, data = protocolParser(fname, fwd=True)
                elif ftype == 'forwardProtocolIP':
                    elec, data = protocolParser(fname, ip=True, fwd=True)
                    self.protocolIPFlag = True
                elif ftype == 'Sting':
                    elec, data = stingParser(fname)
                elif ftype == 'ABEM-Lund':
                    elec, data = ericParser(fname)
                elif ftype == 'Lippmann':
                    elec, data = lippmannParser(fname)
                elif ftype == 'ARES':
                    elec ,data = aresParser(fname)
                elif ftype == 'E4D':
                    elec ,data = srvParser(fname)
                elif ftype == 'BERT':
                    elec, data = bertParser(fname)
                elif ftype == 'Electra':
                    elec, data = electraParser(fname)
                elif ftype == 'DAS-1':
                    elec, data = dasParser(fname)
                else:
                    print("Unrecognised ftype, available types are :")
                    for ftype in avail_ftypes: 
                        print('\t%s'%ftype)
                    raise Exception('Sorry this file type is not implemented yet')

            # assign dataframe and check the types of a,b,m,n (all labels must be string)
            self.df = data.astype({'a':str, 'b':str, 'm':str, 'n':str})
            self.df.reset_index(inplace=True, drop=True) 

            # add error measured to the error columns (so they can be used if no error model are fitted)
            if 'magErr' in self.df.columns:
                self.df['resError'] = self.df['magErr'].copy()
            else:
                self.df['resError'] = np.nan
            if 'phiErr' in self.df.columns:
                self.df['phaseError'] = self.df['phiErr'].copy()
            else:
                self.df['phaseError'] = np.nan

            # set electrode dataframe
            if type(elec) == np.ndarray: # TODO for temporary stability
                self.elec = pd.DataFrame(elec.astype(float), columns=['x','y','z'])
                self.elec['remote'] = False
                self.elec['buried'] = False
                self.elec['label'] = (1 + np.arange(self.elec.shape[0])).astype(str)
            else:
                self.elec = elec
               
            # assign string numbers to electrode labels 
            if string > 0: 
                appendString(self.elec, self.df, string) 
                
            # set sequence according to if electrode labels present or not 
            self.setSeqIds() 

            # convert apparent resistivity to resistance and vice versa
            self.computeK()
            if 'resist' in self.df.columns:
                self.df['app'] = self.df['K']*self.df['resist']
            elif 'app' in self.df.columns:
                self.df['resist'] = self.df['app']/self.df['K']

        # we provide df and elec dataframes
        elif df is not None and elec is not None:
            self.df = df # assuming dataframe is already well formatted

            # check elec dataframe
            if 'remote' not in elec.columns:
                elec['remote'] = False
            if 'buried' not in elec.columns:
                elec['buried'] = False
            if 'label' not in elec.columns:
                elec['label'] = np.arange(elec.shape[0]).astype(str)
            self.elec = elec
            
            # house keeping 
            self.setSeqIds()
            self.computeK()
        
        else:
            raise ValueError('No fname supplied, no df and elec supplied. Returned.')
            return

        # apply basic filtering
        self.filterDefault(False) # dont compute reciprocals yet 
        if compRecip:
            self.computeReciprocal() # compute reciprocals
        
        # create backup copies
        self.dfReset = self.df.copy()
        self.dfPhaseReset = self.df.copy()
        self.dfOrigin = self.df.copy()
        self.isequenceReset = self.isequence.copy() # for retetting filters in the UI
       
     
#     @classmethod
#     def fromDataframe(cls, df, dfelec):
#         """Create a survey class from pandas dataframe.
        
#         Parameters
#         ----------
#         df : pandas.DataFrame
#             Pandas dataframe.
#         dfelec : pandas.DataFrame
#             Pandas dataframe for electrodes.
            
#         Returns
#         -------
#         survey : Survey
#             An instance of the Survey class.
#         """
        
#         def resource_path(relative_path): # a better method to get relative path of files
#             if hasattr(sys, '_MEIPASS'):
#                 return os.path.join(sys._MEIPASS, relative_path)
#             return os.path.join(os.path.abspath("."), relative_path)
#         path = resource_path('examples/dc-2d/protocol.dat')        
        
#         # surrogateFile = '/src/examples/dc-2d/protocol.dat' # tricky
#         # path = os.path.realpath(__file__).replace(
#         #     '\\src\resipy\\Survey.py', surrogateFile).replace(
#         #     '/src/resipy/Survey.py', surrogateFile)

#         survey = cls(fname=path, ftype='ProtocolDC') #surrogate survey
#         survey.df = df
#         if 'remote' not in dfelec.columns:
#             dfelec['remote'] = False
#         if 'buried' not in dfelec.columns:
#             dfelec['buried'] = False
#         if 'label' not in dfelec.columns:
#             dfelec['label'] = np.arange(dfelec.shape[0]).astype(str)
#         survey.elec = dfelec
#         survey.computeReciprocal()
#         survey.filterDefault()
#         survey.dfReset = survey.df.copy()
#         survey.dfPhaseReset = survey.df.copy()
        
#         return survey
    
    
    def __str__(self):
        out = "Survey class with {:d} measurements and {:d} electrodes".format(
            self.df.shape[0], self.elec.shape[0])
        return out
    
    def hasElecString(self):
        """Determine if a electrode strings are present in the electrode labels 
        in self.elec

        Returns
        -------
        bool
            True if strings present in electrode labels 
        """
        df = self.elec
        if 'label' not in df.keys():
            return False
        else:
            label = df['label']
            for l in label:
                if len(l.split()) == 1:
                    return False
        return True
    
    def hasDataString(self):
        """
        Check if string numbers are present in the data. 

        Returns
        -------
        bool
            True if strings present in electrode labels 

        """
        if len(self.df) == 0:
            return False 
        if len(self.df['a'].values[0].split()) <= 1:
            return False 
        return True 
        
    
    def setSeqIds(self): 
        """Convert electrode labels to indexable integers, sets the class wide 
        parameter 'isequence'
        """
        ndata = self.df.shape[0]

        if ndata == 0:
            self.isequence = None 
            return 
        
        if self.hasDataString() != self.hasElecString(): 
            return 

        lookup = {}
        for i,label in enumerate(self.elec.label): 
            lookup[label] = i + 1 

        self.isequence = np.zeros((ndata, 4), dtype=int)
        for i in range(ndata):
            for j,char in enumerate(['a','b','m','n']):
                self.isequence[i,j] = lookup[self.df[char].values[i]]
        

    def convertLocalGrid(self):
        """Converts UTM grid to local grid for mesh stability"""
        electemp = self.elec.copy()
        electemp['order'] = np.arange(1, len(electemp)+1)
        electempXsort = electemp.sort_values(by='x')
        x_local = electempXsort['x'].values - electempXsort['x'].values[0]
        electempXsort['x_orig'] = electempXsort['x'].values # keep the original values in case we want to go back
        electempXsort['x'] = x_local
        electempYsort = electempXsort.sort_values(by='y')
        y_local = electempYsort['y'].values - electempYsort['y'].values[0]
        electempYsort['y_orig'] = electempYsort['y'].values # keep the original values in case we want to go back
        electempYsort['y'] = y_local
        self.elec = electempYsort.sort_values(by='order').drop(['order'], axis=1)
    
    
    def checkTxSign(self, inplace=True):
        """Check the sign of the measurement apparent resistivities. For example
        this is a necessary step for surveys using the syscal instrument 
        as the polarity of measurements is not recorded. 
        
        Parameters
        ----------
        inplace: bool 
            If True, Change the sign of transfer if they are negative. Defualt 
            is True. 
        """
        resist = self.df['resist'].values.copy()
        self.computeK()
        K = self.df['K'].values
        ie = ((K < 0) & (resist > 0)) | ((K > 0) & (resist < 0))
        
        if not inplace: 
            return ie # drop out of function now a return indexing array 
        self.df.loc[ie, 'resist'] = resist[ie]*-1
        if 'recipMean' in self.df.columns: # in case this method is called after importation
            recipMean = self.df['recipMean'].values.copy()
            self.df.loc[ie, 'recipMean'] = recipMean[ie]*-1
        if 'recipMean0' in self.df.columns: # time-lapse: in case this method is called after importation
            recipMean0 = self.df['recipMean0'].values.copy()
            self.df.loc[ie, 'recipMean0'] = recipMean0[ie]*-1

        print('WARNING: change sign of ', np.sum(ie), ' Tx resistance.')
        return ie 
        

    def filterDefault(self, recompute_recip=True):
        """Remove NaN, Inf and duplicates values in the data frame.
        """
        if self.debug:
            def dump(x):
                print(x)
        else:
            def dump(x):
                pass
            
        # remove Inf and NaN
        resist = self.df['resist'].values
        iout = np.isnan(resist) | np.isinf(resist)
        if np.sum(iout) > 0 and self.debug:
            dump('Survey.filterDefault: Number of Inf or NaN : {:d}\n'.format(np.sum(iout)))
        self.filterData(~iout)

        # average duplicates
        # let's first sort quadrupoles as ABMN to check 
        # for duplicates like BAMN, ABNM, BAMN
        shapeBefore = self.df.shape[0]
        ab = self.isequence[:, :2] - 1
        mn = self.isequence[:, 2:] - 1
        ab_sorted = np.sort(ab, axis=1)
        mn_sorted = np.sort(mn, axis=1)
        iab = (ab == ab_sorted).all(axis=1)
        imn = (mn == mn_sorted).all(axis=1)
        # changing sign in case we have BAMN or ABNM
        i2change = np.ones(len(iab))
        i2change[(iab & ~imn) | (~iab & imn)] = -1
        self.df.loc[:, ['ia', 'ib']] = ab_sorted
        self.df.loc[:, ['im', 'in']] = mn_sorted
        self.df.loc[:, 'i2change'] = i2change
        self.df.loc[:, 'iresist'] = self.df['resist']*i2change
        self.df['initial_index'] = np.arange(self.df.shape[0])
        
        # averaging resist, taking first of all other columns
        cols = [col for col in self.df.columns if col not in ['ia', 'ib', 'im', 'in']]
        aggdic = dict(zip(cols, ['first']*len(cols)))
        aggdic['iresist'] = 'mean'
        self.df = self.df.groupby(['ia', 'ib', 'im', 'in']).agg(aggdic).reset_index(drop=True)
        self.df = self.df.sort_values('initial_index').reset_index(drop=True)
        self.setSeqIds() # the groupby disturb the order of the sequence
        self.df['resist'] = self.df['iresist'] * self.df['i2change']
        self.df = self.df.drop(['initial_index', 'iresist'], axis=1)
        ndup = shapeBefore - self.df.shape[0]
        if ndup > 0 and self.debug:
            dump('Survey.filterDefault: {:d} duplicates averaged.\n'.format(ndup))
        
        # # remove duplicates
        # shapeBefore = self.df.shape[0]
        # self.df = self.df.drop_duplicates(subset=['a','b','m','n'], keep = 'first')
        # ndup = shapeBefore - self.df.shape[0]
        # if ndup > 0 and self.debug:
        #     dump('Survey.filterDefault: {:d} duplicates removed.\n'.format(ndup))
        
        # # remove quadrupoles were A or B are also potential electrodes
        # ie1 = self.df['a'].values == self.df['m'].values
        # ie2 = self.df['a'].values == self.df['n'].values
        # ie3 = self.df['b'].values == self.df['m'].values
        # ie4 = self.df['b'].values == self.df['n'].values
        # ie = ie1 | ie2 | ie3 | ie4
        # if np.sum(ie) > 0 and self.debug:
        #     dump('Survey.filterDefault: {:d} measurements with A or B == M or N\n'.format(np.sum(ie)))
        # self.filterData(~ie)
        
        # # remove repeated readings like ABMN, BAMN, ABNM
        # df = self.df[['a', 'b', 'm', 'n']].copy()
        # df[['a', 'b']] = np.sort(df[['a', 'b']].values, axis=1)
        # df[['m', 'n']] = np.sort(df[['m', 'n']].values, axis=1)
        # ie = df.duplicated(subset=['a', 'b', 'm', 'n'])
        # if np.sum(ie) > 0 and self.debug:
        #     dump('Survey.filterDefault: {:d} duplicates ABMN, BAMN or ABNM removed\n'.format(np.sum(ie)))
        # self.filterData(~ie)
        
        # create a backup of the clean dataframe
        self.dfReset = self.df.copy()
        self.dfPhaseReset = self.df.copy()
        self.isequenceReset = self.isequence.copy()
        
        
    def addData(self, fname, ftype='Syscal', parser=None, string=0):
        """Add data to the actual survey (for instance the reciprocal if they
        are not in the same file).
        """
         
        if parser is not None:
            elec, data = parser(fname)
        else:
            if ftype == 'Syscal':
                elec, data = syscalParser(fname)
                self.kFactor = 1.2
            elif ftype =='ProtocolDC':
                elec, data = protocolParser(fname, ip=False)
            elif ftype == 'ResInv':
                elec, data = resInvParser(fname)
            elif ftype == 'BGS Prime':
                elec, data = primeParserTab(fname)
            elif ftype == 'ProtocolIP':
                elec, data = protocolParser(fname, ip=True)
                self.protocolIPFlag = True
            elif ftype == 'Sting':
                elec, data = stingParser(fname)
            elif ftype == 'ABEM-Lund':
                elec, data = ericParser(fname)
            elif ftype == 'Lippmann':
                elec, data = lippmannParser(fname)
            elif ftype == 'ARES':
                elec ,data = aresParser(fname)
            elif ftype == 'E4D':
                elec ,data = srvParser(fname)
            elif ftype == 'BERT':
                elec, data = bertParser(fname)
            elif ftype == 'Electra':
                elec, data = electraParser(fname)
            elif ftype == 'DAS-1':
                elec, data = dasParser(fname)
            else:
                raise Exception('Sorry this file type is not implemented yet or not recognised')
        
        data = data.astype({'a':str, 'b':str, 'm':str, 'n':str})
        
        if type(elec) == np.ndarray: # 
            _elec = pd.DataFrame(elec.astype(float), columns=['x','y','z'])
            _elec['remote'] = False
            _elec['buried'] = False
            _elec['label'] = (1 + np.arange(elec.shape[0])).astype(str)
        else:
            _elec = elec.copy()
            if 'remote' not in _elec.columns:
                _elec['remote'] = False
            if 'buried' not in _elec.columns:
                _elec['buried'] = False
            if 'remote' not in _elec.columns:
                _elec['label'] = (1 + np.arange(elec.shape[0])).astype(str)
                
        
        # check if new electrodes present 
        if string==0:
            for i,label in enumerate(_elec.label):
                if label not in self.elec.label.values: 
                    line = _elec[_elec.index==i]
                    self.elec = pd.concat([self.elec,line]).reset_index(drop=True)
        elif string > 0: 
            data.reset_index(drop=True, inplace=True)
            appendString(_elec, data, string)
            self.elec = pd.concat([self.elec,_elec]).reset_index(drop = True) 
        
        self.df = pd.concat([self.df, data]).reset_index(drop = True) # for being similar to import one df with reciprocals (as new df will have its own indices!)
        self.setSeqIds() # need to recompute sequence 
            
        self.dfOrigin = self.df.copy()
        self.filterDefault()
        self.computeReciprocal()
        
        # backing up data
        self.dfReset = self.df.copy()
        self.dfPhaseReset = self.df.copy()
        self.dfOrigin = self.df.copy()
        self.isequenceReset = self.isequence.copy()

    
    def filterData(self, i2keep):
        """Filter out the data not retained in `i2keep`.
        
        Parameters
        ----------
        i2keep : ndarray of bool
            Index where all measurement to be retained are `True` and the
            others `False`.
        """
        if len(i2keep) != self.df.shape[0]:
            if 'ip' not in self.df.columns:
                raise ValueError('The length of index to be kept (' + str(len(i2keep)) + ')\n'
                                 'does not match the length of the data (' + str(self.df.shape[0]) +').')
            else:
                raise ValueError('The length of index to be kept (' + str(len(i2keep)) + ') '
                                 'does not match the length of the data (' + str(self.df.shape[0]) +').\n'
                                 'Reciprocal Filtering cannot be done after Phase Filtering.\n'
                                 'Reset the filters and redo the filterings, first reciprocity then phase.')
            return
        else:
            if 'irecip' in self.df.columns:
                # get a list of measurement that would be affected by the removal
                recip2reset = self.df[~i2keep]['irecip'].values*-1
            self.df = self.df[i2keep]
            if 'index' in self.df.columns: 
                self.df.reset_index(drop=True, inplace=True)
            else:
                self.df.reset_index(inplace=True)
            
            self.dfPhaseReset = self.df.copy()
            if 'irecip' in self.df.columns:
                ie = np.in1d(self.df['irecip'].values, recip2reset)
                self.df.loc[ie, 'irecip'] = 0 # as their reciprocal is deleted, we set it to 0
                self.df.loc[ie, 'recipError'] = np.nan # they don't contribute to the error model anymore
                self.df.loc[ie, 'recipMean'] = self.df.loc[ie, 'resist'].values
            if self.debug is True:
                print('filterData:', np.sum(~i2keep), '/', len(i2keep), 'quadrupoles removed.')
                
            self.setSeqIds() # not super efficient, but works 
            return np.sum(~i2keep)

    
    
    def filterUnpaired(self):
        """Remove quadrupoles that don't have a reciprocals. This might
        remove dummy measurements added for sequence optimization.
        """
        i2keep = self.df['irecip'] != 0
        if self.debug: 
            print('removeUnpaired:', end='')
        self.filterData(i2keep)
        return np.sum(~i2keep)
    
    def computeReciprocal(self, alg='Bisection Search'):
        """
        Compute Reciprocals and store them in self.df (the dataframe)

        Parameters
        ----------
        alg : str, optional
            Algorithm used to compute reciprocals. Choose between 
            'Bisection Search', 'Pandas Merge' or 'Array Expansion'.
            The default is 'Bisection Search', other string are casted to
            'Array Expansion'.
            
        Notes
        -----
        If the relevant cython code cannot be found then the function falls 
        back to using the pandas merge approach. 

        """
        possible_algs = ['Pandas Merge','Array Expansion','Bisection Search']
        infos = ['Uses Pandas merge function to efficiently search for reciprocals',
                 'Creates expansive arrays to find reciprocals using numpy (can be RAM intense for large Surveys)',
                 'Uses an efficient bisection search algorithm to find reciprocals']
        
        if alg not in possible_algs:
            print('Selected Reciprocal error calculation not recognised, possible algorithms include:')
            for i, palg in enumerate(possible_algs):
                print(palg,' : ',infos[i])
            raise Exception('Reciprocal error algorithm not recognised')
            
        ## need to put something in here in case the cython module can't be found 
        if not fastrecip_installed: 
            alg = 'Pandas Merge' 
            
        if alg == 'Bisection Search':
            self.computeReciprocalC()
        elif alg == 'Pandas Merge':
            self.computeReciprocalP()
        else: 
            self.computeReciprocalN()
    
    def computeReciprocalP(self): # pandas way
        """Compute reciprocal measurements using `pandas.merge()`.
        
        Notes
        -----
        This algorithm create two dataframes with ABMN and MNAB then use the pd.merge()
        function to match normal and their reciprocal. The resulting merged dataframe
        contains indexes of both ABMN and MNAB dataframes and is used to populate the irecip
        column.
        """
        resist = self.df['resist'].values
        phase = -self.kFactor*self.df['ip'].values #converting chargeability to phase shift
        ndata = self.df.shape[0]
        array = self.isequence - 1 
        
        R = np.copy(resist)
        M = np.copy(phase)
        Ri = np.zeros(ndata, dtype=int)
        reciprocalErr = np.zeros(ndata)*np.nan
        reciprocalErrRel = np.zeros(ndata)*np.nan
        reciprocalMean = np.zeros(ndata)*np.nan
        reci_IP_err = np.zeros(ndata)*np.nan
        
        # sort quadrupoles and creates reciprocal array
        sortedArray = np.c_[np.sort(array[:,:2]), np.sort(array[:,2:])]
        
        # build dataframe of normal and reciprocal and merge them
        df1 = pd.DataFrame(sortedArray, columns=['a','b','m','n'])
        df1['index1'] = np.arange(df1.shape[0])
        df2 = pd.DataFrame(sortedArray, columns=['m','n','a','b'])
        df2['index2'] = np.arange(df2.shape[0])
        dfm = pd.merge(df1, df2, on=['a','b','m','n'], how='outer')
        dfm = dfm.dropna()
        
        ###### WIP vvvvvvv
        # df3 = pd.DataFrame(sortedArray, columns=['n','m','a','b'])
        # df3['index3'] = np.arange(df3.shape[0])
        # df4 = pd.DataFrame(sortedArray, columns=['m','n','b','a'])
        # df4['index4'] = np.arange(df4.shape[0])
        # dfs = [df1, df2, df3, df4]
        # dfm = reduce(lambda  left,right: pd.merge(left,right,on=['a','b','m','n'], how='outer'), dfs)
        
        # dfm = dfm.dropna(subset=['index1', 'index2', 'index3', 'index4'], how='all') # drop quad without any reciprocals
        # def countna(x):
        #     arr = np.array(x)
        #     c = np.count_nonzero(np.isnan(arr))
        #     return True if c == 2 else False # there should be only one reciprocal pair for a normal quadrupole 
        # dfm['recippair'] = dfm[['index1', 'index2', 'index3', 'index4']].apply(countna, axis=1)
        # dfm = dfm[dfm['recippair']]
        # dfm = dfm.dropna(axis=1, how='all') #dropping columns with all nans (i.e., no pair indecies) # TODO: fix this for normal only data so it doesn't drop all cols.
        # indexcols = [col for col in dfm.columns if 'index' in col]
        # indexArray = np.sort(dfm[indexcols].values.astype(int), axis=1)
        ###### WIP ^^^^^^
        
        # sort and keep only half
        indexArray = np.sort(dfm[['index1', 'index2']].values.astype(int), axis=1)
        indexArrayUnique = np.unique(indexArray, axis=0)
        inormal = indexArrayUnique[:,0] # I guess it doesn't matter which column is normal and which column is reciprocal; we're using np.abs below
        irecip = indexArrayUnique[:,1]
        
        val = np.arange(ndata) + 1
        Ri[inormal] = val[inormal]
        Ri[irecip] = -val[inormal]

        # establishing information on which swapping of electrode within dipole
        # is expected to cause a change in sign between normal and reciprocal quad.
        ab_diff = (array[:, 0] - sortedArray[:, 0]) != 0  # if != 0, then they were sorted
        mn_diff = (array[:, 2] - sortedArray[:, 2]) != 0 
        idiff = (ab_diff & mn_diff) | (~ab_diff & ~mn_diff)
        swapped = np.ones(array.shape[0], dtype=int)
        swapped[inormal] -= (idiff[inormal] != idiff[irecip])*2
        
        # compute the reciprocal error (using the magnitude of the measurements)
        reciprocalErr[inormal] = R[irecip] - R[inormal]*swapped[inormal]
        reciprocalErr[irecip] = R[irecip] - R[inormal]*swapped[inormal]
        reci_IP_err[inormal] = M[irecip] - M[inormal]
        reci_IP_err[irecip] = M[inormal] - M[irecip]
        
        # compute reciprocal mean with all valid values
        ok1 = ~(np.isnan(R[inormal]) | np.isinf(R[inormal]))
        ok2 = ~(np.isnan(R[irecip]) | np.isinf(R[irecip]))
        
        ie = ok1 & ok2 # both normal and recip are valid
        reciprocalMean[inormal[ie]] = np.mean(np.c_[np.abs(R[inormal]),np.abs(R[irecip])], axis=1)
        reciprocalMean[irecip[ie]] = np.mean(np.c_[np.abs(R[inormal]),np.abs(R[irecip])], axis=1)
        
        ie = ok1 & ~ok2 # only use normal
        reciprocalMean[inormal[ie]] = np.abs(R[inormal[ie]])
        
        ie =  ~ok1 & ok2 # only use reciproal
        reciprocalMean[inormal[ie]] = np.abs(R[irecip[ie]])
        
        reciprocalErrRel = reciprocalErr / reciprocalMean
        
        reciprocalMean = np.sign(resist)*reciprocalMean # add sign
        
        ibad = np.array([np.abs(a) > 0.2 if ~np.isnan(a) else False for a in reciprocalErrRel])
        if self.debug:
            print('{:d}/{:d} reciprocal measurements found.'.format(np.sum(Ri != 0), len(Ri)))
            if np.sum(Ri != 0) > 0: # no need to display that if there is no reciprocal
                print('{:d} measurements error > 20 %'.format(np.sum(ibad)))        
                
        self.df['irecip'] = Ri
        self.df['reciprocalErrRel'] = reciprocalErrRel
        self.df['recipError'] = reciprocalErr
        self.df['recipMean'] = reciprocalMean
        self.df['reci_IP_err'] = reci_IP_err
        # in order to compute error model based on a few reciprocal measurements
        # we fill 'recipMean' column with simple resist measurements for lonely
        # quadrupoles (which do not have reciprocals)
        inotRecip = Ri == 0
        self.df.loc[inotRecip, 'recipMean'] = self.df.loc[inotRecip, 'resist']
        
        return Ri
    
    
    def computeReciprocalN(self): # fast vectorize version
        """Compute reciprocal measurements using numpy array expansion. 
        
        Notes
        -----
        The method first sorts the dipole AB and MN. Then creates a reciprocal
        quadrupole matrix. These matrices are then used with
        numpy.equal to produce a 2D matrix from which reciprocal are extracted.
        """
        resist = self.df['resist'].values
        phase = -self.kFactor*self.df['ip'].values #converting chargeability to phase shift
        ndata = self.df.shape[0]
        array = self.isequence - 1 
        
        R = np.copy(resist)
        M = np.copy(phase)
        ndata = len(R)
        Ri = np.zeros(ndata)
        reciprocalErr = np.zeros(ndata)*np.nan
        reciprocalErrRel = np.zeros(ndata)*np.nan
        reciprocalMean = np.zeros(ndata)*np.nan
        reci_IP_err = np.zeros(ndata)*np.nan
        
        # sort quadrupoles and creates reciprocal array
        sortedArray = np.c_[np.sort(array[:,:2]), np.sort(array[:,2:])]
        recipArray = np.c_[sortedArray[:,2:], sortedArray[:,:2]]
        
        # build matching matrix
        imatch = np.equal(sortedArray[:,None,:], recipArray[None,:,:]).all(2)
        imatch = np.triu(imatch)
        
        inotfound = np.sum(imatch, axis=1) != 1 # all not found (or duplicates)
        imatch[inotfound, :] = False # in case one quad has multiple recip, this discards it
        inormal, irecip = np.where(imatch) # quad with only one reciprocal
        
        val = np.arange(ndata) + 1
        Ri[inormal] = val[inormal]
        Ri[irecip] = -val[inormal]
        
        # check expected polarity of reciprocal measurement?? code might look something like the below 
        # flipsign = (array[inormal,0] != array[irecip,2]) & (array[inormal,1] != array[irecip,3])
        # R[flipsign] *= -1 
        
        # establishing information on which swapping of electrode within dipole
        # is expected to cause a change in sign between normal and reciprocal quad.
        ab_diff = (array[:, 0] - sortedArray[:, 0]) != 0  # if != 0, then they were sorted
        mn_diff = (array[:, 2] - sortedArray[:, 2]) != 0 
        idiff = (ab_diff & mn_diff) | (~ab_diff & ~mn_diff)
        swapped = np.ones(array.shape[0], dtype=int)
        swapped[inormal] -= (idiff[inormal] != idiff[irecip])*2
        
        # compute the reciprocal error (using the magnitude of the measurements)
        reciprocalErr[inormal] = R[irecip] - R[inormal]*swapped[inormal]
        reciprocalErr[irecip] = R[irecip] - R[inormal]*swapped[inormal]
        reci_IP_err[inormal] = M[irecip] - M[inormal]  # TODO should we apply the swapped here too?!
        reci_IP_err[irecip] = M[inormal] - M[irecip]
        
        # compute reciprocal mean with all valid values
        ok1 = ~(np.isnan(R[inormal]) | np.isinf(R[inormal]))
        ok2 = ~(np.isnan(R[irecip]) | np.isinf(R[irecip]))
        
        ie = ok1 & ok2 # both normal and recip are valid
        reciprocalMean[inormal[ie]] = np.mean(np.c_[np.abs(R[inormal]),np.abs(R[irecip])], axis=1)
        reciprocalMean[irecip[ie]] = np.mean(np.c_[np.abs(R[inormal]),np.abs(R[irecip])], axis=1)
        
        ie = ok1 & ~ok2 # only use normal
        reciprocalMean[inormal[ie]] = np.abs(R[inormal[ie]])
        
        ie =  ~ok1 & ok2 # only use reciproal
        reciprocalMean[inormal[ie]] = np.abs(R[irecip[ie]])
        
        reciprocalErrRel = reciprocalErr / reciprocalMean
        
        reciprocalMean = np.sign(resist)*reciprocalMean # add sign
        
        ibad = np.array([np.abs(a) > 0.2 if ~np.isnan(a) else False for a in reciprocalErrRel])
        if self.debug:
            print('{:d}/{:d} reciprocal measurements found.'.format(np.sum(Ri != 0), len(Ri)))
            if np.sum(Ri != 0) > 0: # no need to display that if there is no reciprocal
                print('{:d} measurements error > 20 %'.format(np.sum(ibad)))        
                
        self.df['irecip'] = Ri
        self.df['reciprocalErrRel'] = reciprocalErrRel
        self.df['recipError'] = reciprocalErr
        self.df['recipMean'] = reciprocalMean
        self.df['reci_IP_err'] = reci_IP_err
        
        # in order to compute error model based on a few reciprocal measurements
        # we fill 'recipMean' column with simple resist measurements for lonely
        # quadrupoles (which do not have reciprocals)
        inotRecip = Ri == 0
        self.df.loc[inotRecip, 'recipMean'] = self.df.loc[inotRecip, 'resist']

        return Ri
    
    def computeReciprocalC(self): # fast vectorize version
        """Compute reciprocal measurements. using a bisection seach in cython! 
        Parameters
        ----------
        forceSign: bool, optional
            Force reciprocal and forward measurements to have the same 
            polarity regarding the calculation of the reciprocal errors. 
            Default is False. 
        
        Notes
        -----
        The method first sorts the dipole AB and MN. Then efficiently searches
        for reciprocal pairs with a bisection search. 
        """
        resist = np.asarray(self.df['resist'].values,dtype=float)
        phase = np.asarray(-self.kFactor*self.df['ip'].values,dtype=float) #converting chargeability to phase shift
        ndata = self.df.shape[0]
        array = np.asarray(self.isequence - 1, dtype=int)
        # if platform.system() == 'Windows':
        #     print('converting sequencing to long')
        #     array = np.asarray(self.isequence - 1, dtype=np.longlong)
        
        # print(type(phase))
        unpackme = fastReciprocals(array, resist, phase)
        irecip = unpackme[0] 
        reciprocalErr = unpackme[1] 
        reciprocalErrRel = unpackme[2] 
        reciprocalMean = unpackme[3] 
        reciprocalPhase = unpackme[4] 
        # ifwd = unpackme[5] 
        
        ibad = np.array([np.abs(a) > 0.2 if ~np.isnan(a) else False for a in reciprocalErrRel])
        if self.debug:
            print('{:d}/{:d} reciprocal measurements found.'.format(np.sum(irecip != 0), ndata))
            if np.sum(irecip != 0) > 0: # no need to display that if there is no reciprocal
                print('{:d} measurements error > 20 %'.format(np.sum(ibad)))        
                
        self.df['irecip'] = irecip
        self.df['reciprocalErrRel'] = reciprocalErrRel
        self.df['recipError'] = reciprocalErr
        self.df['recipMean'] = reciprocalMean
        self.df['reci_IP_err'] = reciprocalPhase
        
        return irecip 
    
    
    
    def showErrorDist(self, ax=None):
        """Calculate and plots reciprocal error probablity histogram.
        Good data will have a bell shape (normal) distribution where most datapoints have near
        zero reciprocal error.
        
        Parameters
        ----------
        ax : Matplotlib.Axes
            If specified, the graph will be plotted against it.
        """
        if ax is None:
            fig, ax = plt.subplots()
        
        percentError = 100*self.df['reciprocalErrRel'].replace([np.inf,-np.inf], np.nan).dropna() # nan and inf values must be removed
        ax.hist(percentError,bins=(np.arange(-100,100,0.5)),density=True,alpha=.3,label="Probability")
        errMax = percentError[np.abs(percentError) <= 100] # don't want to show data that has 1000% error
        errPercent = np.max(np.abs(percentError)) + 10 # limits the histogram's X axis
        if errPercent > 100:
            errPercent = 100
            
        try: # sometimes these don't work!
            parametricFit = norm.pdf(np.arange(-100,100,0.5),np.mean(errMax), np.std(errMax))
            KDEfit = gaussian_kde(errMax)
            ax.plot(np.arange(-100,100,0.5),parametricFit,'r--',label="Parametric fit")
            ax.plot(np.arange(-100,100,0.5), KDEfit(np.arange(-100,100,0.5)), 'blue',label="KDE fit")
        except:
            print('Parametric or KDE fit did not work.')
            pass
        
        ax.set_xlim(-1*(int(errPercent)),int(errPercent))
        ax.set_xlabel('Error [%]')
        ax.set_ylabel('Probability')
        ax.legend(loc='best', frameon=True)
        ax.set_title('Error probability distribution')
        
        if ax is None:
            return fig
    
    
    def filterDummy(self):
        """Remove measurements where abs(a-b) != abs(m-n) (likely to be dummy
        measurements added for speed).
        """
        # lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.isequence - 1 
        elecpos = self.elec['x'].values
        AB = np.abs(elecpos[array[:,0]]- elecpos[array[:,1]])
        MN = np.abs(elecpos[array[:,2]] - elecpos[array[:,3]])
        self.filterData(AB == MN)
        
    
    def filterRecip(self, percent=20, debug=True):
        """Filter measurements based on the level reciprocal error. 
        
        Parameters
        ----------
        percent : float, optional
            Percentage level of reciprocal error in which to filter the measurements.
            Percentage Errors > percentage will be removed. By default the value is 
            20.
        debug : bool, optional
            Print output to screen. Default is True. 
        """
        if all(np.isnan(self.df['recipError']) == True):
            raise ValueError("No reciprocal measurements present, cannot filter by reciprocal!")
        reciprocalErrRel = np.abs(self.df['reciprocalErrRel'].replace(np.nan, 0))
        igood = reciprocalErrRel < (percent/100) # good indexes to keep 
        df_temp = self.df.copy()
        # self.df = df_temp[igood] #keep the indexes where the error is below the threshold
        # self.dfPhaseReset = self.df.copy()
        self.filterData(igood)
        if debug:
            numRemoved = len(df_temp)-len(self.df)
            msgDump = "%i measurements with greater than %3.1f%% reciprocal error removed!" % (numRemoved, percent)
            print(msgDump)
            return numRemoved
    
    
    def filterStack(self, percent=2, debug=True):
        """Filter measurements based on the stacking error. 
        
        Parameters
        ----------
        percent : float, optional
            Percentage level of stacking error in which to filter the measurements.
            Percentage Errors > percentage will be removed. By default the value is 2.
        debug : bool, optional
            Print output to screen. Default is True. 
        """
        if 'dev' in self.df.columns:
            dev = self.df['dev'].replace(np.nan, 0)
            igood = dev < (percent)
            df_temp = self.df.copy()
            # self.df = df_temp[igood] #keep the indexes where the error is below the threshold
            # self.dfPhaseReset = self.df.copy()
            self.filterData(igood)
            if debug:
                numRemoved = len(df_temp)-len(self.df)
                msgDump = "%i measurements with greater than %3.1f%% stacking error removed!" % (numRemoved, percent)
                print(msgDump)
                return numRemoved
        else:
            raise ValueError("No stacking error column (dev) found!")
    
    def showInvError(self, index=0, ax=None):
        """Display inversion error by measurment numbers.
        
        Parameters
        ----------
        index : int, optional
            Index of survey (if time-lapse or batch). Default `index == 0`.
        ax : matplotlib axis
            If provided, the graph will be plotted against this axis.
        """
        errors = self.df['resInvError'].values
        errors = errors[~np.isnan(errors)]
        measurement_no = np.arange(1, len(errors) + 1)
        if ax is None:
            fig, ax = plt.subplots()
        ax.scatter(measurement_no, errors)
        ax.set_ylabel('Normalised Error')
        ax.set_xlabel('Measurement Number')
        # add diagnositic lines
        y_pos_limit = (3, 3)
        y_neg_limit = (-3, -3)
        baseline = (0, 0)
        ax.plot((1, measurement_no[-1]), y_pos_limit, 'r--')
        ax.plot((1, measurement_no[-1]), y_neg_limit, 'r--')
        ax.plot((1, measurement_no[-1]), baseline, '--')
        
    def filterInvError(self, vmin=None, vmax=None, debug=False):
        """Filter measurements by inversion error. 
        
        Parameters
        ----------
        vmin : float, optional
            Minimum error.
        vmax : float, optional
            Maximum error.
        debug : bool, optional
            Print output to screen. Default is False.
        """
        df = self.df.copy()
        if 'resInvError' not in df.columns:
            return
        else:
            resInvError = self.df['resInvError'].values
            if vmin is None:
                vmin = np.nanmin(resInvError)
            if vmax is None:
                vmax = np.nanmax(resInvError)
            ikeep = np.isnan(resInvError) | ((resInvError >= vmin) & (resInvError <= vmax))
            self.filterData(ikeep) # use the filter data function to ensure that the sequence is reset 
            
            if debug:
                numRemoved = len(df)-len(self.df)
                msgDump = "%i measurements outside [%s,%s] removed!" % (numRemoved, vmin, vmax)
                print(msgDump)
                return numRemoved
            
    
    def filterContRes(self,vmin=None, vmax=None, debug=True):
         """Filter measurements by contact resistances. 
         
         Parameters
         ----------
         vmin : float, optional
             Minimum value in unit of the cR column.
         vmax : float, optional
             Maximum value in unit of the cR column.
         debug : bool, optional
             Print output to screen. Default is True.
         """
         df = self.df.copy()
         if 'cR' not in df.columns:
             raise Exception('No contact resistance column available')
             
         cR = np.abs(self.df['cR'])
         if vmin is None:
             vmin = np.min(cR)
         if vmax is None:
             vmax = np.max(cR)
         ikeep = (cR >= vmin) & (cR <= vmax)
         # self.df = df[ikeep]
         # self.df.reset_index()
         self.filterData(ikeep)
         
         if debug:
             numRemoved = len(df)-len(self.df)
             msgDump = "%i measurements outside [%s,%s] ohm removed!" % (numRemoved, vmin, vmax)
             print(msgDump)
             return numRemoved
                
        
    def addFilteredIP(self):
        """Add filtered IP data after IP filtering and pre-processing. This is
        because the IP filtering is done on a different dataframe and only
        merged when called this method.
        """
        self.df = pd.merge(self.df, self.filterDataIP[['a','b','m','n']].copy(), how='inner', on=['a','b','m','n'])
        self.setSeqIds()
    
    @staticmethod
    def logClasses3(datax, datay, func, class1=None): # pragma: no cover
        """Perform a log class of datay based on datay and applied function
        func to each bin.
        
        Parameters
        ----------
        datax : array
            x values to be from which the log-classes will be made.
        datay : array
            y values to which the function `func` will be applied inside each
            log-class
        func : function
            Function to be applied to the y value of each log class.
        class1 : array, optional
            Array of values for the log classes. If given, the limit of the 
            classes will be computed as 10**class1
        
        Returns
        -------
        mbins : array
            x-means of each bin.
        vbins : array
            Value of each bin (output of the function `func`).
        nbins : array
            Number of measurements in each bin.
        """
        if class1 is not None:
            bins = 10.0**class1
        else:
            bins = 10.0**np.arange(-4.0,4.0,1.0)
        mbins = np.zeros(len(bins)-1)*np.nan # mean of each bin
        vbins = np.zeros(len(bins)-1)*np.nan # value of each bin after func applied
        nbins = np.zeros(len(bins)-1) # number of points in each bin
        # first remove nan from the datasets
        inan = ~np.isnan(datay)
        if np.sum(~inan) > 0: 
            print('logClasse3: ', np.sum(~inan), ' nan removed')
        datax = datax[inan]
        datay = datay[inan]
        datax = np.abs(datax) # take the ABS of datax
        for i in range(len(bins)-1):
            ie = (datax > bins[i]) & (datax <= bins[i+1])
            mbins[i] = np.mean([bins[i], bins[i+1]])
            vbins[i] = func(datay[ie])
            nbins[i] = np.sum(ie)
        # note : all data might not be in the bins, check with sum(nbins)
        return mbins, vbins, nbins
    
    
    def showError(self, ax=None):
        """Plot the reciprocal errors.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()
        reciprocalMean = self.df['recipMean'].abs().values
        reciprocalErr = self.df['recipError'].abs().values
        ax.loglog(reciprocalMean, reciprocalErr, 'o')
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')  
        ax.set_title('Observed Errors\n')
        if ax is None:
            return fig
        
        
    def showErrorIP(self, ax=None):
        """Plot the reciprocal phase discrepancies.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()
        temp_df_range_filter = self.df.copy().query('reci_IP_err>%s & reci_IP_err<%s' % (-self.phiCbarMax, self.phiCbarMax))
        reciprocalMean = np.abs(temp_df_range_filter['recipMean'].values)
        phase = np.abs(temp_df_range_filter['reci_IP_err'].values)
        ax.semilogx(reciprocalMean, phase, 'o')
        ax.set_xlabel(r'LogR [$\Omega$]')
        ax.set_ylabel(r's($\phi$) [mrad]')
        ax.set_title('Observed Discrepancies\n')
        if ax is None:
            return fig
        
    @staticmethod    
    def R_sqr(y, y_predict): # calculating R squared value to measure fitting accuracy
        rsdl = y - y_predict
        ss_res = np.sum(rsdl**2)
        ss_tot = np.sum((y-np.mean(y))**2)
        R2 = 1-(ss_res/ss_tot)
        R2 = np.around(R2,decimals=4)
        return R2
    
    # @staticmethod
    # def sign_coef(x):
    #     if x>=0:
    #         t = '+'
    #     else:
    #         t = '-'
    #     return t    
    
    @staticmethod
    def numFormating(numList):
        """Formats numbers between -1 to 1 based on their decimals.
        e.g., 0.00001 would be 1e-5 while 0.1 would remain 0.1
        
        Parameters
        ----------
        numList: list of floats
        
        Return
        ------
        formattedNums: string list of formatted floats
        """
        formattedNums = []
        for num in numList:
            if num != 0:
                if np.abs(num) <= 0.001:
                    formattedNums.append('{:.2e}'.format(num))
                else:
                    formattedNums.append('{:.3f}'.format(num))
            else:
                formattedNums.append('0')
        return formattedNums
                
    
    def fitErrorPwlIP(self, ax=None):
        """Plot the reciprocal phase errors with a power-law fit.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()
#        numbins_ip = 16
#        binsize_ip = int(len(self.df['reci_IP_err'])/numbins_ip)
        binsize_ip = 16 # default to 20 sample per bins
        numbins_ip = int(self.df.shape[0]/binsize_ip) # max 20 bins
        if numbins_ip > 20: # we want max 20 bins
            binsize_ip = int(len(self.df['reci_IP_err'])/20) # at least 20 samples per bin
            numbins_ip = 20
        Rn = np.abs(self.df['recipMean'])
        phasedisc = self.df['reci_IP_err']
        error_input_ip = (pd.concat((Rn,phasedisc),axis=1).rename(columns = {'recipMean':'absRn','reci_IP_err':'Phase_dicrep'})).sort_values(by='absRn').reset_index(drop = True).dropna().query('Phase_dicrep>%s & Phase_dicrep<%s' % (-self.phiCbarMax, self.phiCbarMax))# Sorting data based on R. the querry is based on input  phase range
        bins_ip = pd.DataFrame(np.zeros((numbins_ip,2))).rename(columns = {0:'R_mean',1:'Phi_dis_STD'})
        for i in range(numbins_ip): # bining 
            ns=i*binsize_ip
            ne=ns+binsize_ip-1
            bins_ip.iloc[i,0] = np.abs(error_input_ip['absRn'].iloc[ns:ne].mean())
            bins_ip.iloc[i,1] = error_input_ip['Phase_dicrep'].iloc[ns:ne].std()  
        bins_ip = bins_ip.dropna()
#        coefs_ip= np.linalg.lstsq(np.vstack([np.ones(len(bins_ip.iloc[:,0])), np.log(bins_ip.iloc[:,0])]).T, np.log(bins_ip.iloc[:,1]), rcond=None)[0] # calculating fitting coefficients (a,m)
        coefs_ip = polyfit(np.log(bins_ip.iloc[:,0]), np.log(bins_ip.iloc[:,1]), 1)[::-1]
        R_error_predict_ip = np.exp(coefs_ip[0])*(bins_ip.iloc[:,0]**coefs_ip[1]) # error prediction based of fitted power law model       
        ax.semilogx(error_input_ip['absRn'],np.abs(error_input_ip['Phase_dicrep']), '+', label = "Raw")
        ax.semilogx(bins_ip.iloc[:,0],bins_ip.iloc[:,1],'o',label="Bin Means")
        ax.plot(bins_ip.iloc[:,0],R_error_predict_ip,'r', label="Power Law Fit")
        ax.set_ylabel(r's($\phi$) [mrad]')
        ax.set_xlabel(r'$R_{avg}$ [$\Omega$]')      
        ax.legend(loc='best', frameon=True)
        R2_ip = self.R_sqr(np.log(bins_ip.iloc[:,1]),np.log(R_error_predict_ip))
        a1 = np.exp(coefs_ip[0])
        a2 = coefs_ip[1]
        
        self.df['phaseError'] = a1*(np.abs(self.df['recipMean'])**a2)
        self.df['phase'] = -self.kFactor*self.df['ip']
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x)**a2)
        self.phaseErrorModel = errorModel
        # formatting coefs for print
        formattedCoefs = self.numFormating([a1, a2, R2_ip])
        a1s = formattedCoefs[0]
        a2s = formattedCoefs[1]
        R2_ips = formattedCoefs[2]
        print ('Error model is: Sp(m) = {}*R^{} (R^2 = {})'.format(a1s,a2s,R2_ips))
        ax.set_title('Multi bin power-law phase error plot\n' + r's($\phi$) = {}$R^{{{}}}$ (R$^2$ = {})'.format(a1s,a2s,R2_ips))
        if ax is None:
            return fig   

        
    def fitErrorParabolaIP(self, ax=None):
        """Plot the reciprocal phase errors with a parabola fit.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()        
#        numbins_ip = 16
#        binsize_ip = int(len(self.df['reci_IP_err'])/numbins_ip) 
        binsize_ip = 16 # default to 20 sample per bins
        numbins_ip = int(self.df.shape[0]/binsize_ip) # max 20 bins
        if numbins_ip > 20: # we want max 20 bins
            binsize_ip = int(len(self.df['reci_IP_err'])/20) # at least 20 samples per bin
            numbins_ip = 20
        Rn = np.abs(self.df['recipMean'])
        phasedisc = self.df['reci_IP_err']
        error_input_ip = (pd.concat((Rn,phasedisc),axis=1).rename(columns = {'recipMean':'absRn','reci_IP_err':'Phase_dicrep'})).sort_values(by='absRn').reset_index(drop = True).dropna().query('Phase_dicrep>%s & Phase_dicrep<%s' % (-self.phiCbarMax, self.phiCbarMax))# Sorting data based on R. the querry is based on environmental IP
        bins_ip = pd.DataFrame(np.zeros((numbins_ip,2))).rename(columns = {0:'R_mean',1:'Phi_dis_STD'})
        for i in range(numbins_ip): # bining 
            ns=i*binsize_ip
            ne=ns+binsize_ip-1
            bins_ip.iloc[i,0] = np.abs(error_input_ip['absRn'].iloc[ns:ne].mean())
            bins_ip.iloc[i,1] = error_input_ip['Phase_dicrep'].iloc[ns:ne].std()  
        bins_ip = bins_ip.dropna()
        coefs_ip = polyfit(np.log10(bins_ip.iloc[:,0]), bins_ip.iloc[:,1], 2) # calculating fitting coefficients (a, b, c)
        R_error_predict_ip = (coefs_ip[0]*np.log10(bins_ip.iloc[:,0])**2) + (coefs_ip[1]*np.log10(bins_ip.iloc[:,0]) + coefs_ip[2] ) # error prediction based of fitted parabola model       
        ax.semilogx(error_input_ip['absRn'],np.abs(error_input_ip['Phase_dicrep']), '+', label = "Raw")
        ax.semilogx(bins_ip.iloc[:,0],bins_ip.iloc[:,1],'o',label="Bin Means")
        ax.semilogx(bins_ip.iloc[:,0],R_error_predict_ip,'r', label="Parabola Fit")
        ax.set_ylabel(r's($\phi$) [mrad]')
        ax.set_xlabel(r'$R_{avg}$ [$\Omega$]')      
        ax.legend(loc='best', frameon=True)
        R2_ip= self.R_sqr(bins_ip.iloc[:,1],R_error_predict_ip)
        a3 = coefs_ip[0]
        b3 = coefs_ip[1]
        c3 = coefs_ip[2]
        
        self.df['phaseError'] = (a3*np.log10(np.abs(self.df['recipMean']))**2) + (b3*np.log10(np.abs(self.df['recipMean'])) + c3)
        self.df['phase'] = -self.kFactor*self.df['ip']
        def errorModel(df):
            x = df['recipMean'].values
            return (a3*np.log10(np.abs(x))**2) + (b3*np.log10(np.abs(x)) + c3)
        self.phaseErrorModel = errorModel
        # formatting coefs for print
        formattedCoefs = self.numFormating([a3, b3, c3, R2_ip])
        a3s = formattedCoefs[0]
        b3s = formattedCoefs[1]
        c3s = formattedCoefs[2]
        R2_ips = formattedCoefs[3]
        ax.set_title('Multi bin parabola phase error plot\n' + r's($\phi$) = {}$R_{{avg}}^2$ + {}$R_{{avg}}$ + {} (R$^2$ = {})'.format(a3s, b3s, c3s, R2_ips))
        if ax is None:
            return fig   


    def fitErrorPwl(self, ax=None):
        """Fit an power law to the resistivity data.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()        
        if 'recipMean' not in self.df.columns:
            self.computeReciprocal()
        dfg = self.df[self.df['irecip'] > 0][['recipMean','recipError']].copy()
        dfg = dfg.replace([np.inf,-np.inf], np.nan).dropna()
        binsize = 20 # default to 20 sample per bins
        numbins = int(dfg.shape[0]/binsize) # max 20 bins
        if numbins > 20: # we want max 20 bins
            binsize = int(len(dfg['recipMean'])/20) # at least 20 samples per bin
            numbins = 20
        error_input = np.abs(dfg[['recipMean', 'recipError']]).sort_values(by='recipMean').reset_index(drop=True) # Sorting data based on R_avg
        error_input['recipError'] = error_input['recipError']
        bins = np.zeros((numbins,2))
        for i in range(numbins): # bining 
            ns=i*binsize
            ne=ns+binsize-1
            bins[i,0] = error_input['recipMean'].iloc[ns:ne].mean()
            bins[i,1] = error_input['recipError'].iloc[ns:ne].mean()    
#        print(bins)
#        print(np.sum(np.isnan(bins)))
#        print(np.sum(np.isinf(bins)))
#        coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]), rcond=None)[0] # calculating fitting coefficients (a,m)       
        coefs = polyfit(np.log(bins[:,0]), np.log(bins[:,1]), 1)[::-1] #order is of coefs is opposite to lstqd       
        R_error_predict = np.exp(coefs[0])*(bins[:,0]**coefs[1]) # error prediction based of power law model        
        ax.plot(np.abs(dfg['recipMean']),np.abs(dfg['recipError']), '+', label = "Raw")
        ax.plot(bins[:,0],bins[:,1],'o',label="Bin Means")
        ax.plot(bins[:,0],R_error_predict,'r', label="Power Law Fit")
        ax.set_xscale('log')
        ax.set_yscale('log')
        # lines above are work around to https://github.com/matplotlib/matplotlib/issues/5541/
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')      
        ax.legend(loc='best', frameon=True)
        R2 = self.R_sqr(np.log(bins[:,1]), np.log(R_error_predict))
        a1 = np.exp(coefs[0])
        a2 = coefs[1]

        self.df['resError'] = a1*(np.abs(self.df['recipMean'])**a2)
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x)**a2)
        self.errorModel = errorModel
        # formatting coefs for print
        formattedCoefs = self.numFormating([a1, a2, R2])
        a1s = formattedCoefs[0]
        a2s = formattedCoefs[1]
        R2s = formattedCoefs[2]
        if self.debug: 
            print('Error model is R_err = {} R_avg^{} (R^2 = {})'.format(a1s,a2s,R2s))
        ax.set_title('Multi bin power-law resistance error plot\n' + r'$R_{{error}}$ = {}$R_{{avg}}^{{{}}}$ (R$^2$ = {})'.format(a1s,a2s,R2s))
        if ax is None:
            return fig
    
    def fitErrorPwlEnv(self, ax=None):
        """Fit an power law to the resistivity data with an enveloppe fit.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()        
        if 'recipMean' not in self.df.columns:
            self.computeReciprocal()
        dfg = self.df[self.df['irecip'] > 0][['recipMean','recipError']].copy()
        dfg = dfg.replace([np.inf,-np.inf], np.nan).dropna()
        binsize = 20 # default to 20 sample per bins
        numbins = int(dfg.shape[0]/binsize) # max 20 bins
        if numbins > 20: # we want max 20 bins
            binsize = int(len(dfg['recipMean'])/20) # at least 20 samples per bin
            numbins = 20
        error_input = np.abs(dfg[['recipMean', 'recipError']]).sort_values(by='recipMean').reset_index(drop=True) # Sorting data based on R_avg
        error_input['recipError'] = error_input['recipError']
        bins = np.zeros((numbins,2))
        for i in range(numbins): # bining 
            ns=i*binsize
            ne=ns+binsize-1
            bins[i,0] = error_input['recipMean'].iloc[ns:ne].mean()
            bins[i,1] = error_input['recipError'].iloc[ns:ne].quantile(0.95)    
#        print(bins)
#        print(np.sum(np.isnan(bins)))
#        print(np.sum(np.isinf(bins)))
#        coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]), rcond=None)[0] # calculating fitting coefficients (a,m)       
        coefs = polyfit(np.log(bins[:,0]), np.log(bins[:,1]), 1)[::-1] #order is of coefs is opposite to lstqd       
        R_error_predict = np.exp(coefs[0])*(bins[:,0]**coefs[1]) # error prediction based of power law model        
        ax.plot(np.abs(dfg['recipMean']),np.abs(dfg['recipError']), '+', label = "Raw")
        ax.plot(bins[:,0],bins[:,1],'o',label="Bin Means")
        ax.plot(bins[:,0],R_error_predict,'r', label="Power Law Fit")
        ax.set_xscale('log')
        ax.set_yscale('log')
        # lines above are work around to https://github.com/matplotlib/matplotlib/issues/5541/
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')      
        ax.legend(loc='best', frameon=True)
        R2 = self.R_sqr(np.log(bins[:,1]), np.log(R_error_predict))
        a1 = np.exp(coefs[0])
        a2 = coefs[1]

        self.df['resError'] = a1*(np.abs(self.df['recipMean'])**a2)
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x)**a2)
        self.errorModel = errorModel
        # formatting coefs for print
        formattedCoefs = self.numFormating([a1, a2, R2])
        a1s = formattedCoefs[0]
        a2s = formattedCoefs[1]
        R2s = formattedCoefs[2]
        if self.debug: 
            print('Error model is R_err = {} R_avg^{} (R^2 = {})'.format(a1s,a2s,R2s))
        ax.set_title('Multi bin power-law resistance error plot\n' + r'$R_{{error}}$ = {}$R_{{avg}}^{{{}}}$ (R$^2$ = {})'.format(a1s,a2s,R2s))
        if ax is None:
            return fig
        
        
    def fitErrorLin(self, ax=None):
        """Fit a linear relationship to the resistivity data.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if ax is None:
            fig, ax = plt.subplots()        
        if 'recipMean' not in self.df.columns:
            self.computeReciprocal()
        dfg = self.df[self.df['irecip'] > 0][['recipMean','recipError']].copy()
        dfg = dfg.replace([np.inf,-np.inf], np.nan).dropna()
        binsize = 20 # default to 20 sample per bins
        numbins = int(dfg.shape[0]/binsize) # max 20 bins
        if numbins > 20: # we want max 20 bins
            binsize = int(len(dfg['recipMean'])/20) # at least 20 samples per bin
            numbins = 20
        error_input = np.abs(dfg[['recipMean', 'recipError']]).sort_values(by='recipMean').reset_index(drop=True) # Sorting data based on R_avg
        error_input['recipError'] = error_input['recipError']
        bins = np.zeros((numbins,2), dtype=float)
        for i in range(numbins): # bining 
            ns=i*binsize
            ne=ns+binsize-1
            bins[i,0] = error_input['recipMean'].iloc[ns:ne].mean()
            bins[i,1] = error_input['recipError'].iloc[ns:ne].mean()
        # coefs = polyfit(bins[:,0], bins[:,1], 1)
        slope, intercept, r_value, p_value, std_err = linregress(bins[:,0], bins[:,1])
        coefs = [slope, intercept]
        if coefs[1] < 0: # we don't want negative error -> doesn't make sense
#            x = bins[:,0][:,None]
#            slope, _, _, _ = np.linalg.lstsq(x, bins[:,1])
            slope = polyfit(bins[:,0], bins[:,1],0)
            coefs = [slope[0], 0]
            
 
        R_error_predict = ((coefs[0])*(bins[:,0]))+coefs[1] # error prediction based of linear model        
        ax.plot(error_input['recipMean'], error_input['recipError'], '+', label = "Raw")
        ax.plot(bins[:,0],bins[:,1],'o',label="Bin Means")
        ax.plot(bins[:,0],R_error_predict,'r', label="Linear Fit")
        ax.set_xscale('log')
        ax.set_yscale('log')
        # lines above are work around to https://github.com/matplotlib/matplotlib/issues/5541/
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')      
        ax.legend(loc='best', frameon=True)
        R2= self.R_sqr((bins[:,1]),(R_error_predict))
        a1 = coefs[0]
        a2 = coefs[1]
        
        self.df['resError'] = a1*(np.abs(self.df['recipMean']))+a2
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x))+a2
        self.errorModel = errorModel
#        self.errorModel = lambda x : a1*(np.abs(x))+a2
        # formatting coefs for print
        formattedCoefs = self.numFormating([a1, a2, R2])
        a1s = formattedCoefs[0]
        a2s = formattedCoefs[1]
        R2s = formattedCoefs[2]
        if self.debug: 
            print('Error model is R_err = {}*R_avg + {} (R^2 = {})'.format(a1s,a2s,R2s))
        ax.set_title('Multi bin linear resistance error plot\n' + r'$R_{{error}}$ = {}$R_{{avg}}$ + {} (R$^2$ = {})'.format(a1s,a2s,R2s))
        if ax is None:
            return fig                  
        
    
    def fitErrorLME(self, iplot=True, ax=None, rpath=None): # pragma: no cover
        """Fit a linear mixed effect (LME) model by having the electrodes as
        as grouping variables.
        
        Parameters
        ----------
        iplot : bool, optional
            If `True`, then a graph will be plotted.
        ax : matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis,
            otherwise a new figure will be created.
        rpath : str, optional
            Path of the directory with R (for Windows only).
        """
        # MATLAB code: lme4= fitlme(tbl,'recipErr~recipR+(recipR|c1)+(recipR|c2)+(recipR|p1)+(recipR|p2)'); 
        # requires R
        # statmodels or other python packages can't handle variable interaction yet

        OS = platform.system()
        
        if 'recipMean' not in self.df.columns:
            self.computeReciprocal()
        dfg = self.df[self.df['irecip'] > 0][['a','b','m','n','recipMean','recipError']].copy()
        dfg = dfg.replace([np.inf,-np.inf], np.nan).dropna()
        
        recipMean = np.abs(dfg['recipMean'].values)
        recipError = np.abs(dfg['recipError'].values)
        array = dfg[['a','b','m','n']].values.astype(int)        
        data = np.vstack([recipMean, recipError]).T
        data = np.hstack((array, data))
        df = pd.DataFrame(data, columns=['a','b','m','n','recipMean','obsErr'])  
        df.insert(0, 'idx', df.index + 1)
        df = df.astype({"idx": int, "a": int, "b": int, "m": int, "n": int})
        
        outputname = os.path.join(os.path.dirname(os.path.realpath(__file__)),'invdir','protocol-lmeInRecip.dat')
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(len(df)) + '\n')
            with open(outputname, 'a') as f:
                df.to_csv(f, sep='\t', header=False, index=False,float_format='%8.6e')

        outputname = os.path.join(os.path.dirname(os.path.realpath(__file__)),'invdir','protocol-lmeIn.dat')
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(len(self.df)) + '\n')
            with open(outputname, 'a') as f:
                self.df.to_csv(f, sep='\t', header=False, index=True,float_format='%8.6e',columns=['a','b','m','n','recipMean'])
                
        try:        
            if (OS == 'Windows') and (rpath is None):
                R_dir = input("Enter the directory where R.exe is installed: ")
                os.system('R_PATH'+R_dir)
            os.system('Rscript ' + os.path.join(os.path.dirname(os.path.realpath(__file__)),'lmefit.R'))  # Run R
            lmeError = protocolParserLME(os.path.join(os.path.dirname(os.path.realpath(__file__)),'invdir','protocol-lmeOutRecip.dat'))
            df['resError'] = lmeError # fitted results, only with results
            lmeError = protocolParserLME(os.path.join(os.path.dirname(os.path.realpath(__file__)),'invdir','protocol-lmeOut.dat'))
            self.df['resError'] = lmeError # predicted results, entire survey
            
#        df['resError'] = lmeError 
        
#        if 'recipMean' not in self.df.columns:
#            self.computeReciprocal()
#        dfg = self.df[self.df['irecip'] > 0]
#        
#        recipMean = np.abs(dfg['recipMean'].values)
#        recipError = np.abs(dfg['recipError'].values)
#        array = dfg[['a','b','m','n']].values.astype(int)
#        
#        data = np.vstack([recipMean, recipError]).T
#        data = np.hstack((data, array))
#        df = pd.DataFrame(data, columns=['recipMean','obsErr','a','b','m','n'])
#        md = smf.mixedlm('obsErr~recipMean', df, groups=df[['a','b','m','n']])
#        mdf = md.fit()
#        print(mdf.summary())
# 
#        def errorModel(df):
#            return mdf.predict(np.abs(df[['recipMean','a','b','m','n']]))
#        self.errorModel = errorModel
#        self.df['resError'] = self.errorModel(self.df)
 
        
            if ax is None:
                fig, ax = plt.subplots()
            ax.plot(df['obsErr'], df['resError'], 'o')
            ax.plot([np.min(df['obsErr']),np.max(df['obsErr'])], [np.min(df['obsErr']), np.max(df['obsErr'])], 'r-', label='1:1')
            ax.grid()
            ax.legend()
            ax.set_title('Linear Mixed Effect Model Fit')
            ax.set_xlabel(r'Reciprocal Error Observed [$\Omega$]')
            ax.set_ylabel(r'Reciprocal Error Predicted [$\Omega$]')
            ax.set_xscale('log')
            ax.set_yscale('log')
        except Exception as e:
            print('ERROR in Survey.lmefit(): Rscript command might not be available or the lme4 package is not installed.', e)

    
    def showHeatmap(self, ax=None):
        """Plot a phase heatmap (x = M, y = A and value = -phi) based on: 
            Orozco, A. F., K. H. Williams, and A. Kemna (2013), 
            Time-lapse spectral induced polarization imaging of stimulated uranium bioremediation, 
            Near Surf. Geophys., 11(5), 531–544, doi:10.3997/1873-0604.2013020)
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        filterDataIP_plotOrig = self.dfOrigin[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first').copy()
        if self.filt_typ == 'Raw':
            temp_heatmap_recip_filterN = self.dfOrigin[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first')
            dflen = len(self.dfOrigin)
        elif self.filt_typ == 'Filtered':
            if self.filterDataIP.empty:
                temp_heatmap_recip_filterN = self.df[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first')
                dflen = len(self.df)
            else:
                temp_heatmap_recip_filterN = self.filterDataIP[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first')
                dflen = len(self.filterDataIP)
        if self.protocolIPFlag == True:
            temp_heatmap_recip_filterN ['Phase'] = temp_heatmap_recip_filterN ['ip']*-1
        else:
            temp_heatmap_recip_filterN ['Phase'] = temp_heatmap_recip_filterN ['ip']*np.abs(self.kFactor)
        # because 'a' and 'm' are now string, we need to convert them back to int so that
        # the indexes are sorted properly (string sorting gives: '1', '10', '11', '2', ...)
        def tsort(val): # sort string of electrodes
            elec = [int(a.split()[-1]) for a in val]
            if len(val[0].split()) > 1:
                string = [int(a.split()[0]) for a in val]
            else:
                string = np.ones(len(val))
            return np.lexsort((elec, string))
        
        # if len(temp_heatmap_recip_filterN['a'].values[0].split()) > 1: # string is included
        #     temp_heatmap_recip_filterN['A'] = [int(a.split()[0])*100000 + int(a.split()[1]) for a in temp_heatmap_recip_filterN['a']]
        #     temp_heatmap_recip_filterN['M'] = [int(a.split()[0])*100000 + int(a.split()[1]) for a in temp_heatmap_recip_filterN['m']]
        # else:
        #     temp_heatmap_recip_filterN['A'] = temp_heatmap_recip_filterN['a'].astype(int)
        #     temp_heatmap_recip_filterN['M'] = temp_heatmap_recip_filterN['m'].astype(int)
        df = temp_heatmap_recip_filterN.set_index(['m','a']).Phase.unstack(0)
        isortA = tsort(df.index.values)
        isortM = tsort(df.columns.values)
        heat_recip_Filter = df.loc[df.index[isortA], df.columns[isortM]].copy()
        if ax is None:
            fig, ax = plt.subplots()  
        else:
            fig = ax.get_figure()             
        m = ax.imshow(heat_recip_Filter, origin='lower',cmap='jet',vmin=self.phiCbarmin, vmax=self.phiCbarMax)
        # ax.xaxis.set_ticks(np.arange(0,filterDataIP_plotOrig['a'].max()+1,4))
        # ax.yaxis.set_ticks(np.arange(0,filterDataIP_plotOrig['m'].max(),4))
        # ax.xaxis.set_ticks(np.arange(filterDataIP_plotOrig['a']))
        ax.set_yticks(np.arange(0, len(heat_recip_Filter.index), 4))
        ax.set_yticklabels(heat_recip_Filter.index[::4])
        ax.set_xticks(np.arange(0, len(heat_recip_Filter.columns), 4))
        ax.set_xticklabels(heat_recip_Filter.columns[::4])
        ax.set_ylabel('A')#,fontsize = 22)
        ax.set_xlabel('M')#,fontsize = 22)
#        ax.tick_params(labelsize=18)
        ax.set_title('%s\n%s measurements' % (self.filt_typ, dflen))#, fontsize=20)     
        if self.cbar == True:
            cbhnf = fig.colorbar(m, ax=ax)
            cbhnf.set_label(r'-$\phi$ [mrad]')#, fontsize=20)
#            cbhnf.ax.tick_params(labelsize=18)
        if ax is None:
            return fig
        
        
    
    def filterRangeIP(self, phimin, phimax):
        """Filter IP data according to a specified range.
        
        Parameters
        ----------
        phimin : float
            Minimium phase angle [mrad].
        phimax : float
            Maximum phase angle [mrad].
        """
        if self.filterDataIP.empty:
            if self.protocolIPFlag == True:
                self.filterDataIP = self.df.query('ip > %s and ip < %s' % (-phimax, -phimin))
            else:
                self.filterDataIP = self.df.query('ip > %s and ip < %s' % (phimin/np.abs(self.kFactor), phimax/np.abs(self.kFactor)))
        else:
            if self.protocolIPFlag == True:
                self.filterDataIP = self.filterDataIP.query('ip > %s and ip < %s' % (-phimax, -phimin))
            else:
                self.filterDataIP = self.filterDataIP.query('ip > %s and ip < %s' % (phimin/np.abs(self.kFactor), phimax/np.abs(self.kFactor)))
        self.addFilteredIP()
        
    
    def filterRecipIP(self):
        """Removing reciprocal measurements from dataset - only for visualization purposes on heatmap()
        """
        if self.filterDataIP.empty:
            self.filterDataIP = self.df.query('irecip>=0')
        else:
            self.filterDataIP = self.filterDataIP.query('irecip>=0')
        self.addFilteredIP()

    
    def filterNested(self):
        """Removes nested measurements:
            Where M or N are in between A and B
        """
        if self.filterDataIP.empty:
            temp_data = self.df.copy()
            mask = (temp_data.m < temp_data.b) & (temp_data.m > temp_data.a) | (temp_data.n < temp_data.b) & (temp_data.n > temp_data.a)
            temp_data.loc[mask, 'ip'] = np.nan
            self.filterDataIP = temp_data.dropna(subset = ['ip'])
        else:
            temp_data = self.filterDataIP.copy()
            mask = (temp_data.m < temp_data.b) & (temp_data.m > temp_data.a) | (temp_data.n < temp_data.b) & (temp_data.n > temp_data.a)
            temp_data.loc[mask, 'ip'] = np.nan
            self.filterDataIP = temp_data.dropna(subset = ['ip'])
        self.addFilteredIP()
        
    
    def filterNegative(self):
        """Remove negative apparent resistivty values
        """
        df = self.df.copy()
        keep_idx = df['resist']>0 # apparant resistivity values must be bigger than zero
        self.filterData(keep_idx)

    
    def showPseudo(self, ax=None, bx=None, threed=False, **kwargs):
        """Plot pseudo section if 2D survey or just quadrupoles transfer
        resistance otherwise.
        """
        if bx is None:
            bx = self.iBorehole
        if bx is False and threed is False:
            self._showPseudoSection(ax=ax, **kwargs)
        elif bx is False and threed is True:
            self._showPseudoSection3D(ax=ax, **kwargs)
        else:
            if ax is None:
                fig, ax = plt.subplots()
            ax.plot(self.df['resist'].values, '.')
            ax.set_xlabel('Measurements')
            ax.set_ylabel('Transfer Resistance [Ohm]')


    def showPseudoIP(self, ax=None, bx=None, threed=False, **kwargs):
        """Plot pseudo section if 2D survey or just quadrupoles phase otherwise.
        """
        if bx is None:
            bx = self.iBorehole
        if bx is False and threed is False:
            self._showPseudoSectionIP(ax=ax, **kwargs)
        elif bx is False and threed is True:
            self._showPseudoSection3D(ax=ax, column='ip', geom=False, **kwargs)
        else:
            if ax is None:
                fig, ax = plt.subplots()
            ax.plot(self.df['ip'].values, '.')
            ax.set_xlabel('Measurements')
            ax.set_ylabel('Phase [mrad]')
            
            
    def computeK(self):
        """
        Compute geometric factors, if any buried electrodes present then the 
        generalised function for buried electrodes will be used, otherwise 
        all electrodes are assumed to be on a flat surface. 

        Returns
        -------
        K: array like
            Geometric factors for each measurement in Survey.df (dataframe)

        """
        if self.isequence.shape[0] != len(self.df):
            self.setSeqIds()
        if 'buried' in self.elec.columns and any(self.elec['buried']):
            if self.debug:
                print('Computing geometric factors for buried electrodes!')
            return self.computeKborehole()
        else:
            return self.computeKsurface()
    
    
    def computeKsurface(self):
        """Compute geomatrix factor (assuming flat surface) and store it
        in self.df['K'].
        """
        array = self.isequence - 1 
        elec = self.elec[['x','y','z']].values
        
        # Looking for missing topography for remote electrodes
        # TODO: this might not be the right way of fixing this but helps
        if 'remote' in self.elec.columns: # 'remote' column exists in elec df
            iremote = self.elec['remote'].values
        else: # no remote column in elec df yet
            remote_flags = [-9999999, -999999, -99999,-9999,-999,
                        9999999, 999999, 99999] # values asssociated with remote electrodes
            iremote = np.in1d(elec[:,0], remote_flags)
            iremote = np.isinf(self.elec[['x','y','z']].values).any(1) | iremote
        if np.isnan(elec[iremote,2]).any():
            elec[iremote,2] = np.average(elec[~iremote,2])
        
        aposx = elec[:,0][array[:,0]]
        aposy = elec[:,1][array[:,0]]
        aposz = elec[:,2][array[:,0]]
        
        bposx = elec[:,0][array[:,1]]
        bposy = elec[:,1][array[:,1]]
        bposz = elec[:,2][array[:,1]]
        
        mposx = elec[:,0][array[:,2]]
        mposy = elec[:,1][array[:,2]]
        mposz = elec[:,2][array[:,2]]
        
        nposx = elec[:,0][array[:,3]]
        nposy = elec[:,1][array[:,3]]
        nposz = elec[:,2][array[:,3]]
        
        AM = np.sqrt((aposx-mposx)**2 + (aposy-mposy)**2 + (aposz-mposz)**2)
        BM = np.sqrt((bposx-mposx)**2 + (bposy-mposy)**2 + (bposz-mposz)**2)
        AN = np.sqrt((aposx-nposx)**2 + (aposy-nposy)**2 + (aposz-nposz)**2)
        BN = np.sqrt((bposx-nposx)**2 + (bposy-nposy)**2 + (bposz-nposz)**2)
        AM[AM == 0] = np.nan
        BM[BM == 0] = np.nan
        AN[AN == 0] = np.nan
        BN[BN == 0] = np.nan
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        
        self.df['K'] = K
        
    def computeKborehole(self,Gl = None): # ground level 
        """
        Compute geometric factor for a borehole survey assuming a flat 2D 
        surface. Gl = ground level. Calculated according to Eq 4.9 of Binley 
        and Slater (2020). 
        """
      
        # lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        # array = self.df[['a','b','m','n']].replace(lookupDict).values.astype(int)
        array = self.isequence - 1 
        elec = self.elec[['x','y','z']].values
                
        # Looking for missing topography for remote electrodes
        # TODO: this might not be the right way of fixing this but helps
        if 'remote' in self.elec.columns: # 'remote' column exists in elec df
            iremote = self.elec['remote'].values
        else: # no remote column in elec df yet
            remote_flags = [-9999999, -999999, -99999,-9999,-999,
                        9999999, 999999, 99999] # values asssociated with remote electrodes
            iremote = np.in1d(elec[:,0], remote_flags)
            iremote = np.isinf(self.elec[['x','y','z']].values).any(1) | iremote
        if np.isnan(elec[iremote,2]).any():
            elec[iremote,2] = np.average(elec[~iremote,2])
        
        if Gl is None: 
            Gl = np.max(elec[:,2]) # take maximum electrode coordinate elevation as ground level 
            
        Ax = elec[:,0][array[:,0]]
        Ay = elec[:,1][array[:,0]]
        Az = elec[:,2][array[:,0]]
        
        Bx = elec[:,0][array[:,1]]
        By = elec[:,1][array[:,1]]
        Bz = elec[:,2][array[:,1]]
        
        Mx = elec[:,0][array[:,2]]
        My = elec[:,1][array[:,2]]
        Mz = elec[:,2][array[:,2]]
        
        Nx = elec[:,0][array[:,3]]
        Ny = elec[:,1][array[:,3]]
        Nz = elec[:,2][array[:,3]]
        
        Ax_ = Ax.copy()
        Ay_ = Ay.copy()
        Az_ = Az + 2*(Gl-Az) # imaginary elevation 
        
        Bx_ = Bx.copy()
        By_ = By.copy()
        Bz_ = Bz + 2*(Gl-Bz)
        
        rAM = np.sqrt((Ax-Mx)**2 + (Ay-My)**2 + (Az-Mz)**2)
        rAN = np.sqrt((Ax-Nx)**2 + (Ay-Ny)**2 + (Az-Nz)**2)
        rBM = np.sqrt((Bx-Mx)**2 + (By-My)**2 + (Bz-Mz)**2)
        rBN = np.sqrt((Bx-Nx)**2 + (By-Ny)**2 + (Bz-Nz)**2)
        
        rAiM = np.sqrt((Ax_-Mx)**2 + (Ay_-My)**2 + (Az_-Mz)**2)
        rAiN = np.sqrt((Ax_-Nx)**2 + (Ay_-Ny)**2 + (Az_-Nz)**2)
        rBiM = np.sqrt((Bx_-Mx)**2 + (By_-My)**2 + (Bz_-Mz)**2)
        rBiN = np.sqrt((Bx_-Nx)**2 + (By_-Ny)**2 + (Bz_-Nz)**2)
        
        k = (1/rAM)+(1/rAiM)-(1/rAN)-(1/rAiN)-(1/rBM)-(1/rBiM)+(1/rBN)+(1/rBiN)
        K = 4*np.pi/k
        
        self.df['K'] = K 
        return K 
        
        
    def _computePseudoDepth(self,flag3d=False):
        """Compute pseudo-depths. Now returns pseudo depths with respect to the 
        electrode elevations. 
        
        Returns
        -------
        xpos, ypos, zpos all arrays containing position of the pseudo-section.
        """
        array = self.isequence - 1 
        elecm = self.elec[['x','y','z']].values.astype(float).copy() # electrode matrix - should be array of floats so np.inf work properly
        buried = self.elec['buried'].values 
        remote = self.elec['remote'].values
        
        ### first determine if measurements are nested ###
        #find mid points of AB 
        AB = (elecm[array[:,0]] + elecm[array[:,1]]) / 2 # mid points of AB 
        MN = (elecm[array[:,2]] + elecm[array[:,3]]) / 2 # mid points of MN 
        avgZ = np.mean(elecm[:,2][array],axis=1)

        ABrad = np.sqrt(np.sum((elecm[array[:,0]] - AB)**2,axis=1)) # radius of AB circle 
        MNrad = np.sqrt(np.sum((elecm[array[:,2]] - MN)**2,axis=1)) # radius of MN circle 
        
        Amn = np.sqrt(np.sum((elecm[array[:,0]] - MN)**2,axis=1)) # distance of A to mid point of MN 
        Bmn = np.sqrt(np.sum((elecm[array[:,1]] - MN)**2,axis=1)) # distance of B to mid point of MN 
        Nab = np.sqrt(np.sum((elecm[array[:,2]] - AB)**2,axis=1)) # distance of N to mid point of AB 
        Mab = np.sqrt(np.sum((elecm[array[:,3]] - AB)**2,axis=1)) # distance of M to mid point of AB
        
        iABinMN = (Amn < MNrad) & (Bmn < MNrad)
        iMNinAB = (Nab < ABrad) & (Mab < ABrad)
        inested = iABinMN | iMNinAB #if AB encompasses MN or MN encompasses AB 
                       
        # so it will never be taken as minimium
        elecm[remote,:] = np.inf
        
        # compute midpoint position of AB and MN dipoles
        elecx = elecm[:,0]
        elecy = elecm[:,1]

        #CURRENT ELECTRODE MIDPOINTS 
        caddx = np.abs(elecx[array[:,0]]-elecx[array[:,1]])/2
        caddy = np.abs(elecy[array[:,0]]-elecy[array[:,1]])/2
        caddx[np.isinf(caddx)] = 0 
        caddy[np.isinf(caddy)] = 0        
        cmiddlex = np.min([elecx[array[:,0]], elecx[array[:,1]]], axis=0) + caddx
        cmiddley = np.min([elecy[array[:,0]], elecy[array[:,1]]], axis=0) + caddy
        
        #POTENTIAL ELECTRODE MIDPOINTS
        paddx = np.abs(elecx[array[:,2]]-elecx[array[:,3]])/2
        paddy = np.abs(elecy[array[:,2]]-elecy[array[:,3]])/2
        paddx[np.isinf(paddx)] = 0 
        paddy[np.isinf(paddy)] = 0 
        pmiddlex = np.min([elecx[array[:,2]], elecx[array[:,3]]], axis=0) + paddx
        pmiddley = np.min([elecy[array[:,2]], elecy[array[:,3]]], axis=0) + paddy
    
        # for non-nested measurements
        xposNonNested  = np.min([cmiddlex, pmiddlex], axis=0) + np.abs(cmiddlex-pmiddlex)/2
        yposNonNested  = np.min([cmiddley, pmiddley], axis=0) + np.abs(cmiddley-pmiddley)/2
        pcdist = np.sqrt((cmiddlex-pmiddlex)**2 + (cmiddley-pmiddley)**2)

        # zposNonNested = np.sqrt(2)/2*pcdist
        zposNonNested = pcdist/4

        if np.all(cmiddley-pmiddley == 0):
            zposNonNested = 0.25*pcdist
        else: # for 3D arrays where there are mid-line measurements, this works closer to inversion results
            zposNonNested = np.sqrt(2)/2*pcdist 

        # for nested measurements use formula of Dalhin 2006
        xposNested = np.zeros(len(pmiddlex))
        yposNested = np.zeros(len(pmiddlex))
        outerElec1 = np.zeros((len(pmiddlex), 2)) # position of one electrode of outer dipole
        outerElec2 = np.zeros((len(pmiddlex), 2)) # position of one electrode of outer dipole
        # innerMid = np.zeros((len(pmiddlex), 2)) # middle of inner dipole
        if np.sum(iMNinAB) > 0:
            xposNested[iMNinAB] = pmiddlex[iMNinAB]
            yposNested[iMNinAB] = pmiddley[iMNinAB]
            outerElec1[iMNinAB] = np.c_[elecx[array[iMNinAB,0]], elecy[array[iMNinAB,0]]]
            outerElec2[iMNinAB] = np.c_[elecx[array[iMNinAB,1]], elecy[array[iMNinAB,1]]]

        if np.sum(iABinMN) > 0:
            xposNested[iABinMN] = cmiddlex[iABinMN]
            yposNested[iABinMN] = cmiddley[iABinMN]
            outerElec1[iABinMN] = np.c_[elecx[array[iABinMN,2]], elecy[array[iABinMN,2]]]
            outerElec2[iABinMN] = np.c_[elecx[array[iABinMN,3]], elecy[array[iABinMN,3]]]
      
        innerMid = np.c_[pmiddlex, pmiddley] # always use potential dipole
        
        apdist = np.sqrt(np.sum((outerElec1-innerMid)**2, axis=1))
        bpdist = np.sqrt(np.sum((outerElec2-innerMid)**2, axis=1))
        zposNested  = np.min([apdist, bpdist], axis=0)/3
        
        xpos = np.zeros_like(pmiddlex)
        ypos = np.zeros_like(pmiddlex)
        zpos = np.zeros_like(pmiddlex)
        
        xpos[~inested] = xposNonNested[~inested]
        xpos[inested] = xposNested[inested]
        
        ypos[~inested] = yposNonNested[~inested]
        ypos[inested] = yposNested[inested]
        
        zpos[~inested] = zposNonNested[~inested]
        zpos[inested] = zposNested[inested]
        
        # use geometric median instead for buried/unconventional electrodes
        special = np.asarray([False]*array.shape[0],dtype=bool)
        
        if any(buried): 
            elecx = elecm[:,0]
            elecy = elecm[:,1]
            elecz = elecm[:,2]
            for i in range(array.shape[0]):
                if not any(buried[array[i,:]]):
                    continue # skip if not including a buried electrode 
                if any(remote[array[i,:]]):
                    continue # skip if includes a infinite value, thats not good for this method. 
                special[i] = True 
                ex = elecx[array[i,:]]
                ey = elecy[array[i,:]]
                ez = elecz[array[i,:]]
                med = geometricMedian(ex,ey,ez)
                xpos[i] = med[0]
                ypos[i] = med[1]
                zpos[i] = med[2]
            
        if flag3d: 
            # find pseudo depth as if it they were below the surface 
            zpos[np.invert(special)] = avgZ[np.invert(special)] - zpos[np.invert(special)]
        else: 
            # positions will be plotted in terms of pseudo depth, so correct 
            zpos[special] = max(elecm[:,2]) - zpos[special]

        return xpos,ypos,zpos 
    
        
    def _showPseudoSection(self, ax=None, contour=False, log=False, geom=True,
                           vmin=None, vmax=None, column='resist', magFlag=False,
                           darkMode=False, cmap='viridis'):
        """Create a pseudo-section for 2D given electrode positions.
        
        Parameters
        ----------
        ax : matplotlib.Axes, optional
            If specified, the plot will be plotted agains this axis.
        contour : bool, optional
            If `True`, contour will be plotted. Otherwise, only dots.
        log : bool, optional
            If `True`, the log of the resistivity will be taken.
        geom : bool, optional
            If `True`, the geometric factor will be computed and applied
            to the transfer resisttance to obtain apparent resistiivty. If 
            `False`, only the Tx resistance will be plotted.
        vmin : float, optional
            Minimum value for the colorbar.
        vmax : float, optional
            Maximum value for the colorbar.
        magFlag: bool, optional
            If `True` then Tx resistance sign will be checked assuming a flat surface survey.
            `False`, resistance values are given with correct polarity.
        darkmode : bool, optional
            Alters coloring of the plot for a darker appearance
        cmap : str, optional
            Colormap.
        """
        resist = self.df[column].values.copy()
        xpos, _, ypos = self._computePseudoDepth()

        if magFlag: # in case magnitude is provided
            ie = self.checkTxSign(inplace=False)
            resist[ie] *= -1 
        
        if geom and column=='resist': # compute and applied geometric factor
            # self.computeK()
            resist = resist*self.df['K'].values 
            
        label = r'$\rho_a$ [$\Omega.m$]' # default label 
        title = 'Apparent Resistivity\npseudo section'
        if log:
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
            
        if column =='cR': # need to add more labels here probably 
            label = r'$\R_{cr}$ [$\Omega$]'
            title = 'Contact Resistance/nPseudo Section'
            
        if vmin is None:        
            vmin = np.percentile(resist[~np.isnan(resist)],10) # use 10% percentile 
        if vmax is None:
            vmax = np.percentile(resist[~np.isnan(resist)],90) # use 90% percentile 
                               
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
       
        levels = MaxNLocator().tick_values(vmin, vmax)
        # levels = np.linspace(vmin, vmax, 13)
        
        if contour:
            plotPsRes = ax.tricontourf(xpos, ypos, resist, levels=levels, extend='both', cmap=cmap)
        else:
            plotPsRes = ax.scatter(xpos, ypos, c=resist, s=70, vmin=vmin, vmax=vmax, cmap=cmap)#, norm=mpl.colors.LogNorm())
        cbar = fig.colorbar(plotPsRes, ax=ax, fraction=0.046, pad=0.04, 
                            label=label, ticks=levels)
        cbar.set_label(label)
            
        elecColor = 'k' if darkMode is False else 'w'
        if any(self.elec['buried'].values):
            ax.scatter(self.elec['x'], max(self.elec['z'])-self.elec['z'], c=elecColor)
            
        if log:
            val = ['{:.2f}'.format(10**lvl) for lvl in levels]
            cbar.ax.set_yticklabels(val)
        
        ax.invert_yaxis() # to remove negative sign in y axis    
        ax.set_title(title)
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Pseudo depth [m]')
        if ax is None:
            return fig
    
    def _showElecStrings3D(self,ax=None, 
                           strIdx=None, 
                           background_color=(0.8,0.8,0.8),
                           elec_color='k'):
        if ax is None: # make a plotter object if not already given 
            ax = pv.Plotter()
            ax.background_color = background_color
        else: # check the ax argument is for pyvista not matplotlib 
            typ_str = str(type(ax))
            if typ_str.find('pyvista') == -1:
                raise Exception('Error plotting with pyvista, expected a pyvista plotter object but got %s instead'%typ_str)
            ax.set_background(background_color)
            
        def lines_from_points(points):
            """Given an array of points, make a line set, 
            https://docs.pyvista.org/examples/00-load/create-spline.html"""
            poly = pv.PolyData()
            poly.points = points
            cells = np.full((len(points)-1, 3), 2, dtype=np.int_)
            cells[:, 1] = np.arange(0, len(points)-1, dtype=np.int_)
            cells[:, 2] = np.arange(1, len(points), dtype=np.int_)
            poly.lines = cells
            return poly
            
        elec = self.elec[['x','y','z']].values # make electrode numpy array 
        elec = elec[np.invert(self.elec['remote'].values),:] 
        
        if strIdx is not None: #add strings 
            if not isinstance(strIdx,list):
                raise TypeError('strIdx variable is not a list')
            for s in strIdx:
                line = lines_from_points(elec[s])
                tube = line.tube(radius=0.1)
                ax.add_mesh(tube,smooth_shading=True,color=(0.5,0.5,0.5))
        
        pvelec = pv.PolyData(elec)
        ax.add_mesh(pvelec, color=elec_color, point_size=10.,
                    render_points_as_spheres=True)
        
        
    def _showPseudoSection3D(self, ax=None, contour=False, log=False, geom=True,
                           vmin=None, vmax=None, column='resist', 
                           background_color=(0.8,0.8,0.8), elec_color='k',
                           strIdx=None, magFlag=False, darkMode=False, pvshow=True,
                           cmap='viridis'):
        """Create a pseudo-section for 3D surface array.
        
        Parameters
        ----------
        ax : matplotlib.Axes, optional
            If specified, the plot will be plotted against this axis.
        contour : bool, optional
            If `True`, contour will be plotted. Otherwise, only dots. Warning
            this is unstable. 
        log : bool, optional
            If `True`, the log of the resistivity (or desired column values) 
            will be taken.
        geom : bool, optional
            If `True`, the geometric factor will be computed and applied
            to the transfer resisttance to obtain apparent resistiivty. If 
            `False`, only the Tx resistance will be plotted.
        vmin : float, optional
            Minimum value for the colorbar.
        vmax : float, optional
            Maximum value for the colorbar.
        background_color: tuple, optional 
            background color for pyvista plotter object
        elec_color : tuple, string, optional
            color identifier for pyvista plotter object, determines color of 
            electrode points if they can be plotted. 
        strIdx : list, optional 
            returned from Project.detectStrings method. Each entry in list is an 
            array like of ints defining an electrode string. 
        magFlag : bool, optional
            If `True` then Tx resistance sign will be checked assuming a flat surface survey.
            `False`, resistance values are given with correct polarity.
        darkmode : bool, optional
            Alters coloring of pyvista plot for a darker appearance
        pvshow : bool, optional
            If True (default), the `Plotter.show()` is called. Set it to False
            to make 3D subplot with pyvista.
        cmap : str, optional
            Colormap.
        """
        if not pyvista_installed:
            print('pyvista not installed, cannot show 3D pseudo section')
            return
        
        #TODO contour is not working! # using Delaunay 3D instead.
        # if contour:
        #     contour = False
        
        # set colors for dark mode 
        tcolor = 'k'
        if darkMode:
            tcolor = 'w'
            elec_color = 'w'
            background_color = (0.2,0.2,0.2)
            
        if ax is None: # make a plotter object if not already given 
            # ax = BackgroundPlotter()
            ax = pv.Plotter()
            ax.background_color = background_color
        else: # check the ax argument is for pyvista not matplotlib 
            typ_str = str(type(ax))
            if typ_str.find('pyvista') == -1:
                raise Exception('Error plotting with pyvista, expected a pyvista plotter object but got %s instead'%typ_str)
            ax.set_background(background_color)
            
        # if intending to show apparent resistivity the safest option is to 
        # recompute these values 
        if column == 'app':
            column = 'resist'
            geom = True 
            
        resist = self.df[column].values
                
        if magFlag: # for cR3t and its magnitude calculation
            ie = self.checkTxSign(inplace=False) 
            resist[ie] *= -1 
            
        if geom and column=='resist': # compute and applied geometric factor
            self.computeK()
            resist = resist*self.df['K'].values
            resist[np.isinf(resist)] = np.nan # sometimes inf are generated
            # let's set them to nan to prevent colorscale to be meaningless
            column = 'app' # column is now apparent resistivity 
            
        # try and insert some more helpful colorbar labels over the internally used column names 
        if log and column == 'resist':
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = 'log10(Transfer resistance) [Ohm]'
        elif log and column == 'app':
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = 'log10(Apparent Resistivity) [Ohm.m]'
        if column == 'resist':
            label = 'Transfer resistance [Ohm]'
        elif column == 'app':
            label = 'Apparent Resistivity [Ohm.m]'
        elif column == 'ip':
            label = 'Phase [mrad]'
        elif log and column == 'cR':
            resist = np.log10(resist)
            label = 'log10(Contact Resistance) [Ohm]'
        elif column == 'cR':
            label = 'Contact Resistance [Ohm]'
        else:
            label = column 
            
        if vmin is None:
            vmin = np.percentile(resist[~np.isnan(resist)],10) # use 10% percentile 
        if vmax is None:
            vmax = np.percentile(resist[~np.isnan(resist)],90) # use 90% percentile 
            
        elecz = self.elec['z'].values
        if 'remote' in self.elec.columns: 
            elecz = elecz[np.invert(self.elec['remote'].values)]

        dp_x,dp_y,dp_z = self._computePseudoDepth(True) # get dipole depths 
        #Nb pseudo depths are returned as absolute values now 

        points = np.c_[dp_x,dp_y,dp_z]# convert to 3xn matrix. 

        pvpont = pv.PolyData(points)
        pvpont[label] = resist

        if not contour:
            ax.add_mesh(pvpont, point_size=10.,
                        #render_points_as_spheres=True,
                        cmap=cmap, #matplotlib colormap 
                        clim=[vmin,vmax], #color bar limits 
                        #show_scalar_bar=color_bar,#plot the color bar? 
                        #opacity=alpha,
                        scalar_bar_args={'color':tcolor,# 'interactive':True,
                                         'vertical':False})
        else:
            # using Delaunay 3D
            mesh = pvpont.delaunay_3d()
            color_map = plt.cm.get_cmap(cmap, 14) # subdividing colorbar so it look more like contouring!
            ax.add_mesh(mesh,
                        #render_points_as_spheres=True,
                        cmap=color_map, #matplotlib colormap 
                        clim=[vmin,vmax], #color bar limits 
                        #show_scalar_bar=color_bar,#plot the color bar? 
                        #opacity=alpha,
                        scalar_bar_args={'color':tcolor,# 'interactive':True,
                                         'vertical':False})            
                                 
        try:
            self._showElecStrings3D(ax=ax, strIdx = strIdx, 
                                    background_color=background_color,
                                    elec_color=elec_color)
        except Exception as e:
            print("Could not plot 3d electrodes, error = "+str(e))
            
        ax.show_grid(color=tcolor)
        if pvshow:
            ax.show()

    
    def _showPseudoSectionIP(self, ax=None, contour=False, vmin=None, vmax=None, 
                             darkMode=False, cmap='viridis'): #IP pseudo section
        """Create pseudo section of IP data with points (default)
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, the graph is plotted along this axis.
        contour : bool
            If True, use filled contour instead of points in the pseudosection.
        vmin : float, optional
            Miminum value for colorscale.
        vmax : float, optional
            Maximum value for colorscale.
        darkmode : bool, optional
            Alters coloring of the plot for a darker appearance
        cmap : str, optional
            Colormap.
            
        Returns
        -------
        fig : matplotlib figure
            If `ax` is not specified, the method returns a figure.
        """
        xpos, _, ypos = self._computePseudoDepth()
        
        if self.protocolIPFlag == True:
            ip = self.df['ip'].values
        else:
            ip = -self.kFactor*self.df['ip'].values

        label = r'$\phi$ [mrad]'
        
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        
        if contour is False:
            plotPsIP = ax.scatter(xpos, ypos, c=ip, s=70, vmin=vmin, vmax=vmax, cmap=cmap)
            cbar = fig.colorbar(plotPsIP, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(label)
        
        else:
            if vmin is None:
                vmin = np.min(ip)
            if vmax is None:
                vmax = np.max(ip)
            levels = MaxNLocator().tick_values(vmin, vmax)
            # levels = np.linspace(vmin, vmax, 13)
            plotPsIP = ax.tricontourf(xpos, ypos, ip, levels=levels, extend='both', cmap=cmap)
            fig.colorbar(plotPsIP, ax=ax, fraction=0.046, pad=0.04, label=label)
        
        elecColor = 'k' if darkMode is False else 'w'
        if any(self.elec['buried'].values):
            ax.scatter(self.elec['x'], max(self.elec['z'])-self.elec['z'], c=elecColor)
            
        ax.invert_yaxis() # to remove negative sign in y axis
        ax.set_title('Phase Shift\npseudo section')  
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Pseudo depth [m]')
        if ax is None:
            return fig
    
    
    def write2protocol(self, outputname='', err=False, errTot=False,
                       ip=False, res0=False, isubset=None, threed=False,
                       fm0=None):
        """Write a protocol.dat file for R2, cR2, R3t, cR3t.
        
        Parameters
        ----------
        outputname : str, optional
            Path of the output file.
        err : bool, optional
            If `True`, then the `resError` and `phaseError` (if IP present)
            will be used in the protocol.dat file.
        errTot : bool, optional
            If `True`, the modelling error will be added to the error from the
            error model to form the *total error*.
        ip : bool, optional
            If `True` and IP columns will be added to the file.
        res0 : bool, optional
            For time-lapse inversion, the background resistivity will be added 
            if `True`.
        isubset : array_like of bool, optional
            If specified, it will be used to take a subset of the original
            array. It can be used for difference inversion to take measurements
            in common between all the surveys.
        threed : bool, optional
            If `True`, it's for a 3D survey (and then add line numbers) or for
            a 2D survey (default).
       fm0 : numpy.array of float, optional
            Only for 3D timelapse, to compute d-d0+fm0 (if reg_mode == 2). If
            provided, will automatically set res0 to False.
            
        Returns
        -------
        protocol : pandas.DataFrame
            Dataframe which contains the data for the `protocol.dat`. 
        """
        # check if we need to take a subset of the dataframe (e.g. for timelapse)
        if isubset is None:
            # select half of paired and all non-paired quadrupoles
            ie = self.df['irecip'].values >= 0 # reciprocal + non-paired
            df = self.df[ie]
        else:
            df = self.df[isubset]
            # we need to take all measurements is subset is specified as some
            # quadrupoles might have reciprocal in one survey (irecip > 0) but
            # not in the next one (irecip = 0). So we take them all.
                                        
        
        # write quadrupoles
        x = df[['a','b','m','n']].values
        xx = np.c_[1+np.arange(len(x)), x]
        protocol = pd.DataFrame(xx, columns=['num','a','b','m','n'])
        
        # write transfer resistance
        # NOTE for IP, this is the magnitude not the resistance
        # the magnitude is always positive. This is the case as by the way we
        # compute 'recipMean' in self.computeReciprocal()
        if fm0 is not None: # 3D time-lapse with reg_mode == 2 -> d-d0+f(m0)
            res0 = False
            d = df['recipMean'].values
            d0 = df['recipMean0'].values
            fm0 = fm0 # provided as already from quads in common
            protocol['res'] = d-d0+fm0 # according to R3t manual v3.2
        else:
            protocol['res'] = df['recipMean'].values # non-paired will be nan
        
        # write background transfer resistance
        if res0 is True: # background for time-lapse and so
            protocol['res0'] = df['recipMean0'].values
        
        # write phase (and eventually convert chargeability into phase)
        if ip is True:
            if self.protocolIPFlag is True: # if imported from Protocol then it's already in phase
                protocol['phase'] = df['ip'].values 
            else: # otherwise it is in chargeability and we need to convert it
                protocol['phase'] = -self.kFactor*df['ip'].values # "-self.kFactor" is to change m to phi
                
        # write error for DC
        if err is True:
            if np.sum(np.isnan(df['resError'])) == 0: # no NaN inside
                protocol['resError'] = df['resError'].values
                if errTot == True: # we want to add modelling error to that
                    # print('Using total error')
                    if 'modErr' not in df.columns:
                        raise ValueError('ERROR : you must specify a modelling error')
                    else: # if present, compute Guassian propogation of errors (this is not geometric mean)
                        protocol['resError'] = np.sqrt(protocol['resError']**2 + df['modErr'].values**2)
            else:
                raise ValueError('You requested DC error but no error model can be found.')
                
        # write error for IP
        if (ip is True) and (err is True): # ip is present and we want error
            if np.sum(np.isnan(df['phaseError'])) == 0: # no NaN inside
                protocol['phaseError'] = df['phaseError'].values
            else:
                raise ValueError('You requested IP error but none can be found.')

                    
        # if it's 3D, we add the line number (all electrode on line 1)
        if threed:
            if len(protocol['a'].values[0].split()) == 1: # we don't have string number
                for c in ['a','b','m','n']: 
                    protocol[c] = '1 ' + protocol[c]
        
        # write protocol.dat
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(len(protocol)) + '\n')
            with open(outputname, 'a') as f:
                to_csv(protocol, f, sep='\t', header=False, index=False, lineterminator='\n')
        
        return protocol
        
        
    def filterDCA(self, dump=None):
        """Execute DCA filtering. Decay Curve Analysis (DCA) based on.
        Flores Orozco, A., Gallistl, J., Bücker, M., & Williams, K. H. (2017)., 
        Decay curve analysis for data error quantification in time-domain induced polarization imaging., 
        Geophysics, 83(2), 1–48. https://doi.org/10.1190/geo2016-0714.1
        
        Parameters
        ----------
        dump : function, optional
            Callback function to print the progress in percent.
        """
        if self.filterDataIP.empty:
            self.filterDataIP = DCA(self.df, dump=dump)
        else:
            self.filterDataIP = DCA(self.filterDataIP, dump=dump)
        self.addFilteredIP()
        
    
    def filterManual(self, attr='app', ax=None, log=False,
                     label=None, vmin=None, vmax=None, elec=True, 
                     darkMode=False, flag3d=False):
        """Manually filters the data visually. The points manually selected are
        flagged in the `Survey.iselect` vector and can subsequently be removed
        by calling `Survey.filterData(~Survey.iselect)`.
        
        Parameters
        ----------
        attr : str, optional
            Columns of `Survey.df` to use for plotting. Default is `app`
            (apparent resistivity).
        ax : matplotlib axis, optional
            If specified, the graph is plotted along the axis.
        log : bool, optional
            If `True` then all data will be log transformed.
        label : str, optional
            Label of the colorbar. If None, will be given from the 'attr' argument.
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        elec : bool, optional
            If `True`, the electrodes are shown and can be used for filtering.
        darkMode : bool, optional
            If true, electrodes wil be plotted in white, else black
        flag3d: bool, optional
            If true, function exits before plotting. 
        """
        if len(self.df) != self.isequence.shape[0]:
            self.setSeqIds()

        array = self.isequence -1
        if self.df.shape[0] == 0:
            raise ValueError('Unable to plot! Dataset is empty - can be due to filtering out all datapoints')
        
        percFact = 1
        if label is None:
            if attr == 'app':
                # self.computeK() # k is already computed 
                self.df['app'] = self.df['K']*self.df['resist']
                label = r'Apparent Resistivity [$\Omega.m$]'
            elif attr == 'resist':
                label = r'Transfer Resistance [$\Omega$]'
            elif attr == 'reciprocalErrRel':
                label = 'Reciprocal Error [%]'
                percFact = 100
            elif attr == 'dev':
                label = 'Stacking Error (dev) [%]'
            else:
                label = attr
                
        val = percFact*self.df[attr].values

        inan = np.zeros(len(val), dtype=bool)
        val = val.copy()[~inan]
        array = array.copy()[~inan]
        self.iselect = np.zeros(len(inan), dtype=bool)
        
        def setSelect(ie, boolVal):
            ipoints[ie] = boolVal
            self.iselect[~inan] = ipoints
        
        elecpos = self.elec['x'].values.copy()   
        
        elecpos[self.elec['remote']] = np.inf # so it will never be taken as minimium
        xpos, _, ypos = self._computePseudoDepth()
        self.eselect = np.zeros(len(elecpos), dtype=bool)
        
        if log:
            val = np.sign(val)*np.log10(np.abs(val))
        
        def onpick(event): # :pragma: no cover
            if lines[event.artist] == 'data':
                xid, yid = xpos[event.ind[0]], ypos[event.ind[0]]
                isame = (xpos == xid) & (ypos == yid)
                if (ipoints[isame] == True).all():
                    setSelect(isame, False)
                else:
                    setSelect(isame, True)
            if lines[event.artist] == 'elec':
                ie = (array == (event.ind[0]+1)).any(-1)
                if all(ipoints[ie] == True):
                    setSelect(ie, False)
                else:
                    setSelect(ie, True)
                if self.eselect[event.ind[0]] == True:
                    self.eselect[event.ind[0]] = False
                else:
                    self.eselect[event.ind[0]] = True
                elecKilled.set_xdata(elecpos[self.eselect])
                elecKilled.set_ydata(np.zeros(len(elecpos))[self.eselect])
            killed.set_xdata(x[ipoints])
            killed.set_ydata(y[ipoints])
            killed.figure.canvas.draw()           

        # in the ui we dont manually interact with the pseudo section, so the following is redundant 
        if flag3d:
            return 
                     
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        
        # put the numbers right next to the electrodes
        elecNumber = 1 + np.arange(len(elecpos))
#        [ax.text(a, 0, str(b), horizontalalignment='center', verticalalignment='bottom') for a,b in zip(elecpos[::5], elecNumber[::5])]
#        ax2 = ax
        
        # on the axis
        ax.invert_yaxis() # to remove negative sign in y axis
        ax2 = ax.twiny()
        if elec:
            ax.set_xlabel('Electrode number')
        ax.set_ylabel('Pseudo depth [m]')
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.set_xticks(elecpos[::5])
        ax.set_xticklabels(elecNumber[::5])
        
        elecColor = 'ko' if darkMode is False else 'wo'
        if elec:
            caxElec, = ax2.plot(elecpos, np.zeros(len(elecpos)), elecColor, picker=5)

        else:
            caxElec, = ax2.plot([], [], elecColor)
        
        cax = ax2.scatter(xpos, ypos, c=val, marker='o', picker=5, vmin=vmin,
                         vmax=vmax)
        cbar = fig.colorbar(cax, ax=ax2, fraction=0.046, pad=0.04)
        cbar.set_label(label)
        ax2.set_xlabel('Distance [m]')
        ax2.xaxis.set_label_position('bottom')
        ax2.xaxis.set_ticks_position('bottom')
        cax.figure.canvas.mpl_connect('pick_event', onpick)        
        killed, = cax.axes.plot([], [], 'rx')
        elecKilled, = cax.axes.plot([], [], 'rx')
        x = cax.get_offsets()[:,0]
        y = cax.get_offsets()[:,1]        
        ipoints = np.zeros(len(y),dtype=bool)
        lines = {cax:'data', caxElec:'elec', killed:'killed'}
        
        ax.set_xlim(ax2.get_xlim()) # here to get correct limits form ax2
        
        # put the numbers right next to the electrodes
#        elecNumber = 1 + np.arange(len(elecpos))
##        [ax.text(a, 0, str(b)) for a,b in zip(elecpos[::5], elecNumber[::5])]
#        
#        # on the axis
#        ax2 = ax.twiny()
#        ax2.set_xlabel('Electrode number')
#        ax2.set_xticks(elecpos[::5])
#        ax2.set_xticklabels(elecNumber[::5])
#        ax2.set_xlim(ax.get_xlim())

#        ax2xticks = [tick.get_text() for tick in ax.get_xticklabels()]
#        ax2xticks[1] = 1
#        ax2xticks[-2] = nelec
#        ax2.set_xticklabels(ax2xticks)


    
    def filterElec(self, elec=[], debug=True):
        """Filter out specific electrodes given.
        
        Parameters
        ----------
        elec : list
            List of electrode number to be removed.
        debug : bool, optional
            Print output to screen. Default is True.
        """
        ### old way, kinda slow ###
        # for e in elec:
        #     i2keep = (self.df[['a','b','m','n']].values != str(e)).all(1)
        #     self.filterData(i2keep)
        
        array = self.df[['a','b','m','n']].values # get schedule 
        arrayf = array.flatten() #flatten the array
        probef = np.zeros(len(arrayf),dtype=int) # probes to see if electrode is found
        if not isinstance(elec[0],str):#convert to string if not a string
            elec = [str(e) for e in elec]
        
        for i in range(len(arrayf)):#iterate through schedule matrix
            if arrayf[i] in elec:
                probef[i] = 1
        
        probe = probef.reshape(array.shape)
        check = np.sum(probe,axis=1) # keep where electrodes are not found
        i2keep = check == 0 # so where 0 counts are 
        self.filterData(i2keep) # send to filter data 

        if debug:
            numRemoved = np.sum(~i2keep)
            msgDump = '{:d} measurements removed!'.format(numRemoved)
            print(msgDump)
            return numRemoved
    
    
    def filterAppResist(self, vmin=None, vmax=None, debug=True):
        """Filter measurements by apparent resistivity for surface surveys 
        
        Parameters
        ----------
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        debug : bool, optional
            Print output to screen. Default is True.
        """
        df = self.df.copy()
        # self.computeK()
        appRes = self.df['K']*self.df['resist']
        if vmin is None:
            vmin = np.min(appRes)
        if vmax is None:
            vmax = np.max(appRes)
        ikeep = (appRes >= vmin) & (appRes <= vmax)
        self.df = df[ikeep]
        self.isequence = self.isequence[ikeep]
        self.df.reset_index()
        
        if debug:
            numRemoved = len(df)-len(self.df)
            msgDump = "%i measurements outside [%s,%s] removed!" % (numRemoved, vmin, vmax)
            print(msgDump)
            return numRemoved
    
    
    
    def filterTransferRes(self, vmin=None, vmax=None, debug=True):
        """Filter measurements by transfer resistance. 
        
        Parameters
        ----------
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        debug : bool, optional
            Print output to screen. Default is True.
        """
        df = self.df.copy()
        tRes = np.abs(self.df['resist'])
        if vmin is None:
            vmin = np.min(tRes)
        if vmax is None:
            vmax = np.max(tRes)
        ikeep = (tRes >= vmin) & (tRes <= vmax)
        self.df = df[ikeep]
        self.isequence = self.isequence[ikeep]
        self.df.reset_index()
        
        if debug:
            numRemoved = len(df)-len(self.df)
            msgDump = "%i measurements outside [%s,%s] removed!" % (numRemoved, vmin, vmax)
            print(msgDump)
            return numRemoved


    # def shuntIndexes(self, debug=True): # pragma: no cover
    #     """Normalise the indexes the sequence matrix to start at 1.
        
    #     Parameters
    #     ----------
    #     debug : bool, optional
    #         Set to True to print output.
    #     """
    #     df = self.df
    #     sch_mat = np.array((df['a'],df['b'],df['m'],df['n'])).T
    #     imin = np.min(sch_mat)
    #     if imin != 1:
    #         if debug:
    #             print("It looks like scheduling indexing starts at: %i"%imin)
    #             print("...normalising electrode indexing to start at 1.")
    #         #check to see if negative electrode numbers occur. 
    #         if imin < 0:
    #             #check there is a electrode at 0
    #             imin_pos = np.min(np.abs(sch_mat))
    #             if imin_pos == 1:
    #                 if debug:
    #                     print("Any positive electrode indexes now start at %i"%(abs(imin)+1))
    #                 crr_idx = np.argwhere(sch_mat>0)
    #                 sch_mat[crr_idx] -= 1
                    
    #         corrected = sch_mat - (imin - 1)
    #         #return corrected
    #         df['a'] = corrected[:,0]
    #         df['b'] = corrected[:,1]
    #         df['m'] = corrected[:,2]
    #         df['n'] = corrected[:,3]
    #         self.df = df 


    def swapIndexes(self, old_indx, new_indx): # pragma: no cover
        """Replace the electrode number in a sequence matrix with another.
        Survey dataframe is updated after running. 
        
        Parameters
        ----------
        old_idx : int
            Electrode number(or index) to be replaced 
        new_idx : int
            Replacement electrode number
        """
         
        df = self.df.copy()
        sch_mat = np.array((df['a'],df['b'],df['m'],df['n']),dtype=int).T
        idx = np.argwhere(sch_mat == old_indx)
        sch_mat[idx[:,0],idx[:,1]] = new_indx
        nmeas = len(df)
        a = ['0']*nmeas 
        b = ['0']*nmeas 
        m = ['0']*nmeas 
        n = ['0']*nmeas 
        for i in range(len(self.df)): # replace with strings 
            a[i] = str(sch_mat[i,0])
            b[i] = str(sch_mat[i,1])
            m[i] = str(sch_mat[i,2])
            n[i] = str(sch_mat[i,3])
        df['a'] = a
        df['b'] = b
        df['m'] = m
        df['n'] = n
        self.df = df 
        
    
    def normElecIdx(self, debug=True): # pragma: no cover
        """Normalise the electrode indexing sequencing to start at 1 and ascend
        consectively (ie 1 , 2 , 3 , 4 ... )
        
        Function firstly normalises all indexes so that the lowest electrode 
        number is 1. Then removes jumps in the electrode indexing.
        
        Parameters
        ----------
        debug : bool, optional
            Output will be printed to console if `True`. 
        """
        df = self.df.copy()

        sch_mat = np.array((df['a'],df['b'],df['m'],df['n'])).T
        uni_idx = np.unique(sch_mat.flatten()) # returns sorted and unique array of electrode indexes
        comp_idx = np.arange(1,len(uni_idx)+1,1) # an array of values order consectively 
        num_elec = len(uni_idx)  
        min_idx = np.min(uni_idx)
        surrogate = sch_mat.copy()
        
        if min_idx<1:
            print('Zero or Negative electrode indexes detected!')
        
        count = 0 # rolling total for number of indexes which had to be 'corrected' 
        for i in range(num_elec):
            if uni_idx[i] != comp_idx[i]:
                #we need to put the electrode order in sequence 
                off = comp_idx[i] - uni_idx[i]
                if off<0:
                    off+=1 # need to add one to aviod 0 indexed electrodes 
                idx_array = np.argwhere(sch_mat == uni_idx[i])
                new_id = uni_idx[i]+off
                for j in range(len(idx_array)):
                    idx = (idx_array[j][0],idx_array[j][1])
                    surrogate[idx] = new_id
                print('Electrode number %i changed to %i'%(uni_idx[i],new_id))
                count+=1
                
        if count>0: #only correct if chabges detected 
            df['a'] = surrogate[:,0]
            df['b'] = surrogate[:,1]
            df['m'] = surrogate[:,2]
            df['n'] = surrogate[:,3]
            
            self.df = df 
        if debug:
            if count > 0:
                print("%i electrode indexes corrected to be in consective and ascending order"%count)
            else:
                print("Electrode indexing appears to be okay")
                
    def normElecIdxwSeq(self,expected=None): # pragma: no cover
        """Normalise the electrode indexing sequencing to start at 1 and ascend
        consectively (ie 1 , 2 , 3 , 4 ... ). Also checks for any electrodes 
        which are missing out of sequence if an expected sequence is given. 
        
        Function firstly normalises all indexes so that the lowest electrode 
        number is 1. Then removes jumps in the electrode indexing.
        
        Parameters
        ----------
        expected : array like
            Expected sequence. 
        """
        df = self.df.copy()

        sch_mat = np.array((df['a'],df['b'],df['m'],df['n'])).T
        uni_idx = np.unique(sch_mat.flatten()) # returns sorted and unique array of electrode indexes
        
        if expected is None: 
            comp_idx = np.arange(1,len(uni_idx)+1,1) # an array of values order consectively 
        else:
            comp_idx = np.array(expected) # comparison array 
        exp_num_elec = len(comp_idx)#expected number of electrodes 
        min_idx = np.min(uni_idx)
        surrogate = sch_mat.copy()
        
        if min_idx<1:
            print('Zero or Negative electrode indexes detected!')
        
        count = 0 # rolling total for number of indexes which had to be 'corrected' 
        missing = []
        for i in range(exp_num_elec):
            check = uni_idx == comp_idx[i]
            if all(check==False)==True:# then the index is missing 
                print('electrode %i is missing from expected sequence'%comp_idx[i])
                missing.append(comp_idx[i])
        
        missing = np.array(missing)        
        crr_uni_idx = np.sort(np.append(uni_idx,missing))# corrected unique indexes     
        comp_idx = np.arange(1,len(crr_uni_idx)+1,1) # an array of values order consectively 

        for i in range(exp_num_elec):
            if crr_uni_idx[i] != comp_idx[i]:
                #we need to put the electrode order in sequence 
                off = comp_idx[i] - crr_uni_idx[i]
                if off<0:
                    off+=1 # need to add one to aviod 0 indexed electrodes 
                idx_array = np.argwhere(sch_mat == crr_uni_idx[i])
                new_id = crr_uni_idx[i]+off
                ignore=False
                if not len(idx_array)==0:
                    for j in range(len(idx_array)):
                        idx = (idx_array[j][0],idx_array[j][1])
                        surrogate[idx] = new_id
                else:
                    ignore=True
                print('Electrode number %i changed to %i'%(crr_uni_idx[i],new_id),end='')
                if ignore:
                    print(' (but will be ignored during inversion)')
                else:
                    print('')
                count+=1
                
        if count>0: #only correct if chabges detected 
            df['a'] = surrogate[:,0]
            df['b'] = surrogate[:,1]
            df['m'] = surrogate[:,2]
            df['n'] = surrogate[:,3]
            self.df = df 
        
        if count > 0:
            print("%i electrode indexes corrected to be in consective and ascending order"%count)
        else:
            print("Electrode indexing appears to be okay")    
    
    def elec2distance(self):
        """Convert 3d xy data in pure x lateral distance.
        Use for 2D data only!
        """
        elec = self.elec[['x','y','z']].values
        x = elec[:,0]
        y = elec[:,1]
        z = elec[:,2]
        
        idx = np.argsort(x) # sort by x axis
        x_sorted = x[idx]
        y_sorted = y[idx]
        z_sorted = z[idx]
        
        x_abs = np.zeros_like(x, dtype=float)
        #the first entry should be x = 0
        for i in range(len(x)-1):
            delta_x = x_sorted[i] - x_sorted[i+1]
            delta_y = y_sorted[i] - y_sorted[i+1]
            sq_dist = delta_x**2 + delta_y**2
            x_abs[i+1] =  x_abs[i] + np.sqrt(sq_dist)
        
        # return values in the order in which they came
        new_elec = np.zeros_like(elec, dtype=float)        
        for i in range(len(z)):
            put_back = idx[i] # index to where to the value back again
            new_elec[put_back,0] = x_abs[i]
            new_elec[put_back,2] = z_sorted[i]
    
        self.elec[['x','y','z']] = new_elec
        
    def elec2horidist(self):
        """
        Convert 2D xz data into a true horizontal distance. Assumes that survey 
        was done with a tape measure and the X distances are not true horizontal
        distance but distances measured along the ground. 

        """
        elec = self.elec[['x','y','z']].values
        x = elec[:,0]
        y = elec[:,1]
        z = elec[:,2]
        
        idx = np.argsort(x) # sort by x axis
        x_sorted = x[idx]
        z_sorted = z[idx]
        
        x_hor = np.zeros_like(x, dtype=float)
        #the first entry should be x = 0
        for i in range(len(x)-1):
            h = x_sorted[i] - x_sorted[i+1] # hypotenuse 
            dz = z_sorted[i] - z_sorted[i+1]
            sq = np.abs(h**2 - dz**2)
            x_hor[i+1] =  x_hor[i] + np.sqrt(sq)
        
        # return values in the order in which they came
        new_elec = np.zeros_like(elec, dtype=float)
        new_elec[:,1] = y         
        for i in range(len(z)):
            put_back = idx[i] # index to where to the value back again
            new_elec[put_back,0] = x_hor[i]
            new_elec[put_back,2] = z_sorted[i]
    
        self.elec[['x','y','z']] = new_elec
        
    def _rmLineNum(self):
        """Remove line numbers from dataframe.
        """
        print('removing line numbers')
        for i in range(len(self.df)):
            for col in ['a','b','m','n']:
                val = self.df[col][i]
                self.df.loc[i,col] = val.split()[-1]
                
    def _seq2mat(self): # :pragma: no cover
                
        """
        Read in ResIPy scheduling matrix and convert to 4 columns of integers 
        (abmn). 
    
        Returns
        -------
        schedule: nd array 
            N by 4 array of ints 
    
        """
        scheduledf = self.df.copy()
        nmeas = len(self.df)
        
        #go through and extract string and electrode id 
        scheduledf['as'] = np.zeros(nmeas,dtype=int)
        scheduledf['bs'] = np.zeros(nmeas,dtype=int)
        scheduledf['ms'] = np.zeros(nmeas,dtype=int)
        scheduledf['ns'] = np.zeros(nmeas,dtype=int)
        
        scheduledf['ai'] = np.zeros(nmeas,dtype=int)
        scheduledf['bi'] = np.zeros(nmeas,dtype=int)
        scheduledf['mi'] = np.zeros(nmeas,dtype=int)
        scheduledf['ni'] = np.zeros(nmeas,dtype=int)
        
        labels = ['a', 'b', 'm', 'n']
        for i in range(nmeas):
            for l in labels :
                e = scheduledf[l][i].split()
                scheduledf.loc[i,l+'s'] = int(e[0]) # pandas approach 
                scheduledf.loc[i,l+'i'] = int(e[1])
                # scheduledf[l+'s'][i] = int(e[0])
                # scheduledf[l+'i'][i] = int(e[1])
                
                
        stringid = scheduledf[['as','bs','ms','ns']].values
        elecid = scheduledf[['ai','bi','mi','ni']].values
        # stringid = np.zeros((nmeas,4),dtype=int)
        # elecid = np.zeros((nmeas,4),dtype=int)
        for i in range(4):
            stringid[:,i] = scheduledf[labels[i]+'s']
            elecid[:,i] = scheduledf[labels[i]+'i']
    
        stringidf = stringid.flatten()
        elecidf = elecid.flatten()
        strings = np.unique(stringid.flatten())
        correction = np.zeros_like(strings,dtype=int)
        addme = 0
        
        for i in range(len(strings)):
            s = strings[i]
            correction[i] = addme 
            addme += np.max(elecidf[stringidf==s])
            
        #correct the scheduling matrix 
        schedule = np.zeros((nmeas,4),dtype=int)
        
        for i in range(nmeas):
            c = 0
            for l in labels :
                e = scheduledf[l+'i'][i]
                s = scheduledf[l+'s'][i]-1
                schedule[i,c] = e + correction[s]
                c+=1
                
        return schedule 
    
    def addPerError(self,pnct=2.5):
        """Add a flat percentage error to resistivity data.
        
        Parameters
        ----------
        pnct: float
            Error in percent
        """
        res = np.array(self.df['resist'])
        error = (pnct/100)*res
        self.df['resError'] = error + np.array(self.df['resError'])
        
        
    def estimateError(self, a_wgt=0.001, b_wgt=0.02):
        """Estimate reciprocal error data for data with no reciprocals, following
        the same routine present in R2. This allows for the additional inclusion
        of modelling errors. 
        
        Parameters
        ----------
        a_wgt: float, optional
            a_wgt documented in the R2 documentation 
        b_wgt: float, optional 
            b_wgt documented in the R2 documentation  
        """
        res = np.array(self.df['resist'])
        var_res = (a_wgt*a_wgt)+(b_wgt*b_wgt) * (res*res)
        std_res = np.sqrt(var_res)
        self.df['resError'] = std_res
        return std_res 


    def exportSrv(self, fname=None):
        """Export .srv format for which is compatible with E4D. The e4d survey
        file includes the electrode locations, in addition to the scheduling 
        matrix. 
        
        Parameters
        ----------
        fname: string, optional
            Where the output file will be written to. By default the file will 
            take on the name of the survey and is written to the current working 
            directory. 
        """
        if fname is None: # rename output file name to that of the survey name
            fname = self.name + '.srv'
        fh = open(fname,'w')
        numelec = self.elec.shape[0] # number of electrodes 
        elec = self.elec[['x','y','z']].values
        buried = np.ones(elec.shape[0],dtype=int)#buried flag 
        buried[self.elec['buried']]=0
        fh.write('%i number of electrodes\n'%numelec)
        for i in range(numelec):
            line = '{:d} {:f} {:f} {:f} {:d}\n'.format(i+1,
                    elec[i,0],#x coordinate
                    elec[i,1],#y coordinate
                    elec[i,2],#z coordinate
                    buried[i])#buried flag 
            fh.write(line)
        #now write the scheduling matrix to file 
        
        #check for error column 
        check = np.isnan(self.df['resError'].values)
        if any(check==True):
        # if not 'resError' in self.df.columns: # the columns exists
            err = self.estimateError()
        else:
            err = self.df['resError'].values 
            
        if 'modErr' in self.df.columns:
            # if present, compute Guassian propogation of errors (this is not geometric mean)
            err = np.sqrt(err**2 + self.df['modErr'].values**2)
            
        ie = self.df['irecip'].values >= 0 # reciprocal + non-paired
        df = self.df[ie]
        nomeas = len(df) # number of measurements 
        df = df.reset_index().copy()
        err=err[ie]
        fh.write('\n%i number of measurements \n'%nomeas)
        # convert seqeunce to matrix 
        if len(self.df['a'][0].split()) == 2:
            seq = self._seq2mat()
        else:
            seq = self.df[['a','b','m','n']].values
        
        # format >>> m_indx a b m n V/I stdev_V/I
        for i in range(nomeas): 
            line = '{:d} {:d} {:d} {:d} {:d} {:f} {:f}\n'.format(i+1,
                    int(seq[i,0]),
                    int(seq[i,1]),
                    int(seq[i,2]),
                    int(seq[i,3]),
                    df['recipMean'][i],
                    err[i])
            fh.write(line)
        fh.close()

#%% deprecated methods
    def basicFilter(self): # pragma: no cover
        warnings.warn('This function is deprecated, use filterDefault() instead',
                      DeprecationWarning)
        self.filterDefault()
        
    def removeUnpaired(self): # pragma: no cover
        warnings.warn('This function is deprecated, use filterUnpaired() instead',
              DeprecationWarning)
        n = self.filterUnpaired()
        return n
    
    def estError(self, a_wgt=0.01, b_wgt=0.02): # pragma: no cover
        warnings.warn('The function is deprecated, use estimateError() instead.',
                      DeprecationWarning)
        self.estimateError(a_wgt=a_wgt, b_wgt=b_wgt)
            
    
    def filterdip(self, elec):  # pragma: no cover
        warnings.warn('The function is deprecated, use filterElec() instead.',
                      DeprecationWarning)
        index = (self.array == elec[0]).any(-1)
        for i in range(1,len(elec)):
            index = index | (self.array == elec[i]).any(-1)
        n = self.filterData(~index)
        return n


    def dca(self, dump=print): # pragma: no cover
        warnings.warn('The function is deprecated, use filterDCA() instead.',
                      DeprecationWarning)
        self.filterDCA(dump=dump)

        
    def manualFiltering(self, ax=None, figsize=(12,3), contour=False,
                        log=False, geom=False, label='', vmin=None, vmax=None):  # pragma: no cover
        warnings.warn('The function is deprecated, use filterManual() instead.',
                      DeprecationWarning)
        self.filterManual(ax=ax, figsize=figsize, contour=contour,
                          log=log, geom=geom, label=label, vmin=vmin, vmax=vmax)
 
    def pseudo(self, ax=None, bx=None, **kwargs): # pragma: no cover
        warnings.warn('The function is deprecated, use showPseudo() instead.',
                      DeprecationWarning)
        self.showPseudo(ax=ax, bx=bx, **kwargs)
    
    
    def pseudoIP(self, ax=None, bx=None, **kwargs): # pragma: no cover
        warnings.warn('The function is deprecated, use showPseudoIP() instead.',
                      DeprecationWarning)
        self.showPseudoIP(ax=ax, bx=bx, **kwargs)
    
    
    def reciprocal(self): # pragma: no cover
        warnings.warn('This function is deprecated, use computeReciprocal() instead',
              DeprecationWarning)
        out = self.computeReciprocal()
        return out
    
    
    def errorDist(self, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showErrorDist() instead',
              DeprecationWarning)
        self.showErrorDist(ax=ax)


    def removeDummy(self): # pragma: no cover
        warnings.warn('This function is deprecated, use filterDummy() instead',
              DeprecationWarning)
        n = self.filterDummy()
        return n
    
    
    def plotError(self, ax=None): # pragma: no cover
        warnings.warn('The function is deprecated, use showError() instead.',
                      DeprecationWarning)
        self.showError(ax=ax)
    
        
    def phaseplotError(self, ax=None): # pragma: no cover
        warnings.warn('The function is deprecated, use showErrorIP() instead.',
                      DeprecationWarning)
        self.showErrorIP(self, ax=ax)
        
        
    def pwlfitIP(self, ax=None): # pragma: no cover
        warnings.warn('The function is deprecated, use fitErrorPwlIP() instead.',
                      DeprecationWarning)
        self.fitErrorPwlIP(ax=ax)
        
    
    def plotIPFitParabola(self, ax=None): # pragma: no cover
        warnings.warn('The function is deprecated, use fitErrorParabolaIP() instead.',
                      DeprecationWarning)
        self.fitErrorParabolaIP(ax=ax)
              
    
    def pwlfit(self, ax=None): # pragma: no cover
        warnings.warn('The function is deprecated, use fitErrorPwl() instead.',
                      DeprecationWarning)
        self.fitErrorPwl(ax=ax)
        

    def lmefit(self, iplot=True, ax=None, rpath=None): # pragma: no cover
        warnings.warn('The function is deprecated, use fitErrorLME() instead.',
                      DeprecationWarning)
        self.fitErrorLME(iplot=iplot, ax=ax, rpath=rpath)
    
    
    def heatmap(self, ax=None): # pragma: no cover
        warnings.warn('The function is deprecated, use showHeatmap instead.',
                      DeprecationWarning)
        self.showHeatmap()
    
    
    def iprangefilt(self, phimin, phimax): # pragma: no cover
        warnings.warn('The function is deprecated, use showError() instead.',
                      DeprecationWarning)
        self.filterRangeIP(phimin, phimax)
    
    
    def removerecip(self): # pragma: no cover
        warnings.warn('The function is deprecated, use filterRecip() instead.',
                      DeprecationWarning)
        self.filterRecip(self)
    
    
    def removenested(self): # pragma: no cover
        warnings.warn('The function is deprecated, use filterNested() instead.',
                      DeprecationWarning)
        self.filterNested()
        
    def removeneg(self): # pragma: no cover
        warnings.warn('The function is deprecated, use filterNegative() instead.',
                      DeprecationWarning)
        self.filterNegative()
        
        
        
