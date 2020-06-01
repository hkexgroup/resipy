#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:21:54 2018

@author: ResIPy's core developers
"""
import sys
import os
import platform

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
from scipy.stats import norm
from scipy.stats.kde import gaussian_kde

from resipy.parsers import (syscalParser, protocolParserLME, resInvParser,
                     primeParserTab, protocolParser,
                     stingParser, ericParser, lippmannParser, aresParser,
                     srvParser, bertParser)
from resipy.DCA import DCA

# show the deprecation warnings
import warnings
warnings.simplefilter('default', category=DeprecationWarning)

try:#import pyvista if avaiable
    import pyvista as pv
    pyvista_installed = True
except ModuleNotFoundError:
    pyvista_installed = False
    # warnings.warn('pyvista not installed, 3D meshing viewing options will be limited')


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
    """
    def __init__(self, fname, ftype='', name='', spacing=None, parser=None, keepAll=True, debug=True):
        self.elec = pd.DataFrame()
        self.df = pd.DataFrame()
        if name == '':
            name = os.path.basename(os.path.splitext(fname)[0])
        self.name = name
        self.iBorehole = False # True is it's a borehole
        self.protocolIPFlag = False
        self.kFactor = 1
        self.errorModel = None # function instantiated after fitting an error model with reciprocal errors
        self.phaseErrorModel = None # idem for IP
        self.iselect = None # use in filterManual()
        self.eselect = None # idem
        self.debug = debug # plotting all information message by default
        if spacing is not None:
            warnings.warn('The spacing argument is deprecated and will be removed in the next version.',
                      DeprecationWarning)
            
        avail_ftypes = ['Syscal','ProtocolDC','Res2Dinv', 'BGS Prime', 'ProtocolIP',
                        'Sting', 'ABEM-Lund', 'Lippmann', 'ARES', 'E4D', 'BERT']# add parser types here! 

        
        if parser is not None:
            elec, data = parser(fname)
        else:
            if ftype == 'Syscal':
                elec, data = syscalParser(fname)
                self.kFactor = 1.2
            elif ftype =='ProtocolDC':
                elec, data = protocolParser(fname, ip=False)
            elif ftype == 'Res2Dinv':
                elec, data = resInvParser(fname)
            elif ftype == 'BGS Prime':
#                try:
#                    elec, data = primeParser(fname)
#                except:
                elec, data = primeParserTab(fname)
            elif ftype == 'ProtocolIP':
                elec, data = protocolParser(fname, ip=True)
                self.protocolIPFlag = True
            elif ftype == 'forwardProtocolDC':
                elec, data = protocolParser(fname, fwd=True)
            elif ftype == 'forwardProtocolIP':
                elec, data = protocolParser(fname, ip=True, fwd=True)
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
    #        elif (ftype == '') & (fname == '') & (elec is not None) and (data is not None):
    #            pass # manual set up
    #            print('Manual set up, no data will be imported')
            else:
                print("Unrecognised ftype, available types are :",avail_ftypes )
                raise Exception('Sorry this file type is not implemented yet')

        
        self.df = data.astype({'a':str, 'b':str, 'm':str, 'n':str})
        
        # add error measured to the error columns (so they can be used if no error model are fitted)
        if 'magErr' in self.df.columns:
            self.df['resError'] = self.df['magErr'].copy()
        else:
            #pass
            self.df['resError'] = np.nan
        if 'phiErr' in self.df.columns:
            self.df['phaseError'] = self.df['phiErr'].copy()
        else:
            self.df['phaseError'] = np.nan
            
        self.dfOrigin = self.df.copy() # unmodified
        if type(elec) == np.ndarray: # TODO for temporary stability
            dfelec = pd.DataFrame(elec.astype(float), columns=['x','y','z'])
            dfelec['remote'] = False
            dfelec['buried'] = False
            dfelec['label'] = (1 + np.arange(dfelec.shape[0])).astype(str)
        else:
            dfelec = elec
        self.elec = dfelec
        
        self.ndata = self.df.shape[0]
        self.phiCbarmin = 0
        self.phiCbarMax = 25
        self.filt_typ = 'Raw'
        self.cbar = True
        self.filterDataIP = pd.DataFrame()

        if ftype == 'BGS Prime':
            self.checkTxSign()
        if ftype == 'Syscal':
            if np.all(self.df['vp'] >= 0):
                self.checkTxSign()

        # convert apparent resistivity to resistance and vice versa
        self.computeK()
        if 'resist' in self.df.columns:
            self.df['app'] = self.df['K']*self.df['resist']
        elif 'app' in self.df.columns:
            self.df['resist'] = self.df['app']/self.df['K']
        
        # apply basic filtering
        self.filterDefault()
        self.computeReciprocal()
        self.dfReset = self.df.copy()
        self.dfPhaseReset = self.df.copy()
        
            
    @classmethod
    def fromDataframe(cls, df, elec):
        """Create a survey class from pandas dataframe.
        
        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe.
        elec : numpy.ndarray
            A nx3 numpy array witht the XYZ electrodes positions.
            
        Returns
        -------
        s_svy : Survey
            An instance of the Survey class.
        """
        surrogate_file = 'examples/dc-2d/protocol.dat'
        path2resipy = os.path.realpath(__file__).replace('\\src\resipy\\',surrogate_file)
        path = path2resipy # map to the protocol test file
        s_svy = cls(fname = path, ftype='Protocol')#surrogate survey
        s_svy.df = df
        dfelec = pd.DataFrame(elec, columns=['x','y','z'])
        dfelec['remote'] = False
        dfelec['buried'] = False
        dfelec['label'] = np.arange(dfelec.shape[0]).astype(str)
        s_svy.elec = dfelec
        s_svy.ndata = 1#len(data)
        s_svy.computeReciprocal()
        s_svy.filterDefault()
        return s_svy
    
    
    def __str__(self):
        out = "Survey class with {:d} measurements and {:d} electrodes".format(
            self.df.shape[0], self.elec.shape[0])
        return out
        
    
    def checkTxSign(self):
        """Checking the sign of the transfer resistances (2D survey only !).
        """
        resist = self.df['resist'].values.copy()
        self.computeK()
        K = self.df['K'].values
        ie = ((K < 0) & (resist > 0)) | ((K > 0) & (resist < 0))
        self.df.loc[ie, 'resist'] = resist[ie]*-1
        if 'recipMean' in self.df.columns: # in case this method is called after importation
            recipMean = self.df['recipMean'].values.copy()
            self.df.loc[ie, 'recipMean'] = recipMean[ie]*-1
        if 'recipMean0' in self.df.columns: # time-lapse: in case this method is called after importation
            recipMean0 = self.df['recipMean0'].values.copy()
            self.df.loc[ie, 'recipMean0'] = recipMean0[ie]*-1
        print('WARNING: change sign of ', np.sum(ie), ' Tx resistance.')
        

    def filterDefault(self):
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
        if np.sum(iout) > 0:
            dump('Survey.filterDefault: Number of Inf or NaN : {:d}\n'.format(np.sum(iout)))
            self.filterData(~iout)
        
        # remove duplicates
        shapeBefore = self.df.shape[0]
        self.df = self.df.drop_duplicates(subset=['a','b','m','n'], keep = 'first')
        ndup = shapeBefore - self.df.shape[0]
        if ndup > 0:
            dump('Survey.filterDefault: {:d} duplicates removed.\n'.format(ndup))
        
        # remove quadrupoles were A or B are also potential electrodes
        ie1 = self.df['a'].values == self.df['m'].values
        ie2 = self.df['a'].values == self.df['n'].values
        ie3 = self.df['b'].values == self.df['m'].values
        ie4 = self.df['b'].values == self.df['n'].values
        ie = ie1 | ie2 | ie3 | ie4
        if np.sum(ie) > 0:
            dump('Survey.filterDefault: {:d} measurements with A or B == M or N\n'.format(np.sum(ie)))
            self.filterData(~ie)
        
        # we need to redo the reciprocal analysis if we've removed duplicates and ...
        if ndup > 0 or np.sum(ie) > 0:
            self.computeReciprocal()
        
        # remove dummy for 2D case
#        if self.elec[:,1].sum() == 0: # it's a 2D case
#            self.removeDummy() # filter dummy by the rule if n < m then it's a dummy
        
        # create a backup of the clean dataframe
        self.dfReset = self.df.copy()
        self.dfPhaseReset = self.df.copy()
        
        ''' the following piece of code is not useful anymore. The default
        behavior is to keep all measurements except NaN, duplicates, Inf and
        dummy (2D only) and compute error model on the subset which has
        reciprocal then apply the error model to all the quadrupoles
        '''
        
#        # remove measurement without reciprocal
#        if self.keepAll is False:
#            print('ah ah let us make some order here !')
#            irecip = self.df['irecip'].values
#            self.filterData(irecip != 0) # non reciprocal (some are dummy)
#            if self.elec[:,1].sum() == 0: # it's a 2D case
#                self.removeDummy() # filter dummy by the rule if n < m then it's a dummy
#            if np.isnan(np.mean(self.df['recipError'])):# drop NaNs if present
#                self.df = self.df.dropna(subset = ['reciprocalErrRel','recipError','recipMean','reci_IP_err']) # NaN values in error columns cause crash in error analysis and final protocol outcome
#            self.dfReset = self.df.copy()
          
        
    def addData(self, fname, ftype='Syscal', parser=None):
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
            elif ftype == 'Res2Dinv':
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
            else:
                raise Exception('Sorry this file type is not implemented yet')
        self.df = self.df.append(data)
        if ftype == 'BGS Prime':
            self.checkTxSign()
        if ftype == 'Syscal':
            if np.all(self.df['vp'] >= 0):
                self.checkTxSign()
            
        self.dfOrigin = self.df.copy()
        self.ndata = len(self.df)
        self.computeReciprocal()
        self.filterDefault() # we assume the user input reciprocal data not another
        # normal survey

    
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
            self.ndata = len(i2keep)
            if 'irecip' in self.df.columns:
                # get a list of measurement that would be affected by the removal
                recip2reset = self.df[~i2keep]['irecip'].values*-1
            self.df = self.df[i2keep]
            if 'irecip' in self.df.columns:
                ie = np.in1d(self.df['irecip'].values, recip2reset)
                self.df.loc[ie, 'irecip'] = 0 # as their reciprocal is deleted, we set it to 0
                self.df.loc[ie, 'recipError'] = np.nan # they don't contribute to the error model anymore
                self.df.loc[ie, 'recipMean'] = self.df.loc[ie, 'resist'].values
            if self.debug is True:
                print('filterData:', np.sum(~i2keep), '/', len(i2keep), 'quadrupoles removed.')
            return np.sum(~i2keep)

    
    
    def filterUnpaired(self):
        """Remove quadrupoles that don't have a reciprocals. This might
        remove dummy measurements added for sequence optimization.
        """
        i2keep = self.df['irecip'] != 0
        print('removeUnpaired:', end='')
        self.filterData(i2keep)
        return np.sum(~i2keep)


    # def computeReciprocal2(self):
    #     """Compute reciprocal measurements.
        
    #     Notes
    #     -----
    #     The methods create an array where all positive index are normal
    #     measurements and the reciprocal measurements are their negative
    #     counterparts. This array is stored in the main dataframe `Survey.df`
    #     in the columns `irecip`. Measurements with `ìrecip=0` are measurements
    #     without reciprocal.
    #     """
    #     resist = self.df['resist'].values
    #     phase = -self.kFactor*self.df['ip'].values #converting chargeability to phase shift
    #     ndata = self.ndata
    #     array = self.df[['a','b','m','n']].values
        
    #     R = np.copy(resist)
    #     M = np.copy(phase)
    #     ndata = len(R)
    #     Ri = np.zeros(ndata)
    #     reciprocalErr = np.zeros(ndata)*np.nan
    #     reciprocalErrRel = np.zeros(ndata)*np.nan
    #     reciprocalMean = np.zeros(ndata)*np.nan
    #     reci_IP_err = np.zeros(ndata)*np.nan
    #     # search for reciprocal measurement
    #     count=1
    #     notfound=0
    #     for i in range(0,ndata):
    #         rev1=[2,3,0,1]
    #         rev2=[3,2,0,1]
    #         rev3=[2,3,1,0]
    #         rev4=[3,2,1,0]
    #         index1=(array[:,:] == array[i,rev1]).all(1)
    #         index2=(array[:,:] == array[i,rev2]).all(1)
    #         index3=(array[:,:] == array[i,rev3]).all(1)
    #         index4=(array[:,:] == array[i,rev4]).all(1)
    #         index=index1|index2|index3|index4
            
    #         if len(index[index]) == 1:
    #             reciprocalErr[index] = np.abs(R[i])-np.abs(R[index])
    #             reciprocalErrRel[index] = (np.abs(R[i])-np.abs(R[index]))/np.abs(R[i]) # in percent
    #             # flag first reciprocal found like this we can
    #             # delete the second one when we find it
    #             # (no loss of information)
    #             reci_IP_err[index] = M[i]-M[index]
    #             if Ri[i] == 0: # only if Ri == 0 otherwise we will
    #                 # overwrite all the data
    #                 Ri[i] = count # flag the first found
    #             if Ri[index] == 0:
    #                 Ri[index] = -count # flag its reciprocal
    #             # replace reciprocalMean by one measurements if the other one
    #             # is bad (NaN or Inf). Hopefully, the error model will find an
    #             # error to go with
    #             ok1 = ~(np.isnan(R[i]) | np.isinf(R[i]))
    #             ok2 = ~(np.isnan(R[index]) | np.isinf(R[index]))
    #             if ok1 & ok2:
    #                 reciprocalMean[i]=np.mean([np.abs(R[i]),np.abs(R[index])])
    #             elif ok1 & ~ok2:
    #                 reciprocalMean[i] = np.abs(R[i])
    #             elif ~ok1 & ok2:
    #                 reciprocalMean[i] = np.abs(R[index])
    #         else:
    #             #print("no reciprocal found for "+str(array[i,:]))
    #             notfound=notfound+1
    #         count=count+1
    #     if self.debug:
    #         print(str(notfound)+'/'+str(ndata)+' reciprocal measurements NOT found.')
    #     reciprocalMean = np.sign(resist)*reciprocalMean # add sign
    #     ibad = np.array([np.abs(a) > 0.2 if ~np.isnan(a) else False for a in reciprocalErrRel])
    #     if self.debug:
    #         print(str(np.sum(ibad)) + ' measurements error > 20 %')
        
    #     irecip = Ri        
        
    #     self.df['irecip'] = irecip
    #     self.df['reciprocalErrRel'] = reciprocalErrRel
    #     self.df['recipError'] = reciprocalErr
    #     self.df['recipMean'] = reciprocalMean
    #     self.df['reci_IP_err'] = reci_IP_err
    #     # in order to compute error model based on a few reciprocal measurements
    #     # we fill 'recipMean' column with simple resist measurements for lonely
    #     # quadrupoles (which do not have reciprocals)
    #     inotRecip = irecip == 0
    #     self.df.loc[inotRecip, 'recipMean'] = self.df.loc[inotRecip, 'resist']
        
    #     return Ri
    
    
    def computeReciprocal(self): # fast vectorize version
        """Compute reciprocal measurements.
        
        Notes
        -----
        The method first sorts the dipole AB and MN. Then creates a reciprocal
        quadrupole matrix. These matrices are then used with
        numpy.equal to produce a 2D matrix from which reciprocal are extracted.
        """
        resist = self.df['resist'].values
        phase = -self.kFactor*self.df['ip'].values #converting chargeability to phase shift
        ndata = self.ndata
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        
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
        
        reciprocalErr[inormal] = np.abs(R[irecip]) - np.abs(R[inormal])
        reci_IP_err[inormal] = M[irecip] - M[inormal]
        reciprocalErr[irecip] = np.abs(R[irecip]) - np.abs(R[inormal])
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
        parametricFit = norm.pdf(np.arange(-100,100,0.5),np.mean(errMax), np.std(errMax))
        KDEfit = gaussian_kde(errMax)
        ax.plot(np.arange(-100,100,0.5),parametricFit,'r--',label="Parametric fit")
        ax.plot(np.arange(-100,100,0.5), KDEfit(np.arange(-100,100,0.5)), 'k',label="KDE fit")
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
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        elecpos = self.elec['x'].values
        AB = np.abs(elecpos[array[:,0]]- elecpos[array[:,1]])
        MN = np.abs(elecpos[array[:,2]] - elecpos[array[:,3]])
        self.filterData(AB == MN)
        
    
    def filterRecip(self, percent=20, debug=True):
        """Filter measurements based on the level reciprocal error. 
        
        Parameters
        -----------
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
        self.df = df_temp[igood] #keep the indexes where the error is below the threshold
        self.dfPhaseReset = self.df.copy()
        if debug:
            numRemoved = len(df_temp)-len(self.df)
            msgDump = "%i measurements with greater than %3.1f%% reciprocal error removed!" % (numRemoved, percent)
            print(msgDump)
            return numRemoved
    
    
    def filterStack(self, percent=2, debug=True):
        """Filter measurements based on the stacking error. 
        
        Parameters
        -----------
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
            self.df = df_temp[igood] #keep the indexes where the error is below the threshold
            self.dfPhaseReset = self.df.copy()
            if debug:
                numRemoved = len(df_temp)-len(self.df)
                msgDump = "%i measurements with greater than %3.1f%% stacking error removed!" % (numRemoved, percent)
                print(msgDump)
                return numRemoved
        else:
            raise ValueError("No stacking error column (dev) found!")
        
    def addFilteredIP(self):
        """Add filtered IP data after IP filtering and pre-processing. This is
        because the IP filtering is done on a different dataframe and only
        merged when called this method.
        """
        self.df = pd.merge(self.df, self.filterDataIP[['a','b','m','n']].copy(), how='inner', on=['a','b','m','n'])

    
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
        coefs_ip = np.polyfit(np.log(bins_ip.iloc[:,0]), np.log(bins_ip.iloc[:,1]), 1)[::-1]
        R_error_predict_ip = np.exp(coefs_ip[0])*(bins_ip.iloc[:,0]**coefs_ip[1]) # error prediction based of fitted power law model       
        ax.semilogx(error_input_ip['absRn'],np.abs(error_input_ip['Phase_dicrep']), '+', label = "Raw")
        ax.semilogx(bins_ip.iloc[:,0],bins_ip.iloc[:,1],'o',label="Bin Means")
        ax.plot(bins_ip.iloc[:,0],R_error_predict_ip,'r', label="Power Law Fit")
        ax.set_ylabel(r's($\phi$) [mrad]')
        ax.set_xlabel(r'$R_{avg}$ [$\Omega$]')      
        ax.legend(loc='best', frameon=True)
        R2_ip= self.R_sqr(np.log(bins_ip.iloc[:,1]),np.log(R_error_predict_ip))
        a1 = np.exp(coefs_ip[0])
        a2 = coefs_ip[1]
        print ('Error model is: Sp(m) = {:.2f}*R^{:.2f} (R^2 = {:.2f})'.format(a1,a2,R2_ip))
        if a1 > 0.001:
            ax.set_title('Multi bin power-law phase error plot\n' + r's($\phi$) = {:.2f}$R^{{{:.3f}}}$ (R$^2$ = {:.3f})'.format(a1, a2, R2_ip))
        else:
            ax.set_title('Multi bin power-law phase error plot\n' + r's($\phi$) = {:.2e}$R^{{{:.3e}}}$ (R$^2$ = {:.3f})'.format(a1, a2, R2_ip))
        self.df['phaseError'] = a1*(np.abs(self.df['recipMean'])**a2)
        self.df['phase'] = -self.kFactor*self.df['ip']
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x)**a2)
        self.phaseErrorModel = errorModel
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
        coefs_ip = np.polyfit(np.log10(bins_ip.iloc[:,0]), bins_ip.iloc[:,1], 2) # calculating fitting coefficients (a, b, c)
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
        if a3 > 0.001:
            ax.set_title('Multi bin parabola phase error plot\n' + r's($\phi$) = {:.3f}$R_{{avg}}^2${:+.3f}$R_{{avg}}${:+.3f} ($R_{{avg}}^2$ = {:.3f})'.format(a3, b3, c3, R2_ip))
        else:
            ax.set_title('Multi bin parabola phase error plot\n' + r's($\phi$) = {:.2e}$R_{{avg}}^2${:+.2e}$R_{{avg}}${:+.2e} ($R_{{avg}}^2$ = {:.3f})'.format(a3, b3, c3, R2_ip))
        self.df['phaseError'] = (a3*np.log10(np.abs(self.df['recipMean']))**2) + (b3*np.log10(np.abs(self.df['recipMean'])) + c3)
        self.df['phase'] = -self.kFactor*self.df['ip']
        def errorModel(df):
            x = df['recipMean'].values
            return (a3*np.log10(np.abs(x))**2) + (b3*np.log10(np.abs(x)) + c3)
        self.phaseErrorModel = errorModel
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
        dfg = self.df[self.df['irecip'] > 0]
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
        coefs = np.polyfit(np.log(bins[:,0]), np.log(bins[:,1]), 1)[::-1] #order is of coefs is opposite to lstqd       
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
        R2= self.R_sqr(np.log(bins[:,1]),np.log(R_error_predict))
        a1 = np.exp(coefs[0])
        a2 = coefs[1]
#        a3 = np.exp(coefs[0])
#        a4 = coefs[1]
        print('Error model is R_err = {:.2f} R_avg^{:.3f} (R^2 = {:.4f})'.format(a1,a2,R2))
        if a1 > 0.001:
            ax.set_title('Multi bin power-law resistance error plot\n' + r'$R_{{error}}$ = {:.3f}$R_{{avg}}^{{{:.3f}}}$ (R$^2$ = {:.3f})'.format(a1,a2,R2))
        else:
            ax.set_title('Multi bin power-law resistance error plot\n' + r'$R_{{error}}$ = {:.2e}$R_{{avg}}^{{{:.3e}}}$ (R$^2$ = {:.3f})'.format(a1,a2,R2))
        self.df['resError'] = a1*(np.abs(self.df['recipMean'])**a2)
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x)**a2)
        self.errorModel = errorModel
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
        dfg = self.df[self.df['irecip'] > 0]
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
        coefs = np.polyfit(bins[:,0], bins[:,1], 1)
        if coefs[1] < 0: # we don't want negative error -> doesn't make sense
            x = bins[:,0][:,None]
            slope, _, _, _ = np.linalg.lstsq(x, bins[:,1])
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
        print('Error model is R_err = {:.2f}*R_avg + {:.2f} (R^2 = {:.4f})'.format(a1,a2,R2))
        if a1 > 0.001:
            ax.set_title('Multi bin linear resistance error plot\n' + r'$R_{{error}}$ = {:.3f}$R_{{avg}}${:+.3f} (R$^2$ = {:.3f})'.format(a1,a2,R2))
        else:
            ax.set_title('Multi bin linear resistance error plot\n' + r'$R_{{error}}$ = {:.2e}$R_{{avg}}${:+.2e} (R$^2$ = {:.3f})'.format(a1,a2,R2))
        self.df['resError'] = a1*(np.abs(self.df['recipMean']))+a2
        def errorModel(df):
            x = df['recipMean'].values
            return a1*(np.abs(x))+a2
        self.errorModel = errorModel
#        self.errorModel = lambda x : a1*(np.abs(x))+a2
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
        dfg = self.df[self.df['irecip'] > 0]
        
        recipMean = np.abs(dfg['recipMean'].values)
        recipError = np.abs(dfg['recipError'].values)
        array = dfg[['a','b','m','n']].values.astype(int)        
        data = np.vstack([recipMean, recipError]).T
        data = np.hstack((array,data))
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
        temp_heatmap_recip_filterN['A'] = temp_heatmap_recip_filterN['a'].astype(int)
        temp_heatmap_recip_filterN['M'] = temp_heatmap_recip_filterN['m'].astype(int)
        heat_recip_Filter = temp_heatmap_recip_filterN.set_index(['M','A']).Phase.unstack(0)     
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
        ax.grid(False)
        if self.cbar==True:
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
            self.showPseudoSection3D(ax=ax, **kwargs)
        else:
            if ax is None:
                fig, ax = plt.subplots()
            ax.plot(self.df['resist'].values, '.')
            ax.set_xlabel('Measurements')
            ax.set_ylabel('Transfer Resistance [Ohm]')


    def showPseudoIP(self, ax=None, bx=None, **kwargs):
        """Plot pseudo section if 2D survey or just quadrupoles phase otherwise.
        """
        if bx is None:
            bx = self.iBorehole
        if bx is False:
            self._showPseudoSectionIP(ax=ax, **kwargs)
        else:
            if ax is None:
                fig, ax = plt.subplots()
            ax.plot(self.df['ip'].values, '.')
            ax.set_xlabel('Measurements')
            ax.set_ylabel('Phase [mrad]')
    
        
    def computeK(self):
        """Compute geomatrix factor (assuming flat 2D surface) and store it
        in self.df['K'].
        """
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        elec = self.elec[['x','y','z']].values
        
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
        
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        
        self.df['K'] = K
        
        
    def _showPseudoSection(self, ax=None, contour=False, log=False, geom=True,
                           vmin=None, vmax=None, column='resist'):
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
        """
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        elecpos = self.elec['x'].values.copy() # we don't want the x values become np.inf in remote situation as it'll mess up future computeK()
        resist = self.df[column].values.copy()
        
        if geom: # compute and applied geometric factor
            self.computeK()
            resist = resist*self.df['K']

        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
            
        if log:
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
        else:
            label = r'$\rho_a$ [$\Omega.m$]'
                       
        elecpos[self.elec['remote'].values] = np.inf # so it will never be taken as minimium
        
        cadd = np.abs(elecpos[array[:,0]]-elecpos[array[:,1]])/2
        cadd[np.isinf(cadd)] = 0 # they are inf because of our remote
        cmiddle = np.min([elecpos[array[:,0]], elecpos[array[:,1]]], axis=0) + cadd
        
        padd = np.abs(elecpos[array[:,2]]-elecpos[array[:,3]])/2
        padd[np.isinf(padd)] = 0
        pmiddle = np.min([elecpos[array[:,2]], elecpos[array[:,3]]], axis=0) + padd

        xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
        ypos = np.sqrt(2)/2*np.abs(cmiddle-pmiddle)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
       
        if contour is False:
            plotPsRes = ax.scatter(xpos, ypos, c=resist, s=70, vmin=vmin, vmax=vmax)#, norm=mpl.colors.LogNorm())
            cbar = fig.colorbar(plotPsRes, ax=ax, fraction=0.046, pad=0.04, label=label)
            cbar.set_label(label)

        if contour:
            if vmin is None:
                vmin = np.min(resist)
            if vmax is None:
                vmax = np.max(resist)
            levels = np.linspace(vmin, vmax, 13)
            plotPsRes = ax.tricontourf(xpos, ypos, resist, levels = levels, extend = 'both')
            fig.colorbar(plotPsRes, ax=ax, fraction=0.046, pad=0.04, label=label)
        
        ax.invert_yaxis() # to remove negative sign in y axis    
        ax.set_title('Apparent Resistivity\npseudo section')
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Pseudo depth [m]')
        if ax is None:
            return fig
        
        
    def showPseudoSection3D(self, ax=None, contour=False, log=False, geom=True,
                           vmin=None, vmax=None, column='resist', 
                           background_color=(0.8,0.8,0.8), elec_color='k'):
        """Create a pseudo-section for 3D surface array.
        
        Parameters
        ----------
        ax : matplotlib.Axes, optional
            If specified, the plot will be plotted agains this axis.
        contour : bool, optional
            If `True`, contour will be plotted. Otherwise, only dots. Warning
            this is unstable. 
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
        """
        if not pyvista_installed:
            print('pyvista not installed, cannot show 3D psuedo section')
            return 
            
        if ax is None: # make a plotter object if not already given 
            ax = pv.BackgroundPlotter()
            ax.background_color = background_color
        else: # check the ax argument is for pyvista not matplotlib 
            typ_str = str(type(ax))
            if typ_str.find('pyvista') == -1:
                raise Exception('Error plotting with pyvista, show3D (meshTools.py) expected a pyvista plotter object but got %s instead'%typ_str)
            ax.set_background(background_color)
            
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        elec = self.elec[['x','y','z']].values
        resist = self.df[column].values
        
        if geom: # compute and applied geometric factor
            self.computeK()
            resist = resist*self.df['K']

        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
            
        if log:
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
        else:
            label = r'$\rho_a$ [$\Omega.m$]'

        if vmin is None:
            vmin = np.min(resist)
        if vmax is None:
            vmax = np.max(resist)
        
        if self.elec['remote'].sum() > 0: # how to deal with remote electrodes?!! 
            raise Exception('Remote electrodes not currently supported for 3D pseudo sections')
        #     elec = elec[self.iremote!=True,:] # ignore the 
            
        #might be a better way to do this than looping, for example caculating the euclidian matrix (can be RAM limiting)
        #or use SciPy's KDTree
        def find_dist(elec_x,elec_y,elec_z): # find maximum and minimum electrode spacings 
            dist = np.zeros((len(elec_x),len(elec_x)))   
            x1 = np.array(elec_x)
            y1 = np.array(elec_y)
            z1 = np.array(elec_z)
            for i in range(len(elec_x)):
                x2 = elec_x[i]
                y2 = elec_y[i]
                z2 = elec_z[i]
                dist[:,i] = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            return dist.flatten() # array of all electrode distances 
        
        #loop to work out the 3D setup of the apparent resistivity measurements 
        nmeas = array.shape[0]
        dp_x = np.zeros(nmeas)
        dp_y = np.zeros(nmeas)
        dp_z = np.zeros(nmeas)
        for i in range(nmeas):
            dp_elec = elec[array[i]-1,:]
            dp_dist = np.max(find_dist(dp_elec[:,0], dp_elec[:,1], dp_elec[:,2]))
            dp_x[i] = np.mean(dp_elec[:,0])
            dp_y[i] = np.mean(dp_elec[:,1])
            dp_z[i] = np.mean(dp_elec[:,2]) - dp_dist/2
        
        points = np.array([dp_x,dp_y,dp_z]).T
        pvpont = pv.PolyData(points)
        pvpont[label] = resist
        
        if not contour:
            ax.add_mesh(pvpont, point_size=10.,
                        #render_points_as_spheres=True,
                        #cmap=color_map, #matplotlib colormap 
                        clim=[vmin,vmax], #color bar limits 
                        #show_scalar_bar=color_bar,#plot the color bar? 
                        #opacity=alpha,
                        scalar_bar_args={'color':'k',# 'interactive':True,
                                         'vertical':False,
                                         'title_font_size':16,
                                         'label_font_size':14})
        else:
            warnings.warn('3D contours are currently not stable!')
            ax.add_mesh(pvpont.outline())
            levels = np.linspace(vmin, vmax, 13)
            contrmesh = pvpont.contour(levels,scaler=resist)
            ax.add_mesh(contrmesh,
                        #cmap=color_map, #matplotlib colormap 
                        clim=[vmin,vmax], #color bar limits 
                        #show_scalar_bar=color_bar,#plot the color bar? 
                        #opacity=alpha,
                        scalar_bar_args={'color':'k',# 'interactive':True,
                                         'vertical':False,
                                         'title_font_size':16,
                                         'label_font_size':14})
                                 
        try:
            pvelec = pv.PolyData(elec)
            ax.add_mesh(pvelec, color=elec_color, point_size=10.,
                        render_points_as_spheres=True)
        except AttributeError as e:
            print("Could not plot 3d electrodes, error = "+str(e))
        
        ax.show()

    
    def _showPseudoSectionIP(self, ax=None, contour=False, vmin=None, vmax=None): #IP pseudo section
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
            
        Returns
        -------
        fig : matplotlib figure
            If `ax` is not specified, the method returns a figure.
        """
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        elecpos = self.elec['x'].values.copy()
        
        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
        
        if self.protocolIPFlag == True:
            ip = self.df['ip'].values
        else:
            ip = -self.kFactor*self.df['ip'].values

        label = r'$\phi$ [mrad]'
        

        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
        
        elecpos[self.elec['remote']] = np.inf
            
        cadd = np.abs(elecpos[array[:,0]]-elecpos[array[:,1]])/2
        cadd[np.isinf(cadd)] = 0 # they are inf because of our remote
        cmiddle = np.min([elecpos[array[:,0]], elecpos[array[:,1]]], axis=0) + cadd
        
        padd = np.abs(elecpos[array[:,2]]-elecpos[array[:,3]])/2
        padd[np.isinf(padd)] = 0
        pmiddle = np.min([elecpos[array[:,2]], elecpos[array[:,3]]], axis=0) + padd
        
        xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
        ypos = np.sqrt(2)/2*np.abs(cmiddle-pmiddle)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        
        if contour is False:
            plotPsIP = ax.scatter(xpos, ypos, c=ip, s=70, vmin=vmin, vmax=vmax)#, norm=mpl.colors.LogNorm())
#            divider = make_axes_locatable(ax)
#            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(plotPsIP, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(label)
#            ax.set_title('Phase Shift\npseudo section')
        
        else:
            if vmin is None:
                vmin = np.min(ip)
            if vmax is None:
                vmax = np.max(ip)
            levels = np.linspace(vmin, vmax, 13)
            plotPsIP = ax.tricontourf(xpos, ypos, ip, levels = levels, extend = 'both')
            fig.colorbar(plotPsIP, ax=ax, fraction=0.046, pad=0.04, label=label)
#            cbar.set_label(label)
#            ax.set_title('Phase Shift\npseudo section')
        ax.invert_yaxis() # to remove negative sign in y axis
        ax.set_title('Phase Shift\npseudo section')  
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Pseudo depth [m]')
        if ax is None:
            return fig
    
    
    def write2protocol(self, outputname='', err=False, errTot=False,
                       ip=False, res0=False, isubset=None, threed=False):
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
                    print('Using total error')
                    if 'modErr' not in df.columns:
                        raise ValueError('ERROR : you must specify a modelling error')
                    else: # if present, compute geometric mean of the errors
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
                    protocol.loc[:, c] = '1 ' + protocol[c]
            # protocol.insert(1, 'sa', 1)
            # protocol.insert(3, 'sb', 1)
            # protocol.insert(5, 'sm', 1)
            # protocol.insert(7, 'sn', 1)
        
        # write protocol.dat
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(len(protocol)) + '\n')
            with open(outputname, 'a') as f:
                protocol.to_csv(f, sep='\t', header=False, index=False)
        
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
                     label=None, vmin=None, vmax=None, elec=True):
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
            If `True``then all data will be log transformed.
        label : str, optional
            Label of the colorbar. If None, will be given from the 'attr' argument.
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        elec : bool, optional
            If `True`, the electrodes are shown and can be used for filtering.
        """
        lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        array = self.df[['a','b','m','n']].replace(lookupDict).values
        if len(array) == 0:
            raise ValueError('Unable to plot! Dataset is empty - can be due to filtering out all datapoints')
        
        percFact = 1
        if label is None:
            if attr == 'app':
                self.computeK()
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

        self.eselect = np.zeros(len(elecpos), dtype=bool)
        
        if log:
            val = np.sign(val)*np.log10(np.abs(val))

        array = np.sort(array, axis=1) # need to sort the array to make good wenner pseudo section
        
        cadd = np.abs(elecpos[array[:,0]]-elecpos[array[:,1]])/2
        cadd[np.isinf(cadd)] = 0 # they are nan because of our remote
        cmiddle = np.min([elecpos[array[:,0]], elecpos[array[:,1]]], axis=0) + cadd
        
        padd = np.abs(elecpos[array[:,2]]-elecpos[array[:,3]])/2
        padd[np.isinf(padd)] = 0
        pmiddle = np.min([elecpos[array[:,2]], elecpos[array[:,3]]], axis=0) + padd
        
        xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
        ypos = np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
        
        def onpick(event):
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
        
        if elec:
            caxElec, = ax2.plot(elecpos, np.zeros(len(elecpos)), 'ko', picker=5)
        else:
            caxElec, = ax2.plot([], [], 'ko')
        
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
        for e in elec:
            i2keep = (self.df[['a','b','m','n']].values != e).all(1)
            self.filterData(i2keep)
        if debug:
            numRemoved = np.sum(~i2keep)
            msgDump = '{:d} measurements removed!'.format(numRemoved)
            print(msgDump)
            return numRemoved
    
    
    def filterAppResist(self, vmin=None, vmax=None, debug=True):
        """Filter measurements by apparent resistivity for surface surveys 
        Parameters
        -----------
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        debug : bool, optional
            Print output to screen. Default is True.
        """
        df = self.df.copy()
        self.computeK()
        appRes = self.df['K']*self.df['resist']
        if vmin is None:
            vmin = np.min(appRes)
        if vmax is None:
            vmax = np.max(appRes)
        ikeep = (appRes >= vmin) & (appRes <= vmax)
        self.df = df[ikeep]
        self.df.reset_index()
        
        if debug:
            numRemoved = len(df)-len(self.df)
            msgDump = "%i measurements outside [%s,%s] removed!" % (numRemoved, vmin, vmax)
            print(msgDump)
            return numRemoved
    
    def filterTransferRes(self, vmin=None, vmax=None, debug=True):
        """Filter measurements by transfer resistance. 
        Parameters
        -----------
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        debug : bool, optional
            Print output to screen. Default is True.
        """
        df = self.df.copy()
        tRes = self.df['resist']
        if vmin is None:
            vmin = np.min(tRes)
        if vmax is None:
            vmax = np.max(tRes)
        ikeep = (tRes >= vmin) & (tRes <= vmax)
        self.df = df[ikeep]
        self.df.reset_index()
        
        if debug:
            numRemoved = len(df)-len(self.df)
            msgDump = "%i measurements outside [%s,%s] removed!" % (numRemoved, vmin, vmax)
            print(msgDump)
            return numRemoved


    # def shuntIndexes(self, debug=True): # pragma: no cover
    #     """Normalise the indexes the sequence matrix to start at 1.
        
    #     Parameters
    #     -------------
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
         
        df = self.df
        sch_mat = np.array((df['a'],df['b'],df['m'],df['n'])).T
        idx = np.argwhere(sch_mat == old_indx)
        sch_mat[idx[:,0],idx[:,1]] = new_indx
        corrected = sch_mat
        df['a'] = corrected[:,0]
        df['b'] = corrected[:,1]
        df['m'] = corrected[:,2]
        df['n'] = corrected[:,3]
        self.df = df 
        
    
    def normElecIdx(self, debug=True): # pragma: no cover
        """Normalise the electrode indexing sequencing to start at 1 and ascend
        consectively (ie 1 , 2 , 3 , 4 ... )
        
        Function firstly normalises all indexes so that the lowest electrode 
        number is 1. Then removes jumps in the electrode indexing.
        
        Parameters
        -----------
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
        -----------
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
    
    def elec2distance(self): # pragma: no cover
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
            new_elec[put_back,2] =  z_sorted[i]
    
        self.elec.loc[:, ['x','y','z']] = new_elec
        
        
    def addPerError(self,pnct=2.5): # pragma: no cover
        """Add a flat percentage error to resistivity data.
        
        Parameters
        --------------
        pnct: float
            Error in percent
        """
        res = np.array(self.df['resist'])
        error = (pnct/100)*res
        self.df['resError'] = error + np.array(self.df['resError'])
        
        
    def estimateError(self, a_wgt=0.01, b_wgt=0.02): # pragma: no cover
        """Estimate reciprocal error data for data with no reciprocals, following
        the same routine present in R2. This allows for the additional inclusion
        of modelling errors. 
        
        Parameters
        ------------
        a_wgt: float, optional
            a_wgt documented in the R2 documentation 
        b_wgt: float, optional 
            b_wgt documented in the R2 documentation  
        """
        res = np.array(self.df['resist'])
        var_res = (a_wgt*a_wgt)+(b_wgt*b_wgt) * (res*res)
        std_res = np.sqrt(var_res)
        self.df['resError'] = std_res


    def exportSrv(self, fname=None): # pragma: no cover
        """Export .srv format for which is compatible with E4D. The e4d survey
        file includes the electrode locations, in addition to the scheduling 
        matrix. 
        
        Paramters
        ------------
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
        fh.write('%i number of electrodes\n'%numelec)
        for i in range(numelec):
            line = '{:d} {:f} {:f} {:f} {:d}\n'.format(i+1,
                    elec[i,0],#x coordinate
                    elec[i,1],#y coordinate
                    elec[i,2],#z coordinate
                    1)#buried flag 
            fh.write(line)
        #now write the scheduling matrix to file 
        ie = self.df['irecip'].values >= 0 # reciprocal + non-paired
        df = self.df[ie]
        nomeas = len(df) # number of measurements 
        df = df.reset_index().copy()
        fh.write('\n%i number of measurements \n'%nomeas)
        if not 'resError' in self.df.columns: # the columns exists
            self.estError()
        # format >>> m_indx a b m n V/I stdev_V/I
        for i in range(nomeas): 
            line = '{:d} {:d} {:d} {:d} {:d} {:f} {:f}\n'.format(i+1,
                    df['a'][i],
                    df['b'][i],
                    df['m'][i],
                    df['n'][i],
                    df['recipMean'][i],
                    df['resError'][i])
            fh.write(line)

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
        
        
        