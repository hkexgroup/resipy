#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:21:54 2018

@author: jkl
"""
#import sys
import sys, os, platform
#sys.path.append(os.path.relpath('../resipy'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas as pd
from scipy.stats import norm
from scipy.stats.kde import gaussian_kde
#import statsmodels.formula.api as smf

from resipy.parsers import (syscalParser, protocolParser,protocolParserLME,  resInvParser,
                     primeParser, primeParserTab, protocolParserIP,
                     protocol3DParser, forwardProtocolDC, forwardProtocolIP,
                     stingParser, ericParser, lippmannParser, aresParser)
from resipy.DCA import DCA

import warnings
warnings.simplefilter('default', category=DeprecationWarning) # this will show the deprecation warnings

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
    def __init__(self, fname, ftype='', name='', spacing=None, parser=None, keepAll=True):
        self.elec = []
        self.df = pd.DataFrame()
        if name == '':
            name = os.path.basename(os.path.splitext(fname)[0])
        self.name = name
        self.iBorehole = False # True is it's a borehole
        self.protocolIPFlag = False
        self.kFactor = 1
        self.errorModel = None # function instanticated after fitting an error model with reciprocal errors
        self.iselect = None # use in filterManual()
        self.eselect = None # idem
        self.iremote = None # to be set by R2 class when remote detected
        
        avail_ftypes = ['Syscal','Protocol','Res2Dinv', 'BGS Prime', 'ProtocolIP',
                        'Sting', 'ABEM-Lund', 'Lippmann', 'ARES']# add parser types here! 
        
        if parser is not None:
            elec, data = parser(fname)
        else:
            if ftype == 'Syscal':
                elec, data = syscalParser(fname, spacing=spacing)
                self.kFactor = 1.2
            elif ftype =='Protocol':
                elec, data = protocol3DParser(fname)
            elif ftype == 'Res2Dinv':
                elec, data = resInvParser(fname)
            elif ftype == 'BGS Prime':
#                try:
#                    elec, data = primeParser(fname)
#                except:
                elec, data = primeParserTab(fname)
            elif ftype == 'ProtocolIP':
                elec, data = protocolParserIP(fname)
                self.protocolIPFlag = True
            elif ftype == 'forwardProtocolDC':
                elec, data = forwardProtocolDC(fname)
            elif ftype == 'Sting':
                elec, data = stingParser(fname)
            elif ftype == 'ABEM-Lund':
                elec, data = ericParser(fname)
            elif ftype == 'Lippmann':
                elec, data = lippmannParser(fname)
            elif ftype == 'ARES':
                elec ,data = aresParser(fname)
#            elif ftype == 'forwardProtocolIP':
#                self.protocolIPFlag = True
#                elec, data = forwardProtocolIP(fname)
    #        elif (ftype == '') & (fname == '') & (elec is not None) and (data is not None):
    #            pass # manual set up
    #            print('Manual set up, no data will be imported')
            else:
                print("Unrecognised ftype, available types are :",avail_ftypes )
                raise Exception('Sorry this file type is not implemented yet')

        
        self.df = data
        for c in ['a','b','m','n']:
            self.df.loc[:,c] = self.df[c].astype(int)
        
        # add error measured to the error columns (so they can be used if no error model are fitted)
        if 'magErr' in self.df.columns:
            self.df['resError'] = self.df['magErr'].copy()
        else:
            self.df['resError'] = np.nan
        if 'phiErr' in self.df.columns:
            self.df['phaseError'] = self.df['phiErr'].copy()
        else:
            self.df['phaseError'] = np.nan
            
        self.dfOrigin = data.copy() # unmodified
        self.elec = elec
        self.ndata = len(data)
        self.phiCbarmin = 0
        self.phiCbarMax = 25
        self.filt_typ = None
        self.cbar = True
        self.filterDataIP = pd.DataFrame()

        if ftype == 'BGS Prime':
            self.checkTxSign()
        if ftype == 'Syscal':
            if np.all(self.df['vp'] >= 0):
                self.checkTxSign()
            

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
        surrogate_file = '\examples\dc-2d\protocol.dat'
        path2resipy = os.path.realpath(__file__).replace('\\src\resipy\\',surrogate_file)
        path = path2resipy # map to the protocol test file
        s_svy = cls(fname = path, ftype='Protocol')#surrogate survey
        s_svy.df = df
        s_svy.elec = elec
        s_svy.ndata = 1#len(data)
        s_svy.computeReciprocal()
        s_svy.filterDefault()
        return s_svy
    
    
    def __str__(self):
        out = "Survey class with %i measurements and %i electrodes"%(len(self.df),len(self.elec[:,0]))
        return out
 
    
    def calcGeomFactor(self):
        """Calculating geometric factor (K) for the data
        Returns
        -------
        K : array
            Array of geometric factors for each quadrupoles.
        """
        
        elecpos = self.elec[:,0]
        array = self.df[['a','b','m','n']].values.copy().astype(int)
        
        apos = elecpos[array[:,0]-1]
        bpos = elecpos[array[:,1]-1]
        mpos = elecpos[array[:,2]-1]
        npos = elecpos[array[:,3]-1]
        AM = np.abs(apos-mpos)
        BM = np.abs(bpos-mpos)
        AN = np.abs(apos-npos)
        BN = np.abs(bpos-npos)
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        return K
        
    
    def checkTxSign(self):
        """Checking the sign of the transfer resistances (2D survey only !).
        """
        elecpos = self.elec[:,0]
        array = self.df[['a','b','m','n']].values.copy().astype(int)
        resist = self.df['resist'].values.copy()
        
        apos = elecpos[array[:,0]-1]
        bpos = elecpos[array[:,1]-1]
        mpos = elecpos[array[:,2]-1]
        npos = elecpos[array[:,3]-1]
        AM = np.abs(apos-mpos)
        BM = np.abs(bpos-mpos)
        AN = np.abs(apos-npos)
        BN = np.abs(bpos-npos)
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        
        ie = ((K < 0) & (resist > 0)) | ((K > 0) & (resist < 0))
        self.df.loc[ie, 'resist'] = resist[ie]*-1
        print('WARNING: change sign of ', np.sum(ie), ' Tx resistance.')
        

    def filterDefault(self):
        """Remove NaN, Inf and duplicates values in the data frame.
        """
        # remove Inf and NaN
        resist = self.df['resist'].values
        iout = np.isnan(resist) | np.isinf(resist)
        if np.sum(iout) > 0:
            print('Survey.filterDefault: Number of Inf or NaN : ', np.sum(iout))
        print('Inf or NaN: ', end='')
        self.filterData(~iout)
        
        # remove duplicates
        shapeBefore = self.df.shape[0]
        self.df = self.df.drop_duplicates(subset=['a','b','m','n'], keep = 'first')
        ndup = shapeBefore - self.df.shape[0]
        if ndup > 0:
            print('Survey.filterDefault: ', ndup, 'duplicates removed.')
        
        # remove quadrupoles were A or B are also potential electrodes
        ie1 = self.df['a'].values == self.df['m'].values
        ie2 = self.df['a'].values == self.df['n'].values
        ie3 = self.df['b'].values == self.df['m'].values
        ie4 = self.df['b'].values == self.df['n'].values
        ie = ie1 | ie2 | ie3 | ie4
        if np.sum(ie) > 0:
            print('Survey.filterDefault: ', np.sum(ie), 'measurements with A or B == M or N')
        print('strange quadrupoles: ', end='')
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
          
        
    def addData(self, fname, ftype='Syscal', spacing=None, parser=None):
        """Add data to the actual survey (for instance the reciprocal if they
        are not in the same file).
        """
         
        if parser is not None:
            elec, data = parser(fname)
        else:
            if ftype == 'Syscal':
                elec, data = syscalParser(fname, spacing=spacing)
                self.kFactor = 1.2
            elif ftype =='Protocol':
                elec, data = protocol3DParser(fname)
            elif ftype == 'Res2Dinv':
                elec, data = resInvParser(fname)
            elif ftype == 'BGS Prime':
                try:
                    elec, data = primeParser(fname)
                except:
                    elec, data = primeParserTab(fname)
            elif ftype == 'ProtocolIP':
                elec, data = protocolParserIP(fname)
                self.protocolIPFlag = True
            elif ftype == 'Sting':
                elec, data = stingParser(fname)
            elif ftype == 'ABEM-Lund':
                elec, data = ericParser(fname)
            elif ftype == 'Lippmann':
                elec, data = lippmannParser(fname)
            elif ftype == 'ARES':
                elec ,data = aresParser(fname)
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
        
        
    def inferType(self):
        """Define the type of the survey.
        """
        if self.elec[:,2].sum() == 0:
            stype = '2d'
        else:
            stype = 'undefined'
        
        # check if borehole
        # check if IP survey
        return stype
    
    
    def setType(self, stype=None):
        """Set the type of the survey.
        
        Parameters
        ----------
        stype : str
            Type of the survey. Could also be infered using 
            `Survey.inferType()`.
        """
        if stype == None:
            self.stype = self.inferType()
        else:
            self.stype = stype


    def computeReciprocal(self):
        """Compute reciprocal measurements.
        
        Notes
        -----
        The methods create an array where all positive index are normal
        measurements and the reciprocal measurements are their negative
        counterparts. This array is stored in the main dataframe `Survey.df`
        in the columns `irecip`. Measurements with `ìrecip=0` are measurements
        without reciprocal.
        """
        resist = self.df['resist'].values
        phase = -self.kFactor*self.df['ip'].values #converting chargeability to phase shift
        ndata = self.ndata
        array = self.df[['a','b','m','n']].values
        
        R = np.copy(resist)
        M = np.copy(phase)
        ndata = len(R)
        Ri = np.zeros(ndata)
        reciprocalErr = np.zeros(ndata)*np.nan
        reciprocalErrRel = np.zeros(ndata)*np.nan
        reciprocalMean = np.zeros(ndata)*np.nan
        reci_IP_err = np.zeros(ndata)*np.nan
        # search for reciprocal measurement
        count=1
        notfound=0
        for i in range(0,ndata):
            rev1=[2,3,0,1]
            rev2=[3,2,0,1]
            rev3=[2,3,1,0]
            rev4=[3,2,1,0]
            index1=(array[:,:] == array[i,rev1]).all(1)
            index2=(array[:,:] == array[i,rev2]).all(1)
            index3=(array[:,:] == array[i,rev3]).all(1)
            index4=(array[:,:] == array[i,rev4]).all(1)
            index=index1|index2|index3|index4
            
            if len(index[index]) == 1:
                reciprocalErr[index] = np.abs(R[i])-np.abs(R[index])
                reciprocalErrRel[index] = (np.abs(R[i])-np.abs(R[index]))/np.abs(R[i]) # in percent
                # flag first reciprocal found like this we can
                # delete the second one when we find it
                # (no loss of information)
                reci_IP_err[index] = M[i]-M[index]
                if Ri[i] == 0: # only if Ri == 0 otherwise we will
                    # overwrite all the data
                    Ri[i] = count # flag the first found
                if Ri[index] == 0:
                    Ri[index] = -count # flag its reciprocal
                # replace reciprocalMean by one measurements if the other one
                # is bad (NaN or Inf). Hopefully, the error model will find an
                # error to go with
                ok1 = ~(np.isnan(R[i]) | np.isinf(R[i]))
                ok2 = ~(np.isnan(R[index]) | np.isinf(R[index]))
                if ok1 & ok2:
                    reciprocalMean[i]=np.mean([np.abs(R[i]),np.abs(R[index])])
                elif ok1 & ~ok2:
                    reciprocalMean[i] = np.abs(R[i])
                elif ~ok1 & ok2:
                    reciprocalMean[i] = np.abs(R[index])
            else:
                #print("no reciprocal found for "+str(array[i,:]))
                notfound=notfound+1
            count=count+1
        print(str(notfound)+'/'+str(ndata)+' reciprocal measurements NOT found.')
        reciprocalMean = np.sign(resist)*reciprocalMean # add sign
#        ibad = np.abs(reciprocalErrRel) > 0.2 # error if np.nan in this
        ibad = np.array([np.abs(a) > 0.2 if ~np.isnan(a) else False for a in reciprocalErrRel])
        print(str(np.sum(ibad)) + ' measurements error > 20 %')
        
        irecip = Ri        
        
        self.df['irecip'] = irecip
        self.df['reciprocalErrRel'] = reciprocalErrRel
        self.df['recipError'] = reciprocalErr
        self.df['recipMean'] = reciprocalMean
        self.df['reci_IP_err'] = reci_IP_err
        # in order to compute error model based on a few reciprocal measurements
        # we fill 'recipMean' column with simple resist measurements for lonely
        # quadrupoles (which do not have reciprocals)
        inotRecip = irecip == 0
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
        i2keep = np.abs(self.df['a'] - self.df['b']) == np.abs(self.df['m'] - self.df['n'])
        self.filterData(i2keep)
        
    
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
        
        
    def addFilteredIP(self):
        """Add filtered IP data after IP filtering and pre-processing. This is
        because the IP filtering is done on a different dataframe and only
        merged when called this method.
        """
        self.df = pd.merge(self.df, self.filterDataIP[['a','b','m','n']].copy(), how='inner', on=['a','b','m','n'])

    
    @staticmethod
    def logClasses3(datax, datay, func, class1=None):
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
        reciprocalMean = self.df['recipMean'].values
        reciprocalErr = self.df['recipError'].values
        ax.loglog(np.abs(reciprocalMean), reciprocalErr, 'o')
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
    
    @staticmethod
    def sign_coef(x):
        if x>=0:
            t = '+'
        else:
            t = '-'
        return t    
    
    
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
#        coefs= np.linalg.lstsq(np.vstack([bins[:,0], np.ones(len(bins[:,0]))]).T, bins[:,1], rcond=None)[0] # calculating fitting coefficients (a,m) 
        coefs = np.polyfit(bins[:,0], bins[:,1], 1)
#        if coefs[1] < 0: # we don't want negative error -> doesn't make sense
#            slope = np.polyfit(bins[:,0], bins[:,1], 0)
#            coefs = [slope, 0]
        R_error_predict = ((coefs[0])*(bins[:,0]))+coefs[1] # error prediction based of linear model        
        print(np.min(R_error_predict)) # TODO negative error here ! that's why the red fit line goes down
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
        
    
    def fitErrorLME(self, iplot=True, ax=None, rpath=None):
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
            ax.set_xlabel('Reciprocal Error Observed [$\Omega$]')
            ax.set_ylabel('Reciprocal Error Predicted [$\Omega$]')
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
        heat_recip_Filter = temp_heatmap_recip_filterN.set_index(['m','a']).Phase.unstack(0)     
        if ax is None:
            fig, ax = plt.subplots()  
        else:
            fig = ax.get_figure()             
        m = ax.imshow(heat_recip_Filter, origin='lower',cmap='jet',vmin=self.phiCbarmin, vmax=self.phiCbarMax)
        ax.xaxis.set_ticks(np.arange(0,filterDataIP_plotOrig['a'].max()+1,4))
        ax.yaxis.set_ticks(np.arange(0,filterDataIP_plotOrig['m'].max(),4))
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

    
    def showPseudo(self, ax=None, bx=None, **kwargs):
        """Plot pseudo section if 2D survey or just quadrupoles transfer
        resistance otherwise.
        """
        if bx is None:
            bx = self.iBorehole
        if bx is False:
            self._showPseudoSection(ax=ax, **kwargs)
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
            ax.set_ylabel('Chargeability [mV/V]')
    
    
#    def computeK(self):
#        """Compute geomatrix factor (assuming flat 2D surface) and store it
#        in self.df['K'].
#        """
#        array = self.df[['a','b','m','n']].values.astype(int)
#        elecpos = self.elec[:,0]
#
#        apos = elecpos[array[:,0]-1]
#        bpos = elecpos[array[:,1]-1]
#        mpos = elecpos[array[:,2]-1]
#        npos = elecpos[array[:,3]-1]
#        AM = np.abs(apos-mpos)
#        BM = np.abs(bpos-mpos)
#        AN = np.abs(apos-npos)
#        BN = np.abs(bpos-npos)
#        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
#        
#        self.df['K'] = K
        
    def computeK(self):
        """Compute geomatrix factor (assuming flat 2D surface) and store it
        in self.df['K'].
        """
        array = self.df[['a','b','m','n']].values.astype(int)
        elec = self.elec
        aposx = elec[:,0][array[:,0]-1]
        aposy = elec[:,1][array[:,0]-1]
        aposz = elec[:,2][array[:,0]-1]
        
        bposx = elec[:,0][array[:,1]-1]
        bposy = elec[:,1][array[:,1]-1]
        bposz = elec[:,2][array[:,1]-1]
        
        mposx = elec[:,0][array[:,2]-1]
        mposy = elec[:,1][array[:,2]-1]
        mposz = elec[:,2][array[:,2]-1]
        
        nposx = elec[:,0][array[:,3]-1]
        nposy = elec[:,1][array[:,3]-1]
        nposz = elec[:,2][array[:,3]-1]
        
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
        array = self.df[['a','b','m','n']].values.astype(int)
        elecpos = self.elec[:,0].copy()
        resist = self.df[column].values
        
        if geom: # compute and applied geometric factor
#            apos = elecpos[array[:,0]-1]
#            bpos = elecpos[array[:,1]-1]
#            mpos = elecpos[array[:,2]-1]
#            npos = elecpos[array[:,3]-1]
#            AM = np.abs(apos-mpos)
#            BM = np.abs(bpos-mpos)
#            AN = np.abs(apos-npos)
#            BN = np.abs(bpos-npos)
#            K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
            self.computeK()
            resist = resist*self.df['K']

        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
            
        if log:
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
        else:
            label = r'$\rho_a$ [$\Omega.m$]'
                
        # cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
        #     + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        # pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
        #     + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        
        if self.iremote is not None:
            elecpos[self.iremote] = np.inf # so it will never be taken as minimium
        
        cadd = np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        cadd[np.isinf(cadd)] = 0 # they are inf because of our remote
        cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) + cadd
        
        padd = np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        padd[np.isinf(padd)] = 0
        pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) + padd

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
        array = self.df[['a','b','m','n']].values.astype(int)
        elecpos = self.elec[:,0].copy()
        
        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
        
        if self.protocolIPFlag == True:
            ip = self.df['ip'].values
        else:
            ip = -self.kFactor*self.df['ip'].values

        label = r'$\phi$ [mrad]'
        

        # sorting the array in case of Wenner measurements (just for plotting)
        array = np.sort(array, axis=1) # for better presentation
        # cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
        #     + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        # pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
        #     + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        
        if self.iremote is not None:
            elecpos[self.iremote] = np.inf # so it will never be taken as minimium
        
        cadd = np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        cadd[np.isinf(cadd)] = 0 # they are inf because of our remote
        cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) + cadd
        
        padd = np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        padd[np.isinf(padd)] = 0
        pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) + padd
        
        xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
        ypos = np.sqrt(2)/2*np.abs(cmiddle-pmiddle)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        if contour is False:
#            if ax is None:
#                fig, ax = plt.subplots()
#            else:
#                fig = ax.get_figure()
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
            protocol.insert(1, 'sa', 1)
            protocol.insert(3, 'sb', 1)
            protocol.insert(5, 'sm', 1)
            protocol.insert(7, 'sn', 1)
        
        # write protocol.dat
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(len(protocol)) + '\n')
            with open(outputname, 'a') as f:
                protocol.to_csv(f, sep='\t', header=False, index=False)
        
        return protocol
        
        
    def filterDCA(self, dump=print):
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
        dump(100)
        
    
    def filterManual(self, attr='resist', ax=None, log=False, geom=True,
                     label=r'Apprent Resistivity [$\Omega.m$]',
                     vmin=None, vmax=None):
        """Manually filters the data visually. The points manually selected are
        flagged in the `Survey.iselect` vector and can subsequently be removed
        by calling `Survey.filterData(~Survey.iselect)`.
        
        Parameters
        ----------
        attr : str, optional
            Columns of `Survey.df` to use for plotting. Default is `resist` (
            transfer resistance).
        ax : matplotlib axis, optional
            If specified, the graph is plotted along the axis.
        log : bool, optional
            If `True``then all data will be log transformed.
        geom : bool, optional
            If `True` the values will be multiplied by the geometric factor (default).
        label : str, optional
            Label of the colorbar.
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        """
        array = self.df[['a','b','m','n']].values.astype(int)
        if len(array) == 0:
            raise ValueError('Unable to plot! Dataset is empty - can be due to filtering out all datapoints')

        resist = self.df[attr].values
        inan = np.isnan(resist)
        resist = resist.copy()[~inan]
        array = array.copy()[~inan]
        self.iselect = np.zeros(len(inan), dtype=bool)
        
        def setSelect(ie, boolVal):
            ipoints[ie] = boolVal
            self.iselect[~inan] = ipoints
        if self.iremote is not None:
            spacing = np.mean(np.diff(self.elec[~self.iremote,0]))
        else:
            spacing = np.mean(np.diff(self.elec[:,0]))
        nelec = np.max(array)
        elecpos = np.arange(0, spacing*nelec, spacing)
        
        self.eselect = np.zeros(len(elecpos), dtype=bool)
        
        if geom: # compute and applied geometric factor
            # apos = elecpos[array[:,0]-1]
            # bpos = elecpos[array[:,1]-1]
            # mpos = elecpos[array[:,2]-1]
            # npos = elecpos[array[:,3]-1]
            # AM = np.abs(apos-mpos)
            # BM = np.abs(bpos-mpos)
            # AN = np.abs(apos-npos)
            # BN = np.abs(bpos-npos)
            # K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
            self.computeK()
            resist = resist*self.df['K']
            
        if log:
            resist = np.sign(resist)*np.log10(np.abs(resist))

        array = np.sort(array, axis=1) # need to sort the array to make good wenner pseudo section
        # cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
        #     + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        # pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
        #     + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        
        if self.iremote is not None:
            elecpos[self.iremote] = np.inf # so it will never be taken as minimium
        
        cadd = np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        cadd[np.isinf(cadd)] = 0 # they are inf because of our remote
        cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) + cadd
        
        padd = np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        padd[np.isinf(padd)] = 0
        pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) + padd
        
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
        ax.set_xlabel('Electrode number')
        ax.set_ylabel('Pseudo depth [m]')
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_ticks_position('top')
        ax.set_xticks(elecpos[::5])
        ax.set_xticklabels(elecNumber[::5])
        
        
        caxElec, = ax2.plot(elecpos, np.zeros(len(elecpos)), 'ko', picker=5)
        
        cax = ax2.scatter(xpos, ypos, c=resist, marker='o', picker=5, vmin=vmin,
                         vmax=vmax)
        cbar = fig.colorbar(cax, ax=ax2)
        cbar.set_label(label)
        ax2.set_xlabel('Distance [m]')
        ax2.xaxis.set_label_position('bottom')
        ax2.xaxis.set_ticks_position('bottom')
        cax.figure.canvas.mpl_connect('pick_event', onpick)        
        killed, = cax.axes.plot([],[],'rx')
        elecKilled, = cax.axes.plot([],[],'rx')
        x = cax.get_offsets()[:,0]
        y = cax.get_offsets()[:,1]        
        ipoints = np.zeros(len(y),dtype=bool)
        lines = {cax:'data',caxElec:'elec',killed:'killed'}
        
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


    
    def filterElec(self, elec=[]):
        """Filter out specific electrodes given.
        
        Parameters
        ----------
        elec : list
            List of electrode number to be removed.
        
        """
        for e in elec:
            i2keep = (self.df[['a','b','m','n']].values != e).all(1)
            self.filterData(i2keep)
            print(np.sum(~i2keep), '/', len(i2keep), 'quadrupoles removed.')
    
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


    def shuntIndexes(self, debug=True): 
        """Normalise the indexes the sequence matrix to start at 1.
        
        Parameters
        -------------
        debug : bool, optional
            Set to True to print output.
        """
        df = self.df
        sch_mat = np.array((df['a'],df['b'],df['m'],df['n'])).T
        imin = np.min(sch_mat)
        if imin != 1:
            if debug:
                print("It looks like scheduling indexing starts at: %i"%imin)
                print("...normalising electrode indexing to start at 1.")
            #check to see if negative electrode numbers occur. 
            if imin < 0:
                #check there is a electrode at 0
                imin_pos = np.min(np.abs(sch_mat))
                if imin_pos == 1:
                    if debug:
                        print("Any positive electrode indexes now start at %i"%(abs(imin)+1))
                    crr_idx = np.argwhere(sch_mat>0)
                    sch_mat[crr_idx] -= 1
                    
            corrected = sch_mat - (imin - 1)
            #return corrected
            df['a'] = corrected[:,0]
            df['b'] = corrected[:,1]
            df['m'] = corrected[:,2]
            df['n'] = corrected[:,3]
            self.df = df 


    def swapIndexes(self, old_indx, new_indx):
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
        
    
    def normElecIdx(self, debug=True):
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
                
    def normElecIdxwSeq(self,expected=None):
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
    
    def elec2distance(self):
        """Convert 3d xy data in pure x lateral distance.
        Use for 2D data only!
        """
        elec = self.elec.copy()
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
    
        self.elec = new_elec
        
    def addPerError(self,pnct=2.5):
        """Add a flat percentage error to resistivity data.
        
        Parameters
        --------------
        pnct: float
            Error in percent
        """
        res = np.array(self.df['resist'])
        error = (pnct/100)*res
        self.df['resError'] = error + np.array(self.df['resError'])
        
        
    def estimateError(self, a_wgt=0.01, b_wgt=0.02):
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

    def exportSrv(self, fname=None):
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
        fh.write('%i number of electrodes\n'%numelec)
        for i in range(numelec):
            line = '{:d} {:f} {:f} {:f} {:d}\n'.format(i+1,
                    self.elec[i,0],#x coordinate
                    self.elec[i,1],#y coordinate
                    self.elec[i,2],#z coordinate
                    1)#buried flag 
            fh.write(line)
        #now write the scheduling matrix to file 
        nomeas = len(self.df) # number of measurements 
        df = self.df.reset_index().copy()
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
                    df['resist'][i],
                    df['resError'][i])
            fh.write(line)


#%% deprecated methods
    def basicFilter(self):
        warnings.warn('This function is deprecated, use filterDefault() instead',
                      DeprecationWarning)
        self.filterDefault()
        
    def removeUnpaired(self):
        warnings.warn('This function is deprecated, use filterUnpaired() instead',
              DeprecationWarning)
        n = self.filterUnpaired()
        return n
    
    def estError(self, a_wgt=0.01, b_wgt=0.02):
        warnings.warn('The function is deprecated, use estimateError() instead.',
                      DeprecationWarning)
        self.estimateError(a_wgt=a_wgt, b_wgt=b_wgt)
            
    
    def filterdip(self, elec): # deleted specific elec data
        warnings.warn('The function is deprecated, use filterElec() instead.',
                      DeprecationWarning)
        index = (self.array == elec[0]).any(-1)
        for i in range(1,len(elec)):
            index = index | (self.array == elec[i]).any(-1)
        n = self.filterData(~index)
        return n


    def dca(self, dump=print):
        warnings.warn('The function is deprecated, use filterDCA() instead.',
                      DeprecationWarning)
        self.filterDCA(dump=dump)

        
    def manualFiltering(self, ax=None, figsize=(12,3), contour=False,
                        log=False, geom=False, label='', vmin=None, vmax=None):
        warnings.warn('The function is deprecated, use filterManual() instead.',
                      DeprecationWarning)
        self.filterManual(ax=ax, figsize=figsize, contour=contour,
                          log=log, geom=geom, label=label, vmin=vmin, vmax=vmax)
 
    def pseudo(self, ax=None, bx=None, **kwargs):
        warnings.warn('The function is deprecated, use showPseudo() instead.',
                      DeprecationWarning)
        self.showPseudo(ax=ax, bx=bx, **kwargs)
    
    
    def pseudoIP(self, ax=None, bx=None, **kwargs):
        warnings.warn('The function is deprecated, use showPseudoIP() instead.',
                      DeprecationWarning)
        self.showPseudoIP(ax=ax, bx=bx, **kwargs)
    
    
    def reciprocal(self):
        warnings.warn('This function is deprecated, use computeReciprocal() instead',
              DeprecationWarning)
        out = self.computeReciprocal()
        return out
    
    
    def errorDist(self, ax=None):
        warnings.warn('This function is deprecated, use showErrorDist() instead',
              DeprecationWarning)
        self.showErrorDist(ax=ax)


    def removeDummy(self):
        warnings.warn('This function is deprecated, use filterDummy() instead',
              DeprecationWarning)
        n = self.filterDummy()
        return n
    
    
    def plotError(self, ax=None):
        warnings.warn('The function is deprecated, use showError() instead.',
                      DeprecationWarning)
        self.showError(ax=ax)
    
        
    def phaseplotError(self, ax=None):
        warnings.warn('The function is deprecated, use showErrorIP() instead.',
                      DeprecationWarning)
        self.showErrorIP(self, ax=ax)
        
        
    def pwlfitIP(self, ax=None):
        warnings.warn('The function is deprecated, use fitErrorPwlIP() instead.',
                      DeprecationWarning)
        self.fitErrorPwlIP(ax=ax)
        
    
    def plotIPFitParabola(self, ax=None):
        warnings.warn('The function is deprecated, use fitErrorParabolaIP() instead.',
                      DeprecationWarning)
        self.fitErrorParabolaIP(ax=ax)
              
    
    def pwlfit(self, ax=None):
        warnings.warn('The function is deprecated, use fitErrorPwl() instead.',
                      DeprecationWarning)
        self.fitErrorPwl(ax=ax)
        

    def lmefit(self, iplot=True, ax=None, rpath=None):
        warnings.warn('The function is deprecated, use fitErrorLME() instead.',
                      DeprecationWarning)
        self.fitErrorLME(iplot=iplot, ax=ax, rpath=rpath)
    
    
    def heatmap(self, ax=None):
        warnings.warn('The function is deprecated, use showHeatmap instead.',
                      DeprecationWarning)
        self.showHeatmap()
    
    
    def iprangefilt(self, phimin, phimax):
        warnings.warn('The function is deprecated, use showError() instead.',
                      DeprecationWarning)
        self.filterRangeIP(phimin, phimax)
    
    
    def removerecip(self):
        warnings.warn('The function is deprecated, use filterRecip() instead.',
                      DeprecationWarning)
        self.filterRecip(self)
    
    
    def removenested(self):
        warnings.warn('The function is deprecated, use filterNested() instead.',
                      DeprecationWarning)
        self.filterNested()
        
    def removeneg(self):
        warnings.warn('The function is deprecated, use filterNegative() instead.',
                      DeprecationWarning)
        self.filterNegative()
        
        
        