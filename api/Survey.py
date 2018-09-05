#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 11:21:54 2018

@author: jkl
"""
#import sys
import os
#sys.path.append(os.path.relpath('../api'))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import statsmodels.formula.api as smf

from api.parsers import syscalParser, protocolParser, res2invInputParser
from api.DCA import DCA

class Survey(object):
    """ Class that handles geophysical data and some basic functions. One 
    instance is created for each survey.
    """
    def __init__(self, fname, ftype='', name='', spacing=None):
        """ Create a Survey object, that contains the data and some basic 
        procedures.
        
        Parameters
        ----------
        fname : str
            Name of the file where the data are.
        ftype : str
            Type of the data file.
        name : str
            A personal name for the survey.
        spacing : float, optional
            This will be passed to the parser function to determine the
            electrode positions.
        elec : array
            Array with 3 columsn with the electrodes position in X, Y and Z.
        data : pandas.DataFrame
            Dataframe containing the 'a','b','m','n','resist' columns. Note
            that  `resist` is the *transfer resistance* not the apparent
            resistivity.
        """
        self.elec = []
        self.df = pd.DataFrame()
        self.name = name
        
        if ftype == 'Syscal':
            elec, data = syscalParser(fname, spacing=spacing)
        elif ftype =='Protocol':
            elec, data = protocolParser(fname)
        elif ftype == 'Res2Dinv':
            elec, data = res2invInputParser(fname)
#        elif (ftype == '') & (fname == '') & (elec is not None) and (data is not None):
#            pass # manual set up
#            print('Manual set up, no data will be imported')
        else:
            raise Exception('Sorry this file type is not implemented yet')
        
        self.df = data
        self.dfphasereset = pd.DataFrame() #for preserving phase reset ability
        self.dfOrigin = data.copy() # unmodified
        self.elec = elec
        self.ndata = len(data)
        self.filt_typ = None
        self.cbar = True
        self.filterDataIP = pd.DataFrame()
#        self.filterDataIP_plotOrig = data[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first').copy()
#        self.typ = 'R2'
#        self.errTyp = 'obs'
#        self.errTypIP = 'none'
        
        # to discard for the dataframe structure
#        self.array = self.df[['a','b','m','n']].values
#        self.resist = self.df['resist'].values
#        self.dev = self.df['dev'].values
        
        irecip = self.reciprocal()
#        self.mask = np.ones(self.ndata, dtype=bool) # mask of used data
        
        if all(irecip == 0) == False: # contains reciprocal
            self.basicFilter()
        
#        self.typ = 'R2' # or cR2 or R3, cR3
#        self.errTyp = 'none' # type of error to add for DC
#        self.errTypIP = 'none' # type of error to add for IP phase
 
            
    
    def basicFilter(self):
        """ Remove NaN and Inf values in the data.
        """
        resist = self.df['resist'].values
        irecip = self.df['irecip'].values
        iout = np.isnan(resist) | np.isinf(resist)
        if np.sum(iout) > 0:
            print('BAD transfer resistance data : ', np.sum(iout))
        self.filterData(~iout)
        self.filterData(irecip != 0) # filter out dummy and non reciprocal
        self.dfphasereset = self.df.copy()
    
    
    def addData(self, fname, ftype='Syscal'):
        """ Add data to the actual survey (for instance the reciprocal if they
        are not in the same file).
        """
        if ftype == 'Syscal':
            elec, data = syscalParser(fname)
        else:
            raise Exception('Sorry this file type is not implemented yet')
        self.df = self.df.append(data)
        self.ndata = len(self.df)
        self.reciprocal()
        self.basicFilter() # we assume the user input reciprocal data not another
        # normal survey
        
    
    def filterData(self, i2keep):
        """ Filter out the data not retained in `i2keep`.
        
        Parameters
        ----------
        i2keep : ndarray of bool
            Index where all measurement to be retained are `True` and the
            others `False`.
        """
        self.ndata = len(i2keep)
        self.df = self.df[i2keep]
    
    
    def inferType(self):
        """ define the type of the survey
        """
        if self.elec[:,2].sum() == 0:
            stype = '2d'
        else:
            stype = 'undefined'
        
        # check if borehole
        # check if IP survey
        return stype
    
    
    def setType(self, stype=None):
        """ Set the type of the survey.
        
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

        
    def reciprocal(self):
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
        phase = -1.2*self.df['ip'].values #converting chargeability to phase shift
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
        ibad = reciprocalErrRel > 0.2
        print(str(np.sum(ibad)) + ' measurements error > 20 %')
        
        irecip = Ri        
        
        self.df['irecip'] = irecip
        self.df['reciprocalErrRel'] = reciprocalErrRel
        self.df['recipError'] = reciprocalErr
        self.df['recipMean'] = reciprocalMean
        self.df['reci_IP_err'] = reci_IP_err
        
        return Ri
    
    
    def addFilteredIP(self):
        """ Add filtered IP data after IP filtering and pre-processing.
        """
        self.df = pd.merge(self.df, self.filterDataIP[['a','b','m','n']].copy(), how='inner', on=['a','b','m','n'])
        # add Survey.filterDataIP to Survey.df and Survey.dfg
#        self.df['ip'] = np.nan # remove all IP from Survey.df
#        i2keepDf = np.zeros(self.df.shape[0], dtype=bool)
#        i2keepDfg = np.zeros(self.dfg.shape[0], dtype=bool)
#        ipArray = self.filterDataIP[['a','b','m','n']].values
#        dfArray = self.df[['a','b','m','n']].values
#        dfgArray = self.dfg[['a','b','m','n']].values
#        for i in range(self.filterDataIP.shape[0]):
##            print(i, end='')
#            rev1=[0,1,2,3]
#            rev2=[0,1,3,2]
#            rev3=[1,0,2,3]
#            rev4=[1,0,3,2]
#            index1=(dfArray == ipArray[i,rev1]).all(1)
#            index2=(dfArray == ipArray[i,rev2]).all(1)
#            index3=(dfArray == ipArray[i,rev3]).all(1)
#            index4=(dfArray == ipArray[i,rev4]).all(1)
#            index=index1|index2|index3|index4
#            
#            if np.sum(index) > 0:
#                self.df.loc[index, 'ip'] = self.filterDataIP.iloc[i]['ip']
#                i2keepDf[index] = True
#            
#            rev1=[0,1,2,3]
#            rev2=[0,1,3,2]
#            rev3=[1,0,2,3]
#            rev4=[1,0,3,2]
#            index1=(dfgArray == ipArray[i,rev1]).all(1)
#            index2=(dfgArray == ipArray[i,rev2]).all(1)
#            index3=(dfgArray == ipArray[i,rev3]).all(1)
#            index4=(dfgArray == ipArray[i,rev4]).all(1)
#            index=index1|index2|index3|index4
#            
#            if np.sum(index) > 0:
#                self.dfg.loc[index, 'ip'] = self.filterDataIP.iloc[i]['ip']
#                i2keepDfg[index] = True
#            
#        
##        print('df deleted = ', np.sum(~i2keepDf), 'dfg : ', np.sum(~i2keepDfg))
#        self.filterData(i2keepDf)
#        self.dfg = self.dfg[i2keepDfg] 
    
    @staticmethod
    def logClasses3(datax, datay, func, class1=None):
        """ Perform a log class of datay based on datay and applied function
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
    
    def plotError(self, ax=None):
        """ Plot the reciprocal errors.
        
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
#        ax.set_xlabel('Reciprocal Mean [$\Omega$]')
#        ax.set_ylabel('Reciprocal Error [$\Omega$]')
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')  
        ax.set_title('Observed Errors\n')
        if ax is None:
            return fig

    def phaseplotError(self, ax=None): #plotting phase discrepancies over R
        """ Plot the reciprocal phase discrepancies.
        
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
        reciprocalMean = np.abs(self.df['recipMean'].values)
        phase = np.abs(-1.2*self.df['reci_IP_err'].values)
        ax.semilogx(reciprocalMean, phase, 'o')
        ax.set_xlabel(r'LogR [$\Omega$]')
        ax.set_ylabel(r's($\phi$) [mRad]')
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
        
    def plotIPFit(self, ax=None):
        """ Plot the reciprocal phase errors.
        
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
        numbins_ip = 16
        binsize_ip = int(len(self.df['reci_IP_err'])/numbins_ip) 
        Rn = np.abs(self.df['resist'])
        phasedisc = self.df['reci_IP_err']
        error_input_ip = (pd.concat((Rn,phasedisc),axis=1).rename(columns = {'resist':'absRn','reci_IP_err':'Phase_dicrep'})).sort_values(by='absRn').reset_index(drop = True).dropna().query('Phase_dicrep>-25 & Phase_dicrep<25')# Sorting data based on R. the querry is based on environmental IP
        bins_ip = pd.DataFrame(np.zeros((numbins_ip,2))).rename(columns = {0:'R_mean',1:'Phi_dis_STD'})
        for i in range(numbins_ip): # bining 
            ns=i*binsize_ip
            ne=ns+binsize_ip-1
            bins_ip.iloc[i,0] = np.abs(error_input_ip['absRn'].iloc[ns:ne].mean())
            bins_ip.iloc[i,1] = error_input_ip['Phase_dicrep'].iloc[ns:ne].std()  
        bins_ip = bins_ip.dropna()
        coefs_ip= np.linalg.lstsq(np.vstack([np.ones(len(bins_ip.iloc[:,0])), np.log(bins_ip.iloc[:,0])]).T, np.log(bins_ip.iloc[:,1]))[0] # calculating fitting coefficients (a,m)
        R_error_predict_ip = np.exp(coefs_ip[0])*(bins_ip.iloc[:,0]**coefs_ip[1]) # error prediction based of fitted power law model       
        ax.semilogx(error_input_ip['absRn'],np.abs(error_input_ip['Phase_dicrep']), '+', label = "Raw")
        ax.semilogx(bins_ip.iloc[:,0],bins_ip.iloc[:,1],'o',label="bin means")
        ax.plot(bins_ip.iloc[:,0],R_error_predict_ip,'r', label="Power law fit")
        ax.set_ylabel(r's($\phi$) [mRad]')
        ax.set_xlabel(r'R [$\Omega$]')      
        ax.legend(loc='best', frameon=True)
        R2_ip= self.R_sqr(np.log(bins_ip.iloc[:,1]),np.log(R_error_predict_ip))
        a1 = np.around(np.exp(coefs_ip[0]),decimals=3)
        a2 = np.around(coefs_ip[1], decimals=3)
        a3 = np.around(np.exp(coefs_ip[0]),decimals=1)
        a4 = np.around(coefs_ip[1], decimals=1)
        print ('Error model is: Sp(m) = %s*%s^%s (R^2 = %s) \nor simply Sp(m) = %s*%s^%s' % (a1,'R',a2,R2_ip,a3,'R',a4))
        ax.set_title('Multi bin phase error plot\na = %s, b = %s (R$^2$ = %s)' % (a1,a2,R2_ip))
        self.df['PhaseError'] = a1*(np.abs(self.df['resist'])**a2)
        self.df['Phase'] = -1.2*self.df['ip']
        if ax is None:
            return fig   

    def pwlfit(self, ax=None):
        """ Fit an power law to the resistivity data.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
#        self.errTyp = 'pwl'
        if ax is None:
            fig, ax = plt.subplots()        
        numbins = 20
        if 'recipMean' not in self.df.columns:
            self.reciprocal()
        dfg = self.df[self.df['irecip'] > 0]
        binsize = int(len(dfg['recipMean'])/numbins) 
        error_input = np.abs(dfg[['recipMean', 'recipError']]).sort_values(by='recipMean') # Sorting data based on R_avg
        bins = np.zeros((numbins,2))
        for i in range(numbins): # bining 
            ns=i*binsize
            ne=ns+binsize-1
            bins[i,0] = error_input['recipMean'].iloc[ns:ne].mean()
            bins[i,1] = error_input['recipError'].iloc[ns:ne].mean()    
        coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]))[0] # calculating fitting coefficients (a,m)       
        R_error_predict = np.exp(coefs[0])*(bins[:,0]**coefs[1]) # error prediction based of power law model        
        ax.loglog(np.abs(dfg['recipMean']),np.abs(dfg['recipError']), '+', label = "Raw")
        ax.loglog(bins[:,0],bins[:,1],'o',label="bin means")
        ax.plot(bins[:,0],R_error_predict,'r', label="Power law fit")
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')      
        ax.legend(loc='best', frameon=True)
        R2= self.R_sqr(np.log(bins[:,1]),np.log(R_error_predict))
        a1 = np.around(np.exp(coefs[0]),decimals=3)
        a2 = np.around(coefs[1], decimals=3)
        a3 = np.around(np.exp(coefs[0]),decimals=1)
        a4 = np.around(coefs[1], decimals=1)
        print ('Error model is: R_err = %s*%s^%s (R^2 = %s) \nor simply R_err = %s*%s^%s' % (a1,'(R_n/r)',a2,R2,a3,'(R_n/r)',a4))
        ax.set_title('Multi bin power-law plot\n' + r'$\alpha =  %s, \beta = %s$ (R$^2$ = %s)' % (a1,a2,R2))           
        self.df['pwlError'] = a1*(np.abs(self.df['recipMean'])**a2)
        self.errorModel = lambda x : a1*(np.abs(x)**a2)
        if ax is None:
            return fig

    def linfit(self, ax=None):
        """ Fit a linear relationship to the resistivity data.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        
        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
#        self.errTyp = 'lin'
        # # linear fit
        # if 'recipMean' not in self.dfg.columns:
        #     self.reciprocal()
        # recipMean = np.abs(self.dfg['recipMean'].values)
        # recipError = np.abs(self.dfg['recipError'].values)

        # # logspace classes, fit, R^2
        # xm,ystdErr,nbs=self.logClasses3(recipMean,recipError,np.mean, class1=np.arange(-2,2,0.5))# clength=10)
        # inan = ~np.isnan(ystdErr)
        # slope,offset = np.polyfit(np.log10(xm[inan]),np.log10(ystdErr[inan]),1)
        # predy=10**(offset+slope*np.log10(xm[inan]))
        # # R2 makes sens in linear fitting space only and only for the classes
        # r2=1-np.sum((np.log10(predy)-np.log10(ystdErr[inan]))**2)/np.sum((np.log10(predy)-np.mean(np.log10(ystdErr[inan])))**2)
        # print('Simple linear log fit : \n\t offset = {0:0.3f}\n\t slope = {1:.3f}\nR^2 = {2:.4f}'.format(offset, slope, r2))        
        
        # predictRecipError = 10**offset*recipMean**slope
        # self.dfg['linError'] = predictRecipError
        
        # if iplot:
        #     if ax is None:
        #         fig, ax = plt.subplots()
        #     ax.loglog(recipMean, recipError, '+', label='raw')
        #     ax.loglog(xm[inan],ystdErr[inan],'o', label='bin means')
        #     ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        #     ax.set_ylabel('Reciprocal Error [$\Omega$]')
        #     ax.loglog(xm[inan],predy,'r-',label='linear regression')
        #     ax.set_title('loglog graph (R$^2$ = {0:.3f})'.format(r2))
        #     ax.legend()
            
#            fig,ax=plt.subplots()
#            isort = np.argsort(recipMean)
#            ax.plot(recipMean,recipError,'o', label='raw errors')
#            ax.plot(recipMean,predictRecipError,'o', label='predicted')
#            ax.plot(recipMean[isort],2*predictRecipError[isort],'-', label='2*prediction')
#            ax.legend()
#            ax.set_title('Linear fit prediction')
#            ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
#            ax.set_ylabel('Reciprocal Error [$\Omega$]')
##            fig.show()
#            figs.append(fig)
        if ax is None:
            fig, ax = plt.subplots()        
        numbins = 20
        if 'recipMean' not in self.df.columns:
            self.reciprocal()
        dfg = self.df[self.df['irecip'] > 0]
        binsize = int(len(dfg['recipMean'])/numbins) 
        error_input = np.abs(dfg[['recipMean', 'recipError']]).sort_values(by='recipMean') # Sorting data based on R_avg
        bins = np.zeros((numbins,2))
        for i in range(numbins): # bining 
            ns=i*binsize
            ne=ns+binsize-1
            bins[i,0] = error_input['recipMean'].iloc[ns:ne].mean()
            bins[i,1] = error_input['recipError'].iloc[ns:ne].mean()    
        coefs= np.linalg.lstsq(np.vstack([bins[:,0], np.ones(len(bins[:,0]))]).T, bins[:,1])[0] # calculating fitting coefficients (a,m) 
        R_error_predict = ((coefs[0])*(bins[:,0]))+coefs[1] # error prediction based of linear model        
        ax.loglog(np.abs(dfg['recipMean']),np.abs(dfg['recipError']), '+', label = "Raw")
        ax.loglog(bins[:,0],bins[:,1],'o',label="bin means")
        ax.loglog(bins[:,0],R_error_predict,'r', label="Linear fit")
        ax.set_ylabel(r'$R_{error} [\Omega]$')
        ax.set_xlabel(r'$R_{avg} [\Omega]$')      
        ax.legend(loc='best', frameon=True)
        R2= self.R_sqr((bins[:,1]),(R_error_predict))
        a1 = np.around((coefs[0]),decimals=3)
        a2 = np.around(coefs[1], decimals=3)
        a3 = np.around((coefs[0]),decimals=1)
        a4 = np.around(coefs[1], decimals=1)
        print ('Error model is: R_err = %s*%s+%s (R^2 = %s) \nor simply R_err = %s*%s+%s' % (a1,'(R_n/r)',a2,R2,a3,'(R_n/r)',a4))
        ax.set_title('Multi bin Linear plot\n' + r'$m =  %s, b = %s$ (R$^2$ = %s)' % (a1,a2,R2))     
        self.df['linError'] = a1*(np.abs(self.df['recipMean']))+a2
        self.errorModel = lambda x : a1*(np.abs(x))+a2
        if ax is None:
            return fig                  
        
    def linfitStd(self, iplot=False):
        # linear fit with std
        ie = self.irecip > 0
        recipMean0 = np.abs(self.reciprocalMean[ie1])
        recipMean = np.abs(self.reciprocalMean[ie])
        recipError = np.abs(self.reciprocalErr[ie])
        
        # logspace classes, fit, R^2
        xm,ystdErr,nbs=self.logClasses(recipMean,recipError,np.std)
        inan = ~np.isnan(ystdErr)
        slope,offset = np.polyfit(np.log10(xm[inan]),np.log10(ystdErr[inan]),1)
        predy=10**(offset+slope*np.log10(xm[inan]))
        # R2 makes sens in linear fitting space only and only for the classes
        r2=1-np.sum((np.log10(predy)-np.log10(ystdErr[inan]))**2)/np.sum((np.log10(predy)-np.mean(np.log10(ystdErr[inan])))**2)
        print('Simple linear log fit (with std) : \n\t offset = {0:0.3f}\n\t slope = {1:.3f}\nR^2 = {2:.4f}'.format(offset, slope, r2))        
        
        predictRecipError = 10**offset*recipMean**slope
        self.linstdError = predictRecipError
        
        
        if iplot:
            figs = []
            fig,ax=plt.subplots()
#            ax.loglog(recipMean, recipError, '+o') # doesn't make sens because
            # we are plotting the STD of the errors and not the absolute errors
            ax.loglog(xm[inan],ystdErr[inan],'o')
            ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
            ax.set_ylabel('Standard deviation $\sigma$ [$\Omega$]')
            ax.loglog(xm[inan],predy,'r-',label='linear regression')
            ax.set_title('loglog graph (R$^2$ = {0:.3f})'.format(r2))
            fig.show()
            figs.append(fig)
            
            fig,ax=plt.subplots()
            isort = np.argsort(recipMean)
            ax.plot(recipMean,recipError,'o', label='raw errors')
            ax.plot(recipMean,predictRecipError,'o', label='predicted')
            ax.plot(recipMean[isort],2*predictRecipError[isort],'-', label='2*prediction')
            ax.legend()
            ax.set_title('Raw class with stdErr (not in log scale)')
            ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
            ax.set_ylabel('Standard deviation $\sigma$ [$\Omega$]')
            fig.show()
            figs.append(fig)
            
            return figs
        
    def lmefit(self, iplot=True, ax=None):
#        self.errTyp = 'lme'
        print('NOT IMPLEMENTED YET')
#        # fit linear mixed effect model
#        # NEED filterData() before
#        recipMean = np.abs(self.dfg['recipMean'].values)
#        recipError = np.abs(self.dfg['recipError'].values)
#        irecip = self.df['irecip'].values
#        array = self.df[['a','b','m','n']].values
#        
#        ie = irecip > 0
#        data = np.vstack([recipMean, recipError]).T
#        data = np.hstack((data, array[ie]))
#        df = pd.DataFrame(data, columns=['avgR','obsErr','c1','c2','p1','p2'])
#        md = smf.mixedlm('obsErr~avgR', df, groups=df[['c1', 'c2', 'p1', 'p2']])
#        mdf = md.fit()
#        #print(np.min(df['avgR']))
#        #print(np.min(df['obsErr']))
#        #print(np.min(mdf.predict()))
#        print(mdf.summary())
#        
#
#        self.dfg['lmeError'] = mdf.predict()
#        
#        if iplot:
#            if ax is None:
#                fig, ax = plt.subplots()
#
#            ax.loglog(df['obsErr'], mdf.predict(), 'o')
#            ax.loglog([np.min(df['obsErr']),np.max(df['obsErr'])], [np.min(df['obsErr']), np.max(df['obsErr'])], 'r-', label='1:1')
#            ax.grid()
#            ax.legend()
#            ax.set_title('Linear Mixed Effect Model Fit')
#            ax.set_xlabel('Reciprocal Error Observed [$\Omega$]')
#            ax.set_ylabel('Reciprocal Error Predicted [$\Omega$]')

    def heatmap(self,ax=None):
        """ Plot a heatmap based on 
        Orozco, A. F., K. H. Williams, and A. Kemna (2013), Time-lapse spectral induced polarization imaging of stimulated uranium bioremediation, Near Surf. Geophys., 11(5), 531–544, doi:10.3997/1873-0604.2013020)
        
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
#            temp_heatmap_recip_filterN = self.filterDataIP_plot
            if self.filterDataIP.empty:
                temp_heatmap_recip_filterN = self.df[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first')
                dflen = len(self.df)
            else:
                temp_heatmap_recip_filterN = self.filterDataIP[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first')
                dflen = len(self.filterDataIP)
        temp_heatmap_recip_filterN ['Phase'] = temp_heatmap_recip_filterN ['ip']*1.2
        heat_recip_Filter = temp_heatmap_recip_filterN.set_index(['m','a']).ip.unstack(0)     
        if ax is None:
            fig, ax = plt.subplots()  
        else:
            fig = ax.get_figure()             
        m = ax.imshow(heat_recip_Filter, origin='lower',cmap='jet',vmin=0, vmax=25)
        ax.xaxis.set_ticks(np.arange(0,filterDataIP_plotOrig['a'].max()+1,4))
        ax.yaxis.set_ticks(np.arange(0,filterDataIP_plotOrig['m'].max(),4))
        ax.set_ylabel('A',fontsize = 22)
        ax.set_xlabel('M',fontsize = 22)
        ax.tick_params(labelsize=18)
        ax.set_title('%s\n%s measurements' % (self.filt_typ, dflen), fontsize=20)     
        ax.grid(False)
        if self.cbar==True:
            cbhnf = fig.colorbar(m, ax=ax)
            cbhnf.set_label(r'-$\phi$ [mRad]', fontsize=20)
            cbhnf.ax.tick_params(labelsize=18)
        if ax is None:
            return fig
    
    def iprangefilt(self, phimin, phimax):
        """ Filter IP data according to a specified range.
        
        Parameters
        ----------
        phimin : float
            Minimium phase angle [mrad].
        phimax : float
            Maximum phase angle [mrad].
        """
        if self.filterDataIP.empty:
            self.filterDataIP = self.df.query('ip > %s and ip < %s' % (phimin/1.2, phimax/1.2))
        else:
            self.filterDataIP = self.filterDataIP.query('ip > %s and ip < %s' % (phimin/1.2, phimax/1.2))
        self.addFilteredIP()
            
#        temp_data = self.filterDataIP_plotOrig
#        mask = (temp_data.ip < self.phimin) | (temp_data.ip > self.phimax)
#        temp_data.loc[mask, 'ip'] = np.nan
#        self.filterDataIP_plot = temp_data
        
    def removerecip(self):
        
        if self.filterDataIP.empty:
            self.filterDataIP = self.df.query('irecip>=0')
        else:
            self.filterDataIP = self.filterDataIP.query('irecip>=0')
#        self.filterDataIP_plot = self.filterDataIP[['a','m','ip']].drop_duplicates(subset=['a','m'], keep = 'first')
        self.addFilteredIP()

    def removenested(self):
        if self.filterDataIP.empty:
#            self.filterDataIP = self.df.query('m>a & m>b & n>a & n>b')
            temp_data = self.df.copy()
            mask = (temp_data.m < temp_data.b) & (temp_data.m > temp_data.a) | (temp_data.n < temp_data.b) & (temp_data.n > temp_data.a)
            temp_data.loc[mask, 'ip'] = np.nan
            self.filterDataIP = temp_data.dropna()
        else:
#            self.filterDataIP = self.filterDataIP.query('m>a & m>b & n>a & n>b')
            temp_data = self.filterDataIP.copy()
            mask = (temp_data.m < temp_data.b) & (temp_data.m > temp_data.a) | (temp_data.n < temp_data.b) & (temp_data.n > temp_data.a)
            temp_data.loc[mask, 'ip'] = np.nan
            self.filterDataIP = temp_data.dropna()
        self.addFilteredIP()

    def pseudo(self, ax=None, contour=False, log=True, geom=True):
        ''' create true pseudo graph with points and no interpolation
        '''
        array = self.df[['a','b','m','n']].values.astype(int)
        elecpos = self.elec[:,0]
        resist = self.df['resist'].values
        
        if geom: # compute and applied geometric factor
            apos = elecpos[array[:,0]-1]
            bpos = elecpos[array[:,1]-1]
            mpos = elecpos[array[:,2]-1]
            npos = elecpos[array[:,3]-1]
            AM = np.abs(apos-mpos)
            BM = np.abs(bpos-mpos)
            AN = np.abs(apos-npos)
            BN = np.abs(bpos-npos)
            K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
            resist = resist*K
            
        if log:
            resist = np.sign(resist)*np.log10(np.abs(resist))
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
        else:
            label = r'$\rho_a$ [$\Omega.m$]'
        
        cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
            + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
            + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
        ypos = - np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
        
        if contour is False:
            if ax is None:
                fig, ax = plt.subplots()
            else:
                fig = ax.get_figure()
            cax = ax.scatter(xpos, ypos, c=resist, s=70)#, norm=mpl.colors.LogNorm())
            cbar = fig.colorbar(cax, ax=ax)
            cbar.set_label(label)
            ax.set_title('Pseudo Section')
    #        fig.suptitle(self.name, x= 0.2)
#            fig.tight_layout()
        
        if contour:
            from matplotlib.mlab import griddata
            def grid(x, y, z, resX=100, resY=100):
                "Convert 3 column data to matplotlib grid"
                xi = np.linspace(min(x), max(x), resX)
                yi = np.linspace(min(y), max(y), resY)
                Z = griddata(x, y, z, xi, yi, interp='linear')
                X, Y = np.meshgrid(xi, yi)
                return X, Y, Z
            X, Y, Z = grid(xpos, ypos, resist)
            if ax is None:
                fig, ax = plt.subplots()
            figsize=(10,6)
            cax = ax.contourf(X,Y,Z)
#            cbar = fig.colorbar(cax, ax=ax)
#            cbar.set_label(label)
            ax.set_title('Pseudo Section')
#            fig.suptitle(self.name, x= 0.2)
#            fig.tight_layout()

        if ax is None:
            return fig
    
    def pseudoIP(self, ax=None, contour=False): #IP pseudo section
        """ Create pseudo section of IP data with points (default)
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, the graph is plotted along this axis.
        contour : bool
            If True, use filled contour instead of points in the pseudosection.
        
        Returns
        -------
        fig : matplotlib figure
            If `ax` is not specified, the method returns a figure.
        """
        array = self.df[['a','b','m','n']].values
        elecpos = self.elec[:,0]
        ip = self.df['ip'].values            

        label = r'$\phi$ [mRad]'
        
        cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
            + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
        pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
            + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
        xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
        ypos = - np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
        
        if contour is False:
            if ax is None:
                fig, ax = plt.subplots()
            else:
                fig = ax.get_figure()
            cax = ax.scatter(xpos, ypos, c=ip, s=70)#, norm=mpl.colors.LogNorm())
            cbar = fig.colorbar(cax, ax=ax)
            cbar.set_label(label)
            ax.set_title('IP pseudo Section')
    #        fig.suptitle(self.name, x= 0.2)
#            fig.tight_layout()
        
        if contour:
            from matplotlib.mlab import griddata
            def grid(x, y, z, resX=100, resY=100):
                "Convert 3 column data to matplotlib grid"
                xi = np.linspace(min(x), max(x), resX)
                yi = np.linspace(min(y), max(y), resY)
                Z = griddata(x, y, z, xi, yi, interp='linear')
                X, Y = np.meshgrid(xi, yi)
                return X, Y, Z
            X, Y, Z = grid(xpos, ypos, ip)
            if ax is None:
                fig, ax = plt.subplots()
            figsize=(10,6)
            cax = ax.contourf(X,Y,Z)
            cbar = fig.colorbar(cax, ax=ax)
            cbar.set_label(label)
            ax.set_title('IP pseudo Section')
#            fig.suptitle(self.name, x= 0.2)
#            fig.tight_layout()

        if ax is None:
            return fig
    
    def write2protocol(self, outputname='', errTyp='none', errTot=False, ip=False, errTypIP='none'):
        """ Write a protocol.dat file for R2 or cR2.
        
        Parameters
        ----------
        outputname : str, optional
            Path of the output file.
        errTyp : str, optional
            If `none` no error columns will be added. Other options are : 
                'lin','lme','pwl'or 'obs'.
        errTot : bool, optional
            If `True`, the modelling error will be added to the error from the
            error model to form the *total error*.
        ip : bool, optional
            If `True` and IP columns will be added to the file.
        errTypIP : str, optional
            Needs `ip=True` to be taken into account. Specify the string of
            the error model to apply for IP data.
        """
        ie = self.df['irecip'].values > 0 # consider only mean measurement (not reciprocal)
        haveReciprocal = all(self.df['irecip'].values == 0)
        if haveReciprocal == False: # so we have reciprocals
            x = self.df[['a','b','m','n']].values[ie,:].astype(int)
            xx = np.c_[1+np.arange(len(x)), x]
            protocol = pd.DataFrame(xx, columns=['num','a','b','m','n'])
            dfg = self.df[self.df['irecip'] > 0] 
            protocol['R'] = dfg['recipMean'].values
            if ip == True:
#                protocol['R'] = np.abs(protocol['R'])
                protocol['Phase'] = -1.2*dfg['ip'].values # "-1.2" factor is for IRIS syscal instrument
            if errTyp != 'none':
                if errTyp == 'obs':
                    protocol['error'] = self.df['recipError'].values[ie]
                if errTyp =='lme':
                    protocol['error'] = dfg['lmeError'].values
                if errTyp == 'lin':
                    protocol['error'] = dfg['linError'].values
                if errTyp == 'pwl':
                    protocol['error'] = dfg['pwlError'].values
        #            error = self.linStdError
                if errTot == True:
                    if len(self.modError) == 0:
                        print('ERROR : you must specify a modelling error')
                    else:
                        protocol['error'] = np.sqrt(protocol['error']**2 + self.modError[ie]**2)
            if errTypIP != 'none':  # or == 'pwlip'
                if 'PhaseError' not in self.df.columns: # TO BE DELETED
                    dfg['PhaseError'] = 0.1 # TO BE DELTED
                protocol['ipError'] = self.df['PhaseError'].values[ie]
                
        else: # why don't they have reciprocals my god !!
            x = self.df[['a','b','m','n']].values.astype(int)
            xx = np.c_[1+np.arange(len(x)), x]
            protocol = pd.DataFrame(xx, columns=['num','a','b','m','n'])
            protocol['R'] = self.df['resist'].values
        
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(len(protocol)) + '\n')
            with open(outputname, 'a') as f:
                protocol.to_csv(f, sep='\t', header=False, index=False)
                

    def manualFilter(self, ax=None):
        """ Manually filters the data visually.
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, the graph is plotted along the axis.
        
        Returns
        -------
            If `ax` is not None, a matplotlib figure is returned.
        """
        array = self.df[['a','b','m','n']].values
        if 'recipError' in self.df.columns:
            resist = self.df['recipError'].values
            print('set to reciprocal')
        else:
            resist = np.ones(self.df.shape[0])
        spacing = np.mean(np.diff(self.elec[:,0]))
        pseudo(array, resist, spacing, ax=ax, label='Reciprocal Error [Ohm.m]')
        
    def dca(self, dump=print):
        ''' execute DCA filtering
        '''
        if self.filterDataIP.empty:
            self.filterDataIP = DCA(self.df, dump=dump)
        else:
            self.filterDataIP = DCA(self.filterDataIP, dump=dump)
        self.addFilteredIP()
        dump(100)
        
#TODO put this in function as Survey.manualFilter() ? or make it more general
def pseudo(array, resist, spacing, name='', ax=None, figsize=(12,3), contour=False, log=True, geom=False, label=''):
    #figsize=(12,3)
    """ create true pseudo graph with points and no interpolation
    """
    geom = False
    log = False
    nelec = np.max(array)
    elecpos = np.arange(0, spacing*nelec, spacing)
    resist = resist
    
    if geom: # compute and applied geometric factor
        apos = elecpos[array[:,0]-1]
        bpos = elecpos[array[:,1]-1]
        mpos = elecpos[array[:,2]-1]
        npos = elecpos[array[:,3]-1]
        AM = np.abs(apos-mpos)
        BM = np.abs(bpos-mpos)
        AN = np.abs(apos-npos)
        BN = np.abs(bpos-npos)
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        resist = resist*K
        
    if log:
        resist = np.sign(resist)*np.log10(np.abs(resist))
#        label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
#    else:
#        label = r'$\rho_a$ [$\Omega.m$]'
    
    cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
        + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
    pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
        + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
    xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
    ypos = - np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
    
    
    def onpick(event):
        # TODO single doesn't want to change the electrode selection
        if lines[event.artist] == 'data':
            print('onpick event', event.ind[0])
            print(iselect[event.ind[0]])
            if iselect[event.ind[0]] == True:
                print('set to false')
                iselect[event.ind[0]] = False
            else:
                iselect[event.ind[0]] = True
        
        if lines[event.artist] == 'elec':
            print('onpick2', event.ind[0])
            ie = (array == (event.ind[0]+1)).any(-1)
            if all(iselect[ie] == True):
                iselect[ie] = False
            else:
                iselect[ie] = True
            if eselect[event.ind[0]] == True:
                eselect[event.ind[0]] = False
            else:
                eselect[event.ind[0]] = True
            elecKilled.set_xdata(elecpos[eselect])
            elecKilled.set_ydata(np.zeros(len(elecpos))[eselect])
        print('update canvas')
        killed.set_xdata(x[iselect])
        killed.set_ydata(y[iselect])
        killed.figure.canvas.draw()
            
            
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    caxElec, = ax.plot(elecpos, np.zeros(len(elecpos)), 'ko', picker=5)
    cax = ax.scatter(xpos, ypos, c=resist, marker='o', picker=5)
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label)
    cax.figure.canvas.mpl_connect('pick_event', onpick)
    #for i in range(int(len(array)/2)):
    #    ax.text(xpos[i], ypos[i], str(array[i,:]), fontsize=8)
    
    line = cax
    killed, = line.axes.plot([],[],'rx')
    elecKilled, = line.axes.plot([],[],'rx')
#    x = np.array(line.get_xdata())
#    y = np.array(line.get_ydata())
    x = line.get_offsets()[:,0]
    y = line.get_offsets()[:,1]
    
    # TODO we need to get those variables out of the function
    iselect = np.zeros(len(y),dtype=bool)
    eselect = np.zeros(len(elecpos), dtype=bool)
    
    lines = {line:'data',caxElec:'elec',killed:'killed'}





"""        
    def addModError(self, fname):
        # add the modelling error for both normal and reciprocal quadrupoles
        x = np.genfromtxt(fname)
        err = np.zeros(self.ndata)*np.nan
        for i in range(0, self.ndata):
            index = (self.array[i,:] == x[:,:-1]).all(1)
            if np.sum(index) == 1:
                err[i] = x[index,-1]
            else:
                print('ERROR : same quadrupoles has', np.sum(index),' modelling error')
        print(np.sum(self.reciprocalErr < err), '/', self.ndata,
                     'modelling errors bigger than observed errors')
        self.modError = err
        
    def plotError(self):
        fig, ax = plt.subplots()
        ax.loglog(np.abs(self.reciprocalMean), self.reciprocalErr, 'o')
        ax.set_xlabel('Reciprocal Mean [$\Omega$]')
        ax.set_ylabel('Reciprocal Error [$\Omega$]')
        return fig
    
    def toTable(self, outputname=''):
        # return formated csv with normal and reciprocal
        narray = self.array[self.irecip > 0]
        nresist = self.resist[self.irecip > 0]
        nerror = self.reciprocalErr[self.irecip > 0]
        rarray = self.array[self.irecip < 0]
        rresist = self.resist[self.irecip < 0]
        rerror = self.reciprocalErr[self.irecip <0]
        # align them
        ipos = self.irecip[self.irecip > 0]
        ineg = self.irecip[self.irecip < 0]
        isortp = np.argsort(ipos)
        isortn = np.argsort(np.abs(ineg))
        # sort them
        table = np.hstack([narray[isortp,:],
                           nresist[isortp, None],
                           nerror[isortp, None],
                           rarray[isortn,:],
                           rresist[isortn, None],
                           rerror[isortn, None]])
        
        if outputname != '':
            np.savetxt(outputname, table, fmt='%d %d %d %d %f %f %d %d %d %d %f %f')
        
        return table
    
    
    def estimateError(self):
        # calculate reciprocal error and reciprocal mean
        #Ri=self.reciprocal()
        Ri = self.irecip        
        
        # create a logical index
        self.Ri=Ri>0
                    
        # remove data where mean(normal,reciprocal) > 100 
        #reciprocal error
        goodRecip=np.abs(self.reciprocalErr/self.reciprocalMean)<1    
        print("reciprocal error > 100% : "+str(len(self.R[~goodRecip])))        
        self.R[~goodRecip]=np.nan
        self.reciprocalMean[~goodRecip]=np.nan 
        self.reciprocalErr[~goodRecip]=np.nan
        
        # to use % error, use self.reciprocalErrRel
        self.xm,self.ystdErr,nbs=self.logClasses(np.abs(self.reciprocalMean),np.abs(self.reciprocalErr),np.std)
        # sometimes all is no comprised into the bin log scales so
        # it could happen that some measurements are missing (few)
        
        #print 'nbs=',np.sum(nbs)
        #print 'number of data per class : ', nbs
        nbs=nbs+1 # to not have weight of zero
        nbs=nbs+1*np.arange(len(nbs),0,-1) # more weight on small R
        nbs=nbs/np.sum(nbs) # to have relative weigth
        inan=~np.isnan(self.ystdErr)
        
        #print len(self.reciprocalErr[~np.isnan(self.reciprocalErr)])
        #print len(goodRecip)
        
        # --------------------------linear model fitting
        # (with standard deviation of error class !)
        fig,ax=plt.subplots()
        ax.plot(self.xm[inan],self.ystdErr[inan],'bo')
        ax.set_title('Raw class with stdErr (not in log scale)')
        ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        ax.set_ylabel('Standard deviation $\sigma$ [$\Omega$]')
        #fig.show()
        
        fig,ax=plt.subplots()
        ax.loglog(self.xm[inan],self.ystdErr[inan],'go')
        ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        ax.set_ylabel('Standard deviation $\sigma$ [$\Omega$]')
        # linear fit in log space
        slope,offset=np.polyfit(np.log10(self.xm[inan]),np.log10(self.ystdErr[inan]),1)
        print('linear fit in log space : \n\t slope = ' + str(slope) + '\n\t offset = ' + str(offset))
        
        predy=10**(offset+slope*np.log10(self.xm[inan]))
        # R2 makes sens in linear fitting space only
        r2=1-np.sum((np.log10(predy)-np.log10(self.ystdErr[inan]))**2)/np.sum((np.log10(predy)-np.mean(np.log10(self.ystdErr[inan])))**2)
        print('R^2 = '+str(r2))
        ax.loglog(self.xm[inan],predy,'r-',label='linear regression')
        ax.set_title('loglog graph (R$^2$ = {0:.3f})'.format(r2))
        # TODO add nb of points in the class on graph and equation
        # coming back to non log space
        a=10**offset
        b=slope
        #fig.show()
        
        def efunc(x,a,b):
            return a*x**b
        
        self.errModel=efunc

        # graph with 2*std as in BertFile
        fig, ax = plt.subplots()
        ax.plot(self.reciprocalMean, self.reciprocalErr, 'ro')
        ax.plot(self.xm[inan], predy*2, 'b-')
        ax.plot(self.xm[inan], -predy*2, 'b-')
        ax.plot(self.xm[inan], predy, 'g-')
        ax.plot(self.xm[inan], -predy, 'g-')
        ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        ax.set_ylabel('Reciprocal Error [$\Omega$]')
        #fig.show()
        
        
        '''
        # exponential model fitting (curve_fit)
        def func(x,a,b):
            return a*x**b
        
        from scipy.optimize import curve_fit
        popt0=[offset,slope]
        popt,pcov=curve_fit(func,self.xm[inan],self.ystdErr[inan])
        print(popt)
        fig,ax=plt.subplots()
        ax.loglog(self.xm[inan],self.ystdErr[inan],'bo')
        ax.loglog(self.xm[inan],func(self.xm[inan],*popt),'r') # *popt expand terms
        #ax.loglog(self.xm[inan],popt[0]*self.xm[inan]**popt[1],'p-',label='exponential regression')
        fig.show()
        
        # exponential model fitting (minimize)
        def objfunc(param): # in log space
            return np.sum((self.ystdErr[inan]-param[0]*self.xm[inan]**param[1])**2)
                
        from scipy.optimize import minimize
        bnds=((0.0001,1.0),(0.0001,1.0))
        res=minimize(objfunc,[0.1,0.1],method='TNC',options={'xtol':1e-8})
        print(res.x)
        fig,ax=plt.subplots()
        ax.loglog(self.xm[inan],self.ystdErr[inan],'bo')
        ax.loglog(self.xm[inan],res.x[0]*self.xm[inan]**res.x[1],'r')
        #ax.loglog(self.xm[inan],popt[0]*self.xm[inan]**popt[1],'p-',label='exponential regression')
        fig.show()
        '''
        
        # ------------------simple fit with reciprocal Error (not stdErr)
        # reciprocalErr VS meanResistivity (linear fit in logscale)
        fig, ax = plt.subplots()
        ax.loglog(self.reciprocalMean, self.reciprocalErr, 'b+', label='raw')
        ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        ax.set_ylabel('Reciprocal Error [$\Omega$]')
        
        # and now in logclasses
        xlog, ylog, nbs = self.logClasses(np.abs(self.reciprocalMean), np.abs(self.reciprocalErr),np.mean)
        inan = ~np.isnan(ylog)
        slope, offset = np.polyfit(np.log10(xlog[inan]), np.log10(ylog[inan]), 1)
        ylogPredicted = 10**(offset + slope*np.log10(xlog))
        print(ylogPredicted.shape)
        r2=1-np.sum((np.log10(ylogPredicted[inan])-np.log10(ylog[inan]))**2)/np.sum((np.log10(ylogPredicted[inan])-np.mean(np.log10(ylog[inan])))**2)
        print('Simple linear log fit : \n\t offset = {0:0.3f}\n\t slope = {1:.3f}\nR^2 = {2:.4f}'.format(offset, slope, r2))        
        #fig, ax = plt.subplots()
        ax.loglog(xlog[inan], ylog[inan], 'go', label='binned')
        ax.loglog(xlog[inan], ylogPredicted[inan], 'r-')
        ax.set_title('Linear fit in logscale (R$^2$ = {0:.3f})'.format(r2))
        #ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        #ax.set_ylabel('Reciprocal Error [$\Omega$]')
        ax.legend()
        #fig.show()
        
        # graph with 2*std as in BertFile
        ypredict = 10**(offset + slope*np.log10(self.reciprocalMean))
        self.ypredict = ypredict
        self.iout = np.zeros(len(self.reciprocalMean), dtype=bool)
        self.iout[np.abs(self.reciprocalErr) > ypredict*2] = True
        
        fig, ax = plt.subplots()
        ax.plot(self.reciprocalMean, self.reciprocalErr, 'r+')
        ax.plot(self.reciprocalMean[self.iout], self.reciprocalErr[self.iout], 'k+')
        ax.plot(xlog[inan], ylogPredicted[inan]*2, 'b-', label='2x')
        ax.plot(xlog[inan], -ylogPredicted[inan]*2, 'b-')
        ax.plot(xlog[inan], ylogPredicted[inan], 'g-', label='1x')
        ax.plot(xlog[inan], -ylogPredicted[inan], 'g-')
        ax.set_title('Linear fit prediction')
        ax.set_xlabel('Mean of transfer resistance [$\Omega$]')
        ax.set_ylabel('Reciprocal Error [$\Omega$]')
        ax.legend()
        #fig.show()
        
        return a,b # from linear fitting
    
    def output(self):
        ie = self.irecip > 0
        return np.vstack([self.array[ie], self.reciprocalMean[ie], self.predictErr[ie]])
                          
        
    def predictError(self,a,b):
        # predict error based on error model
        self.predictErr=self.errModel(np.abs(self.resist),a,b)
        
        fig,ax=plt.subplots()
        ax.set_title('predicted Error')
        ax.loglog(self.xm,self.ystdErr,'bo', label='measured')
        ax.loglog(np.abs(self.resist),self.predicErr,'ro', label='predicted')
        ax.legend()
        ax.set_xlabel('Transfer Resistance [$\Omega$]')
        ax.set_ylabel('Error [%]')
        return fig
    
    
    def write2protocol(self, outputname='', errTyp='lme', errTot=False):
        ie = self.irecip > 0 # consider only mean measurement (not reciprocal)
        n = np.sum(ie)
        if errTyp=='obs':
            error = self.reciprocalErr[ie]
        if errTyp =='lme':
            error = self.lmeError
        if errTyp == 'lin':
            error = self.linError
#            error = self.linStdError
        if errTot == True:
            if len(self.modError) == 0:
                print('ERROR : you must specify a modelling error')
            else:
                error = np.sqrt(error**2 + self.modError[ie]**2)
        arr = self.array[ie,:]
        res = self.reciprocalMean[ie]
        
        if outputname != '':
            with open(outputname, 'w') as f:
                f.write(str(n) + '\n')
                for i in range(n):
                    f.write('{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}\t{4:.0f}\t{5:f}\t{6:f}\n'.format(
                            i+1, arr[i,0], arr[i,1], arr[i,2], arr[i,3], res[i], error[i]))
        return np.hstack([arr, res.reshape((len(res), 1)), error.reshape((len(error), 1))])
        
    
    def write2txt(self, outputname='', error=False, errTyp='obs'):
        ''' write protocol.dat files with errors or not
        '''
        if outputname=='':
            outputname='../output/'+self.name+'.txt'
        with open(outputname,'w') as f:
            f.write(str(self.ndata))
            f.write('\n')
            for i in range(0,self.ndata):
                f.write('{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}\t{4:.0f}'.format(i+1,self.array[i,0],self.array[i,1],self.array[i,2],self.array[i,3]))
                if error:
                    if errTyp == 'obs':
                        err = self.reciprocalErr
#                    elif errTyp == 'rel':
#                        err = 100*self.reciprocalErrRel # in percent
                    elif errTyp == 'lme':
                        err = self.lmeError
                    elif errTyp == 'lin':
                        err = self.ypredict
                    
                    f.write('\t{0:f}\t{1:f}\n'.format(self.resist[i], err[i])) # if self.reciprocalMean -> lost sign !
                else:
                    f.write('\t{0:f}\n'.format(self.resist[i]))
    
    def write2bert(self, outputname='', position=None):
        if outputname=='':
            outputname='../output/'+self.name+'.dat'
        if position is None:
            position = self.electrodes
        with open(outputname,'w') as f:
            f.write(str(len(position)) + ' # Number of electrodes\n# x y z\n')
            for j in range(0, position.shape[0]):
                f.write('{0:.3f}\t{1:.3f}\t{2:.3f}\n'.format(position[j,0], position[j,1], position[j,2]))
            
            if (len(self.reciprocalErr)>0)&(np.sum(self.m) > 0):
                print('PRINT with IP and ERROR')
                f.write(str(self.ndata) + '# Number of data\n# a b m n R err ip\n')
                for i in range(0,self.ndata):
                    f.write('{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}\t{4:f}\t{5:f}\t{6:f}\n'.format(self.array[i,0],self.array[i,1],self.array[i,2],self.array[i,3],self.resist[i], self.reciprocalErr[i], self.m[i]))
            
            elif (len(self.reciprocalErr)>0)&(np.sum(self.m) == 0):
                print('PRINT with ERROR')
                f.write(str(self.ndata) + '# Number of data\n# a b m n R err\n')
                for i in range(0,self.ndata):
                    f.write('{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}\t{4:f}\t{5:f}\n'.format(self.array[i,0],self.array[i,1],self.array[i,2],self.array[i,3],self.resist[i], self.reciprocalErr[i]))
            
            else:
                print('PRINT default')
                f.write(str(self.ndata) + '# Number of data\n# a b m n R\n')
                for i in range(0,self.ndata):
                    f.write('{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}\t{4:f}\n'.format(self.array[i,0],self.array[i,1],self.array[i,2],self.array[i,3],self.resist[i]))
    
    
    def keepRecip(self):
        ''' filter out all data without reciprocal (dummy measurements)
        ''''
        irecip = self.reciprocal()
        ikeep = irecip > 0
        self.filterData(ikeep)
                                                   
    def append(self, k):
        self.array = np.r_[self.array, k.array]
        self.resist = np.r_[self.resist, k.resist]
        self.i = np.r_[self.i, k.i]
        self.vp = np.r_[self.vp, k.vp]
        self.dev = np.r_[self.dev, k.dev]
        self.sp = np.r_[self.sp, k.sp]
        self.m = np.r_[self.m, k.m]
        self.ndata = len(self.array)
        self.irecip = self.reciprocal()
    
        
    def filterData(self,i2keep):
        # applied filter to self.array before writing measurment to
        # inversion files
        
        newarray=np.copy(self.array)
        newresist=np.copy(self.resist)
        #newcontactr=np.copy(self.contactR)

        #i2keep=self.irecip&self.is2b&self.Ri
        #i2keep=np.where(i2keep)
        
        self.array=newarray[i2keep,:] #ffilter(self.array,i2keep)
        self.resist=newresist[i2keep] #ffilter(self.resist,i2keep)
        #self.contactR=newcontactr[i2keep]
        self.i = self.i[i2keep]
        self.vp = self.vp[i2keep]
        self.irecip = self.irecip[i2keep]
        
        if len(self.reciprocalErr)>0:
            newreciperr=np.copy(self.reciprocalErr)
            newrecipmean=np.copy(self.reciprocalMean)
            newreciperrrel=np.copy(self.reciprocalErrRel)
            self.reciprocalErrRel=newreciperrrel[i2keep]
            self.reciprocalErr=newreciperr[i2keep]
            self.reciprocalMean=newrecipmean[i2keep]
        
        if len(self.modError) > 0:
            self.modError = self.modError[i2keep]
            
        if self.fip>0:
            newip=np.copy(self.ip)
            self.ip=newip[i2keep,:] #ffilter(self.ip,i2keep)
        
        self.ndata=len(self.array)
        print(str(len(np.where(~i2keep)[0])) + ' data deleted | ' + str(self.ndata) + ' data kept.')
        
        self.irecip = self.reciprocal() # to rebuild self.recipMean and recipError
        
        
    def filterdip(self,elec): # deleted specific elec data
        #index=(self.array==any(elec)).any(-1)
        index = (self.array == elec[0]).any(-1)
        for i in range(1,len(elec)):
            index = index | (self.array == elec[i]).any(-1)
        self.filterData(~index)
    
    def plotData(self):
        fig, ax = plt.subplots()
        ax.plot(self.resist, 'o')
        index = self.dev > 1
        print(self.array[index,:])
        ax.plot(np.where(index)[0], self.resist[index], 'o', label='dev > 1')
        ax.legend()
        ax.set_ylabel('Transfer Resistance [$\Omega$]')
        ax.set_xlabel('Measurements')
        fig.show()
        return fig


    def errorDist(self):
        ''' plot the reciprocal relative error distribution and fit a 
        gaussian model
        '''
        fig, ax = plt.subplots()
        logclassx, logclassy, nbs = self.logClasses(self.reciprocalErrRel,
                    self.reciprocalErrRel, np.mean)
        p = ax.plot(np.log10(np.abs(logclassx)), nbs, 'o-', label='reciprocalErrRel')
        inan = ~np.isnan(nbs)
        xx = np.log10(logclassx)[inan]
        popt, pcov = curve_fit(gauss, xx, nbs[inan], [5, 20,-2, 0.5])
        ax.plot(xx, gauss(xx, popt[0], popt[1], popt[2], popt[3]), '--',
                          color = p[0].get_color(), label='({0:.2f}, {1:.2f})'.format(popt[2], popt[3]))
        ax.set_title('PDF')
        ax.set_xlabel('log$_{10}$(Relative Reciprocal Error) [%]')
        ax.set_ylabel('Frequency')
        ax.grid()
        ax.legend()
        
        return fig
    


    
"""
       



        
#%% test code
#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#s = Survey('api/test/syscalFile.csv', ftype='Syscal')
#s = Survey('api/test/rifleday8.csv', ftype='Syscal')
#s.dca()
#s.addFilteredIP()
#s.pwlfit()
#s.plotIPFit()
#s.write2protocol('kk.txt', ip=True, errTyp='pwl', errTypIP='pwl')
#s = Survey('test/syscalFileNormalOnly.csv', ftype='Syscal')
#s.addData('test/syscalFileReciprocalOnly.csv', ftype='Syscal')
#fig, ax = plt.subplots()
#fig.suptitle('kkkkkkkkkkkkkk')
#s.plotError(ax=ax)
#s.manualFilter()
#s.pseudo(contour=True)
#s.linfit()
#s.pwlfit()
#s.plotIPFit()
#s.write2protocol('api/test/protocol.dat', errTyp='lin', ip=True, errTypIP='pwl')
#s.dca()
        
