#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 21:01:53 2017

@author: Sina
"""

import numpy as np
import os
import pandas as pd
from scipy.stats.kde import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn.metrics import mean_squared_error as mse
#%% data directory
##### no GUI #####
#fdir = '/Volumes/Other/Dropbox/Hydrogeophysics Workshop/Codes/data_test/TusseyMtn/' #for mac
directory = 'F:\Sina\Google Drive\Hydrogeophysics Workshop\Codes\data_test\leadingRidge\delete'
#directory = 'F:\Sina\Dropbox\Hydrogeophysics Workshop\Codes\data_test\TusseyMtn'
b="\ ".strip()
fdir = directory+b
os.chdir(fdir)
coordinates = pd.read_csv('Electrode_coordinates_leading_ridge.csv')
#coordinates = pd.read_csv('Electrode_coordinates_TusseyMtn.csv')
numnorm = 17 #number of normal measurements
numrecip = 5 #number of reciprocal measurements
#%% Universal functions
def res(data): # calculating Resistance (ohm)
    R = data[:,5]/data[:,6]
    return R
def fit(x,a,m): # power-law fit (based on previous studies)
    return a*(x**m)
def fit_lin(x,m,b):
    return m*x+b
def positive_test (Dcurve,DecayTime): #calculating decay curve trend - positive trends are bad data
    DC_slope = np.zeros(np.shape(Dcurve)[0])
    for i in range(np.shape(Dcurve)[0]):     
        DecayCoefs = np.linalg.lstsq(np.vstack([DecayTime, np.ones(len(DecayTime))]).T,Dcurve.iloc[i,:])[0]
        DC_slope[i]=DecayCoefs[0]
    return DC_slope
def R_sqr (y,y_predict): # calculating R squared value to measure fitting accuracy
    	rsdl = y - y_predict
    	ss_res = np.sum(rsdl**2)
    	ss_tot = np.sum((y-np.mean(y))**2)
    	R2 = 1-(ss_res/ss_tot)
    	R2 = np.around(R2,decimals=4)
    	return R2
def linear_coefs (x,y): #linear fit parameteres for decay curve
    data = np.concatenate((np.log(x)[:,None],np.log(y)[:,None]),axis=1)
    data_no_nan = data[~np.isnan(data).any(axis=1)]
    coefs = np.linalg.lstsq(np.vstack([data_no_nan[:,0], np.ones(len(data_no_nan[:,0]))]).T,data_no_nan[:,1])[0]
    return coefs
def find_common (N,R): #finding common measurments between normal and reciprocal datasets
    	df_1 = N
    	df_2 = R
    	common_1 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Mr', 'Nr', 'Ar', 'Br'])
    	common_2 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Nr', 'Mr', 'Ar', 'Br'])
    	common_2['Rr'] = common_2['Rr']*-1
    	common_3 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Mr', 'Nr', 'Br', 'Ar'])
    	common_3['Rr'] = common_3['Rr']*-1
    	common_4 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Nr', 'Mr', 'Br', 'Ar'])
    	common = pd.concat([common_1, common_2, common_3, common_4])
    	return common
def commons_normal (data1,data2,num_data): #finding common measurements among time-lapse normal datasets (for time-lapse inversion)
    shape_array = pd.DataFrame(np.zeros([len(data1),2])).rename(columns={0:'length',1:'order'})
    for i in range(len(data1)):
        shape_array.iloc[i,:]= [len(data1[i]),i]
    common_index = shape_array.sort_values(by='length').reset_index(drop = True).iloc[0,1]
    common_final = []
    for i in range(num_data):
        temp_common = pd.merge(data2[i], data1[int(common_index)][['An', 'Bn', 'Mn', 'Nn']], how = 'inner').drop_duplicates(subset = ['An', 'Bn', 'Mn', 'Nn'])
        common_final.append(temp_common)
    min_num_array = pd.DataFrame(np.zeros([len(common_final),2])).rename(columns={0:'length',1:'order'})
    for i in range(len(common_final)):
        min_num_array.iloc[i,:]= [len(common_final[i]),i]
    min_nums = min_num_array.sort_values(by='length').reset_index(drop = True).iloc[0,0]
    return (common_final,min_nums)
#%% building raw dataframes
normal_sets = []
recip_sets = []
for i in range(numnorm): #loading normal datasets numbered as "norm i.txt" i is an integer
    m=i+1
    readN = np.loadtxt('norm (%s).txt' % (m), delimiter = ' ', skiprows=1, usecols = (5,6,7,8,10,13,14))
    dfrawN = pd.DataFrame(np.concatenate((readN,res(readN)[:,None]),axis=1)).rename(columns = {0:'An',1:'Bn',2:'Mn',3:'Nn',4:'Dev_n',5:'V_n',6:'I_n',7:'Rn'}).query('(I_n>10) & (V_n>-15000) & (V_n<15000)') #initial voltage and current filtering
    normal_sets.append(dfrawN)
    
for i in range(numrecip): #loading reciprocal datasets numbered as "norm i.txt" i is an integer
    mn=i+1
    readR = np.loadtxt('recip (%s).txt' % (mn), delimiter = ' ', skiprows=1, usecols = (5,6,7,8,10,13,14))
    dfrawR = pd.DataFrame(np.concatenate((readR,res(readR)[:,None]),axis=1)).rename(columns= {0:'Ar',1:'Br',2:'Mr',3:'Nr',4:'Dev_r',5:'V_r',6:'I_r',7:'Rr'}).query('(I_r>10) & (V_r>-15000) & (V_r<15000)')
    recip_sets.append(dfrawR)

#%% Removing duplicates and finding common measurements for normal and reciprocal measurements + removing NaN values
reciprocity_1 = find_common(normal_sets[0],recip_sets[0]).drop_duplicates(subset=['An','Bn', 'Mn', 'Nn'], keep = 'first').drop_duplicates(subset=['Ar','Br', 'Mr', 'Nr'], keep = 'first').replace([np.inf, -np.inf], np.nan).dropna(subset=['Rn','Rr'])
reciprocity_2 = find_common(normal_sets[13],recip_sets[1]).drop_duplicates(subset=['An','Bn', 'Mn', 'Nn'], keep = 'first').drop_duplicates(subset=['Ar','Br', 'Mr', 'Nr'], keep = 'first').replace([np.inf, -np.inf], np.nan).dropna(subset=['Rn','Rr'])
reciprocity_3 = find_common(normal_sets[14],recip_sets[2]).drop_duplicates(subset=['An','Bn', 'Mn', 'Nn'], keep = 'first').drop_duplicates(subset=['Ar','Br', 'Mr', 'Nr'], keep = 'first').replace([np.inf, -np.inf], np.nan).dropna(subset=['Rn','Rr'])
reciprocity_4 = find_common(normal_sets[15],recip_sets[3]).drop_duplicates(subset=['An','Bn', 'Mn', 'Nn'], keep = 'first').drop_duplicates(subset=['Ar','Br', 'Mr', 'Nr'], keep = 'first').replace([np.inf, -np.inf], np.nan).dropna(subset=['Rn','Rr'])
reciprocity_5 = find_common(normal_sets[16],recip_sets[4]).drop_duplicates(subset=['An','Bn', 'Mn', 'Nn'], keep = 'first').drop_duplicates(subset=['Ar','Br', 'Mr', 'Nr'], keep = 'first').replace([np.inf, -np.inf], np.nan).dropna(subset=['Rn','Rr'])

#%% appending [random] reciprocal datasets
reciprocity_sets = []
reciprocity_sets.append(reciprocity_1)
reciprocity_sets.append(reciprocity_2)
reciprocity_sets.append(reciprocity_3)
reciprocity_sets.append(reciprocity_4)
reciprocity_sets.append(reciprocity_5)
#%% calculating Resistance errors (between normal and reciprocal dataset) 
##########POWER LAW error model##############
filtered_R = []
for i in range(len(reciprocity_sets)):
    reciprocity_sets[i]['avg_R_n/r'] = reciprocity_sets[i][['Rn','Rr']].mean(axis=1) # average Resistance ((Rn+Rr)/2)
    reciprocity_sets[i]['|R_err|'] = np.abs(reciprocity_sets[i]['Rn']-reciprocity_sets[i]['Rr']) # absolute Reciprocal error (|Rn - Rr|)
    reciprocity_sets[i]['percent_error'] = 100*(reciprocity_sets[i]['Rn']-reciprocity_sets[i]['Rr'])/reciprocity_sets[i]['avg_R_n/r'] # error percentage (100*(Rn - Rr)/R_avg)
    ##%% entire data set error histograms:
    error_dist = plt.hist(reciprocity_sets[i]['percent_error'],50) # clearly shows most of data has an error around 0%
    plt.ylabel('Frequency')
    plt.xlabel('error (%)')
    plt.title('error distribution')
    plt.show(error_dist)
    print('clearly shows most of data has an error around "0" percent')
    ## Filtering data <20%
    reciprocity_sets[i] = reciprocity_sets[i][np.abs(reciprocity_sets[i]['percent_error'])<=20] #Filtering based only on R
    filtered_R.append(reciprocity_sets[i])
    ##%% error analysis and error model for data with less than 20% reciprocal error
    numbins = 12 # I prefer 20 bins
    binsize = int(len(reciprocity_sets[i])/numbins) # roughly each been will hav 61 measurements
    error_input = np.abs(reciprocity_sets[i][['avg_R_n/r','|R_err|']]).sort_values(by='avg_R_n/r') # Sorting data based on R_avg
    bins = np.zeros((numbins,2))
    for j in range(numbins): # bining 
        ns=j*binsize
        ne=ns+binsize-1
        bins[j,0] = error_input['avg_R_n/r'].iloc[ns:ne].mean()
        bins[j,1] = error_input['|R_err|'].iloc[ns:ne].mean()    
    coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]))[0] # calculating fitting coefficients (a,m)
    R_error_predict = fit(bins[:,0],np.exp(coefs[0]),coefs[1]) # error prediction based of fitted model
    fig = plt.figure()
    errmodel_plot = plt.loglog(bins[:,0],bins[:,1],'o',label="error")
    fit_line = plt.plot(bins[:,0],R_error_predict,'r', label="fit")
    plt.ylabel('R_error (ohm)')
    plt.xlabel('Average R_n/r (ohm)')
    plt.title('Resistance error model')
    plt.legend(loc='best')
    R2= R_sqr(np.log(bins[:,1]),np.log(R_error_predict))
    RMSE = mse(np.log(bins[:,1]),np.log(R_error_predict))**0.5
    a1 = np.around(np.exp(coefs[0]),decimals=3)
    a2 = np.around(coefs[1], decimals=3)
    a3 = np.around(np.exp(coefs[0]),decimals=1)
    a4 = np.around(coefs[1], decimals=1)
    print ('Error model is: R_err = %s*%s^%s (R^2 = %s) \nor simply R_err = %s*%s^%s' % (a1,'(R_n/r)',a2,R2,a3,'(R_n/r)',a4))
    ax = fig.add_subplot(111)
    ax.text(0.75, 0.1,'R_err = %s*%s^%s \n(R^2 = %s)' % (a1,'(R_n/r)',a2,R2), verticalalignment='center', horizontalalignment='center', transform=ax.transAxes,fontsize=10)
    plt.savefig('%serror_model_%s.png' % (fdir,i),dpi=800)
    plt.show()
    ### Narrowing down the error distribution and probablity within 100% of error
    error_dist_n = plt.hist(reciprocity_sets[i]['percent_error'],bins=(np.arange(-100,100,0.5)),normed=True,alpha=.3,label="Probablity")
    err_5 = reciprocity_sets[i]['percent_error'][np.abs(reciprocity_sets[i]['percent_error'])<=5] # narrowing the fit error for better visualization
    fit_err_dist = mlab.normpdf(np.arange(-100,100,0.5), np.mean(err_5), np.std(err_5)) # parametric gaussian fit (normal pdf)
    KDEpdf = gaussian_kde(err_5) # non parametrict gaussian fit (Kernel density estimation)
    fit_err_dist_plot = plt.plot(np.arange(-100,100,0.5),fit_err_dist,'r--',label="Parametric fit")
    KDE = plt.plot(np.arange(-100,100,0.5),KDEpdf(np.arange(-100,100,0.5)),'k',label="KDE fit")
    plt.xlim(-20,20)
    plt.ylabel('Probablity')
    plt.xlabel('error (%)')
    plt.title('error probablity distribution')
    plt.legend(loc='best')
    plt.savefig('%serror_Histogram_%s.png' % (fdir,i),dpi=800)
    plt.show()
    print('The probablity of data to have ~0 percent error is more than %s \nTherefor we simply can ignore data that has more than 10 percent error and the dataset would stay intact' % (np.around(np.max(error_dist_n[0]), decimals=3)))
#%% entire data error model (power law)
merged_data = pd.concat(filtered_R)
numbins = 12 # I prefer 20 bins
binsize = int(len(merged_data)/numbins) # roughly each been will hav 61 measurements
error_input = np.abs(merged_data[['avg_R_n/r','|R_err|']]).sort_values(by='avg_R_n/r') # Sorting data based on R_avg
bins = np.zeros((numbins,2))
for j in range(numbins): # bining 
    ns=j*binsize
    ne=ns+binsize-1
    bins[j,0] = error_input['avg_R_n/r'].iloc[ns:ne].mean()
    bins[j,1] = error_input['|R_err|'].iloc[ns:ne].mean()    
coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]))[0] # calculating fitting coefficients (a,m)
R_error_predict = fit(bins[:,0],np.exp(coefs[0]),coefs[1]) # error prediction based of fitted model
fig = plt.figure()
errmodel_plot = plt.loglog(bins[:,0],bins[:,1],'o',label="error")
fit_line = plt.plot(bins[:,0],R_error_predict,'r', label="fit")
plt.ylabel('R_error (ohm)')
plt.xlabel('Average R_n/r (ohm)')
plt.title('Resistance error model')
plt.legend(loc='best')
R2= R_sqr(np.log(bins[:,1]),np.log(R_error_predict))
RMSE = mse(np.log(bins[:,1]),np.log(R_error_predict))**0.5
a1 = np.around(np.exp(coefs[0]),decimals=3)
a2 = np.around(coefs[1], decimals=3)
a3 = np.around(np.exp(coefs[0]),decimals=1)
a4 = np.around(coefs[1], decimals=1)
print ('Error model is: R_err = %s*%s^%s (R^2 = %s) \nor simply R_err = %s*%s^%s' % (a1,'(R_n/r)',a2,R2,a3,'(R_n/r)',a4))
ax = fig.add_subplot(111)
ax.text(0.75, 0.1,'R_err = %s*%s^%s \n(R^2 = %s)' % (a1,'(R_n/r)',a2,R2), verticalalignment='center', horizontalalignment='center', transform=ax.transAxes,fontsize=10)
plt.savefig('%serror_model_total.png' % (fdir),dpi=800)
plt.show()
## error histograms:
error_dist_n = plt.hist(merged_data['percent_error'],bins=(np.arange(-100,100,0.5)),normed=True,alpha=.3,label="Probablity")
err_5 = merged_data['percent_error'][np.abs(merged_data['percent_error'])<=5] # narrowing the fit error for better visualization
fit_err_dist = mlab.normpdf(np.arange(-100,100,0.5), np.mean(err_5), np.std(err_5)) # parametric gaussian fit (normal pdf)
KDEpdf = gaussian_kde(err_5) # non parametrict gaussian fit (Kernel density estimation)
fit_err_dist_plot = plt.plot(np.arange(-100,100,0.5),fit_err_dist,'r--',label="Parametric fit")
KDE = plt.plot(np.arange(-100,100,0.5),KDEpdf(np.arange(-100,100,0.5)),'k',label="KDE fit")
plt.xlim(-20,20)
plt.ylabel('Probablity')
plt.xlabel('error (%)')
plt.title('error probablity distribution')
plt.legend(loc='best')
plt.savefig('%serror_Histogram_total.png' % (fdir),dpi=800)
plt.show()
print('The probablity of data to have ~0 percent error is more than %s \nTherefor we simply can ignore data that has more than 10 percent error and the dataset would stay intact' % (np.around(np.max(error_dist_n[0]), decimals=3)))

#%% entire data error model (linear model)
#merged_data = pd.concat(reciprocity_sets)
#merged_data = merged_data.query('(percent_error>-20) & (percent_error<20)')
##merged_data['avg_R_n/r'] = merged_data[['Rn','Rr']].mean(axis=1) # average Resistance ((Rn+Rr)/2)
##merged_data['|R_err|'] = np.abs(merged_data['Rn']-merged_data['Rr']) # absolute Reciprocal error (|Rn - Rr|)
##merged_data['percent_error'] = 100*(merged_data['Rn']-merged_data['Rr'])/merged_data['avg_R_n/r'] # error percentage (100*(Rn - Rr)/R_avg)
###%% error analysis
#numbins = 12 # I prefer 20 bins
#binsize = int(len(merged_data)/numbins) # roughly each been will hav 61 measurements
#error_input = np.abs(merged_data[['avg_R_n/r','|R_err|']]).sort_values(by='avg_R_n/r') # Sorting data based on R_avg
#bins = np.zeros((numbins,2))
#for j in range(numbins): # bining 
#    ns=j*binsize
#    ne=ns+binsize-1
#    bins[j,0] = error_input['avg_R_n/r'].iloc[ns:ne].mean()
#    bins[j,1] = error_input['|R_err|'].iloc[ns:ne].mean()    
#coefs= np.linalg.lstsq(np.vstack([bins[:,0], np.ones(len(bins[:,0]))]).T, bins[:,1])[0] # calculating fitting coefficients (a,m)
#R_error_predict = fit_lin(bins[:,0],coefs[0],coefs[1]) # error prediction based of fitted model
#fig = plt.figure()
#errmodel_plot = plt.plot(bins[:,0],bins[:,1],'o',label="error")
#fit_line = plt.plot(bins[:,0],R_error_predict,'r', label="fit")
#plt.ylabel('R_error (ohm)')
#plt.xlabel('Average R_n/r (ohm)')
#plt.title('Resistance error model')
#plt.legend(loc='best')
#R2= R_sqr(bins[:,1],R_error_predict)
#RMSE = mse(bins[:,1],R_error_predict)**0.5
#a1 = np.around(coefs[0],decimals=3)
#a2 = np.around(coefs[1], decimals=3)
#a3 = np.around(coefs[0],decimals=1)
#a4 = np.around(coefs[1], decimals=1)
#print ('Error model is: R_err = %s*%s+%s (R^2 = %s) \nor simply R_err = %s*%s^%s' % (a1,'(R_n/r)',a2,R2,a3,'(R_n/r)',a4))
##plt.text(0.2,0.002,'R_err = %s*%s^%s \n(R^2 = %s)' % (a3,'(R_n/r)',a4,R2),fontsize=10)
#ax = fig.add_subplot(111)
#ax.text(0.75, 0.1,'R_err = %s*%s+%s \n(R^2 = %s)' % (a1,'(R_n/r)',a2,R2), verticalalignment='center', horizontalalignment='center', transform=ax.transAxes,fontsize=10)
#plt.savefig('%serror_model_total.png' % (fdir),dpi=800)
#plt.show()

#%% finding common measurements among normal data sets
common_normal,nums = commons_normal(filtered_R,normal_sets,17)
print(nums)
n,f=1,0
common_good =[]
for i in range(len(common_normal)):
    if len(common_normal[i])>=600:
        temp_common_good = common_normal[i]
        common_good.append(temp_common_good)
while n!=f:
    common_good,n = commons_normal(common_good,common_good,len(common_good))
    temp = pd.concat(common_good)
    f=len(temp)/len(common_good)
    if f>n-0.2 and f<n+0.2:
        break
#%% Calculating K (geometric factor) and rho
common_good_rho = []
for i in range(len(common_good)):
    df_with_coordinates_0 = pd.merge(common_good[i],coordinates,left_on =['An'],right_on=['electrode'],how='inner').drop(['electrode'],1).rename(columns = {'x':'Ax','y':'Ay','z':'Az'})
    df_with_coordinates_1 = pd.merge(df_with_coordinates_0,coordinates,left_on =['Bn'],right_on=['electrode'],how='inner').drop(['electrode'],1).rename(columns = {'x':'Bx','y':'By','z':'Bz'})
    df_with_coordinates_2 = pd.merge(df_with_coordinates_1,coordinates,left_on =['Mn'],right_on=['electrode'],how='inner').drop(['electrode'],1).rename(columns = {'x':'Mx','y':'My','z':'Mz'})
    df_with_coordinates_final = pd.merge(df_with_coordinates_2,coordinates,left_on =['Nn'],right_on=['electrode'],how='inner').drop(['electrode'],1).rename(columns = {'x':'Nx','y':'Ny','z':'Nz'})
    df_with_coordinates_final['AM'] = np.sqrt(((df_with_coordinates_final['Ax'].values-df_with_coordinates_final['Mx'].values)**2)+((df_with_coordinates_final['Ay'].values-df_with_coordinates_final['My'].values)**2)+((df_with_coordinates_final['Az'].values-df_with_coordinates_final['Mz'].values)**2))
    df_with_coordinates_final['AN'] = np.sqrt(((df_with_coordinates_final['Ax'].values-df_with_coordinates_final['Nx'].values)**2)+((df_with_coordinates_final['Ay'].values-df_with_coordinates_final['Ny'].values)**2)+((df_with_coordinates_final['Az'].values-df_with_coordinates_final['Nz'].values)**2))
    df_with_coordinates_final['BN'] = np.sqrt(((df_with_coordinates_final['Bx'].values-df_with_coordinates_final['Nx'].values)**2)+((df_with_coordinates_final['By'].values-df_with_coordinates_final['Ny'].values)**2)+((df_with_coordinates_final['Bz'].values-df_with_coordinates_final['Nz'].values)**2))
    df_with_coordinates_final['BM'] = np.sqrt(((df_with_coordinates_final['Bx'].values-df_with_coordinates_final['Mx'].values)**2)+((df_with_coordinates_final['By'].values-df_with_coordinates_final['My'].values)**2)+((df_with_coordinates_final['Bz'].values-df_with_coordinates_final['Mz'].values)**2))
    df_with_coordinates_final['K'] = (2*np.pi)/((1/df_with_coordinates_final['AM'].values)-(1/df_with_coordinates_final['AN'].values)-(1/df_with_coordinates_final['BM'].values)+(1/df_with_coordinates_final['BN'].values))
    df_with_coordinates_final['rho'] = (df_with_coordinates_final['V_n'].values*df_with_coordinates_final['K'].values)/df_with_coordinates_final['I_n'].values
    df_with_coordinates_final = df_with_coordinates_final[df_with_coordinates_final['rho']>0]
    data_final = pd.concat((df_with_coordinates_final.iloc[:,:8],df_with_coordinates_final['rho']),axis=1)
    common_good_rho.append(data_final)
#%% finding common measurements among normal data sets (final step after calulting the rho)
n,f=1,0
while n!=f:
    common_good_rho,n = commons_normal(common_good_rho,common_good_rho,len(common_good_rho))
    temp = pd.concat(common_good_rho)
    f=len(temp)/len(common_good_rho)
    if f>n-0.2 and f<n+0.2:
        break
#%% plotting time-lapse changes and building time-lapse ready txt files
bkg  = common_good_rho[0]
for i in range(len(common_good_rho)):
    changes = pd.DataFrame((common_good_rho[i]['rho'].values-bkg['rho'].values)*100/bkg['rho'].values)
    plt.plot(np.arange(len(changes)),changes)
    plt.ylim(-50,3000)
    plt.xlabel('measurement no.')
    plt.ylabel('100*(Rt-R0)/R0')
    plt.savefig('%schanges_%s.png' % (fdir,i),dpi=800)
    plt.show()
    locations_e4d = pd.DataFrame(pd.concat((coordinates[['electrode','x','y']],pd.DataFrame(np.zeros(len(coordinates))).astype(int),pd.DataFrame(np.ones(len(coordinates))).astype(int)),axis=1))
    data_e4d = pd.concat((pd.DataFrame(np.arange(1,len(common_good_rho[i])+1)).rename(columns = {0:'num'}).astype(int),common_good_rho[i][['An','Bn','Mn','Nn']].astype(int),common_good_rho[i]['Rn'],pd.DataFrame(0.001+np.abs(common_good_rho[i]['Rn'])*0.04).rename(columns = {'Rn':'err'})),axis=1)
    save_file = open('%sset_%s.srv' % (fdir,i), 'w')
    save_file.write(str(len(locations_e4d)) + '\n')
    save_file.close()
    with open('%sset_%s.srv' % (fdir,i),mode='a') as file:
        locations_e4d.to_csv(file,sep='\t',header = None,index=None)
        file.write(str(len(data_e4d)) + '\n')
        data_e4d.to_csv(file,sep='\t',header = None,index=None)