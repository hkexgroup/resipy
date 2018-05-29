import numpy as np
import os
import pandas as pd
from scipy.stats.kde import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
# Finding common measurements between Normal and Reciprocal datasets and doing the reciprocal error analysis (power law)
# Dataset structure: Pandas, Normal dataset columns: An, Bn, Mn, Nn, Rn and Reciprocal dataset columns: Ar, Br, Mr, Nr, Rr
# needed functions below
def fit(x,a,m): # power-law fit (based on previous studies)
    return a*(x**m)
def R_sqr (y,y_predict): # calculating R squared value to measure fitting accuracy
    	rsdl = y - y_predict
    	ss_res = np.sum(rsdl**2)
    	ss_tot = np.sum((y-np.mean(y))**2)
    	R2 = 1-(ss_res/ss_tot)
    	R2 = np.around(R2,decimals=4)
def linear_coefs (x,y): #linear fit parameteres for decay curve
    data = np.concatenate((np.log(x)[:,None],np.log(y)[:,None]),axis=1)
    if np.log(y).sum() !=0:
        data_no_nan = data[~np.isnan(data).any(axis=1)]
        coefs = np.linalg.lstsq(np.vstack([data_no_nan[:,0], np.ones(len(data_no_nan[:,0]))]).T,data_no_nan[:,1])[0]
        return coefs
    else:
        coefs = np.array([0,0])
        return
def find_common (N,R):
	df_1 = N
	df_2 = R
	commons_1 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Mr', 'Nr', 'Ar', 'Br'])
	commons_2 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Nr', 'Mr', 'Ar', 'Br'])
	commons_2['Rr'] = commons_2['Rr']*-1
	commons_3 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Mr', 'Nr', 'Br', 'Ar'])
	commons_3['Rr'] = commons_3['Rr']*-1
	commons_4 = pd.merge(df_1, df_2, how = 'inner', left_on = ['An', 'Bn', 'Mn', 'Nn'], right_on = ['Nr', 'Mr', 'Br', 'Ar'])
	commons = pd.concat([commons_1, commons_2, commons_3, commons_4])
	return commons

#%% calculating Resistance errors (between normal and reciprocal dataset) (Reference: Robinson, J. L., L. D. Slater, and K. V. R. Schäfer (2012), Evidence for spatial variability in hydraulic redistribution within an oak-pine forest from resistivity imaging, J. Hydrol., 430–431, 69–79, doi:10.1016/j.jhydrol.2012.02.002.)
common = find_common(dataN,dataR)
common = common.replace([np.inf, -np.inf], np.nan).dropna(subset=['Rn','Rr'])
common['avg_R_n/r'] = common[['Rn','Rr']].mean(axis=1) # average Resistance ((Rn+Rr)/2)
common['|R_err|'] = np.abs(common['Rn']-common['Rr']) # absolute Reciprocal error (|Rn - Rr|)
common['percent_error'] = 100*(common['Rn']-common['Rr'])/common['avg_R_n/r'] # error percentage (100*(Rn - Rr)/R_avg)
common = common.dropna(subset=['percent_error']).drop_duplicates(subset=['An', 'Bn', 'Mn', 'Nn'], keep = 'first')
#%% error histograms:
error_dist = plt.hist(common['percent_error'],50) # clearly shows most of data has an error around 0%
plt.ylabel('Frequency')
plt.xlabel('error (%)')
plt.title('error distribution')
plt.show(error_dist)
plt.gcf().clear()
print('clearly shows most of data has an error around "0" percent')
#%% error analysis (power law)
numbins = 20 # I prefer 20 bins
temp_filtered = common[np.abs(common['percent_error'])<=int(Arrayinfo[5])]
binsize = int(len(temp_filtered)/numbins) 
error_input = np.abs(temp_filtered[['avg_R_n/r','|R_err|']]).sort_values(by='avg_R_n/r') # Sorting data based on R_avg
bins = np.zeros((numbins,2))
for i in range(numbins): # bining 
	ns=i*binsize
	ne=ns+binsize-1
	bins[i,0] = error_input['avg_R_n/r'].iloc[ns:ne].mean()
	bins[i,1] = error_input['|R_err|'].iloc[ns:ne].mean()    
coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]))[0] # calculating fitting coefficients (a,m)
R_error_predict = fit(bins[:,0],np.exp(coefs[0]),coefs[1]) # error prediction based of fitted model
fig = plt.figure(figsize=(8,6))
errmodel_plot = plt.loglog(bins[:,0],bins[:,1],'o',label="error", markersize = 15)
fit_line = plt.plot(bins[:,0],R_error_predict,'r', label="fit", linewidth = 5)
plt.ylabel(r'Log$R_{error} [\Omega]$', fontsize=28)
plt.xlabel(r'Log$R_{avg} [\Omega]$', fontsize=28)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
#    plt.title('Resistance error model')
#    plt.legend(loc='best',fontsize = 22)
leg = plt.legend(loc='best', fontsize =24, frameon=True)
leg.set_frame_on(True)#.set_edgecolor(0.5)
leg.get_frame().set_edgecolor("black")
leg.get_frame().set_linewidth(2)
R2= R_sqr(np.log(bins[:,1]),np.log(R_error_predict))
RMSE = mse(np.log(bins[:,1]),np.log(R_error_predict))**0.5
a1 = np.around(np.exp(coefs[0]),decimals=3)
a2 = np.around(coefs[1], decimals=3)
a3 = np.around(np.exp(coefs[0]),decimals=1)
a4 = np.around(coefs[1], decimals=1)
print ('Error model is: R_err = %s*%s^%s (R^2 = %s) \nor simply R_err = %s*%s^%s' % (a1,'(R_n/r)',a2,R2,a3,'(R_n/r)',a4))
ax = fig.add_subplot(111)
plt.title(r'$\alpha = %s, \beta = %s$' % (a1,a2), fontsize = 28) # error parameters will be ritten on top
# error parameters will be written as in-set equation vvvvvv
#    ax.text(0.75, 0.1,'R_err = %s*%s^%s \n(R^2 = %s)' % (a1,'(R_n/r)',a2,R2), verticalalignment='center', horizontalalignment='center', transform=ax.transAxes,fontsize=18)
ax.tick_params(axis='both', labelsize=22)
ax.grid(False)
#    plt.text(0.2,0.002,'R_err = %s*%s^%s \n(R^2 = %s)' % (a3,'(R_n/r)',a4,R2),fontsize=10)
plt.savefig('%serror_model.png' % (fdir),dpi=800)
plt.show()
plt.gcf().clear()
%## Narrowing down the error distribution and probablity within 100% of error (better visualization)
error_dist_n = plt.hist(common['percent_error'],bins=(np.arange(-100,100,0.5)),normed=True,alpha=.3,label="Probablity")
err_5 = common['percent_error'][np.abs(common['percent_error'])<=5] # narrowing the fit error for better visualization
fit_err_dist = mlab.normpdf(np.arange(-100,100,0.5), np.mean(err_5), np.std(err_5)) # parametric gaussian fit (normal pdf)
KDEpdf = gaussian_kde(err_5) # non parametrict gaussian fit (Kernel density estimation)
fit_err_dist_plot = plt.plot(np.arange(-100,100,0.5),fit_err_dist,'r--',label="Parametric fit")
KDE = plt.plot(np.arange(-100,100,0.5),KDEpdf(np.arange(-100,100,0.5)),'k',label="KDE fit")
plt.xlim(-1*(int(Arrayinfo[5])),int(Arrayinfo[5]))
plt.ylabel('Probablity')
plt.xlabel('error (%)')
plt.title('error probablity distribution')
plt.legend(loc='best')
plt.savefig('%serror_Histogram.png' % (fdir),dpi=800)
plt.show()
print('The probablity of data to have ~0 percent error is more than %s \nTherefor we simply can ignore data that has more than 10 percent error and the dataset would stay intact' % (np.around(np.max(error_dist_n[0]), decimals=3)))
