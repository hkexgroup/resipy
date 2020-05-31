# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 09:17:45 2018

@author: Sina
"""

import numpy as np
import pandas as pd
# import warnings
# warnings.filterwarnings("ignore")
#%%
def positive_test (Dcurve,DecayTime): 
    """Calculating TDIP chargeability decay curve trend: 
        positive (increasing over time) trends are bad data
    """
    DC_slope = np.zeros(np.shape(Dcurve)[0])
    for i in range(np.shape(Dcurve)[0]):     
        DecayCoefs = np.linalg.lstsq(np.vstack([DecayTime, np.ones(len(DecayTime))]).T,Dcurve.iloc[i,:])[0]
        DC_slope[i]=DecayCoefs[0]
    return DC_slope

def linear_coefs (x,y): #linear fit parameteres for decay curve
    data = np.concatenate((np.log(x)[:,None],np.log(y)[:,None]),axis=1)
    if np.log(y).sum() !=0:
        data_no_nan = data[~np.isnan(data).any(axis=1)]
        coefs = np.linalg.lstsq(np.vstack([data_no_nan[:,0], np.ones(len(data_no_nan[:,0]))]).T,data_no_nan[:,1])[0]
        return coefs
    else:
        coefs = np.array([0,0])
        return

def DCA(data_in, dump=None): 
    """Decay Curve Analysis (Only for Syscal files):
        calculating master decay curve based on individual decay curves, 
        then compares individual decay curves with a master decay curve (avg(all good curves)) 
        and remove data with STD > 2 * STD(dataset).
        
    Reference:
    ----------
    Flores Orozco, A., Gallistl, J., Bücker, M., & Williams, K. H. (2017)., 
    Decay curve analysis for data error quantification in time-domain induced polarization imaging., 
    Geophysics, 83(2), 1–48. https://doi.org/10.1190/geo2016-0714.1)
    """
    if dump is None:
        def dump(x):
            pass
    data = data_in.copy()
    decayN = data[['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
                'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']]
    DecayTime_int = data['TM1'][0]
    DecayTime = np.arange(int(DecayTime_int),np.shape(decayN)[1]*(int(DecayTime_int)+1),int(DecayTime_int))
    data['DC_slope'] = positive_test(decayN.copy(),DecayTime) #decay curve trend - positive trends are bad data
    
    if data['ip'].mean() == 0: 
        print('\nNo reciprocal IP data available (fast reciprocal measurement)')
        filtered_R_IP = data.query('(DC_slope<0)').rename(columns = {'a':'An', 'b':'Bn'})
    else:
        filtered_R_IP = data.query('(DC_slope<0)').rename(columns = {'a':'An', 'b':'Bn'})
        
    DC_fit_linear = np.zeros([len(filtered_R_IP['ip']),2])
    for i in range(len(filtered_R_IP['ip'])): #calculating decay curve fit parameteres - m=at^b (m: chargeability, t: time, a and b: fitting parameters)
        DC_fit_linear[i,:] = linear_coefs(DecayTime,filtered_R_IP[['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
                    'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']].iloc[i,:])            
    filtered_R_IP['DC_fit_a'],filtered_R_IP['DC_fit_b']=np.exp(DC_fit_linear[:,1]),DC_fit_linear[:,0]
    filtered_R_IP = filtered_R_IP.dropna(subset=['DC_fit_a','DC_fit_b']).query('(DC_fit_a>0) & (DC_fit_b<0) | (DC_fit_a<0) & (DC_fit_b>0)') #filtering meaningless decay curves 
    filtered_R_IP = filtered_R_IP.dropna(subset=['DC_fit_a','DC_fit_b']) #for test, remove it
    fit_DC = pd.DataFrame(filtered_R_IP['DC_fit_a'][:,None]*DecayTime**filtered_R_IP['DC_fit_b'][:,None]).rename(columns = {0:'Fit_DC_1', 1:'Fit_DC_2', 2:'Fit_DC_3', 3:'Fit_DC_4',4:'Fit_DC_5', 5:'Fit_DC_6', 6:'Fit_DC_7', 7:'Fit_DC_8',
                        8:'Fit_DC_9', 9:'Fit_DC_10',10:'Fit_DC_11',11:'Fit_DC_12',12:'Fit_DC_13',13:'Fit_DC_14',14:'Fit_DC_15',15:'Fit_DC_16',
                        16:'Fit_DC_17',17:'Fit_DC_18',18:'Fit_DC_19',19:'Fit_DC_20'}) #building fitted decay curve       
    filtered_R_IP = pd.concat([filtered_R_IP.reset_index(drop=True),fit_DC.reset_index(drop=True)], axis=1) #conatenating fitted decay curve dataframe with measured decaycurve        
    filtered_R_IP['DC_rmsd'] = ((((fit_DC.values-filtered_R_IP[['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
                'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']].values)**2).sum(axis=1))/np.shape(fit_DC)[1])**0.5 #calculating decay curve RMSD with fitted curve
    ####  Building master decay curves for each IP window
    filtered_R_IP['weight'] = 1/filtered_R_IP['DC_rmsd']
    temp_DC_master = pd.concat([filtered_R_IP[['An','Bn']],(pd.DataFrame(fit_DC.values*filtered_R_IP['weight'][:,None]))],axis=1)
    master_DC = pd.DataFrame((temp_DC_master.groupby(['An','Bn']).sum()).values/(filtered_R_IP.groupby(['An','Bn'])['weight']).sum()[:,None])
    groups_temp = filtered_R_IP.groupby(['An','Bn'],sort=True)
    groups_keys = []
    for key in groups_temp.groups.keys():
        groups_keys.append(key)
    groups_keys = pd.DataFrame(np.reshape(groups_keys, (len(groups_keys),2))).rename(columns = {0:'An', 1:'Bn'})
    groups_keys = pd.concat([groups_keys,pd.DataFrame(np.arange(0,len(groups_keys)))], axis=1).rename(columns = {0:'i'})
    groups_keys_sorted = groups_keys.sort_values(['An','Bn']).reset_index(drop=True)
    master_DC_indexed = pd.concat([master_DC,groups_keys_sorted['i']], axis=1).sort_values(['i']).reset_index(drop=True)
    shift = np.linspace(-10, 10, 20)
    appended_groups = []
    groups_list = [groups_temp.get_group(x) for x in groups_temp.groups]  
    for i in range(len(filtered_R_IP.groupby(['An','Bn']))):
        group = groups_list[i]
        temp_K=np.zeros([len(group)])
        for j in range(len(group)):
            errors = np.array([np.mean((master_DC_indexed.iloc[i,:-1][:,None] - group[['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
                                        'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']].iloc[j,:][:,None] + m)) for m in shift])
            temp_K[j] = np.polyfit(shift, errors, 2)[2]
        # group.loc[:, 'K_std'] = temp_K
        # appended_groups.append(group)
        dfK = pd.DataFrame(temp_K).rename(columns={0:'K_std'})
        appended_groups.append(pd.concat([group.reset_index(drop=True),dfK.reset_index(drop=True)], axis=1))
        percent_progress = i*100/len(filtered_R_IP.groupby(['An','Bn']))
        dump(percent_progress)
        print('\r%s%s -Done' % (int(percent_progress),'%'), end='')
    dump(100)
    print('\r100% -Done - finished!')
    appended_groups = pd.concat(appended_groups)
    K_std = np.std(appended_groups['K_std'])
    final_data = appended_groups.loc[np.abs(appended_groups['K_std'])<(2*K_std)]
    final_data_dropped = final_data.drop(['DC_slope', 'DC_fit_a', 'DC_fit_b', 'Fit_DC_1', 
                                   'Fit_DC_2', 'Fit_DC_3', 'Fit_DC_4', 'Fit_DC_5', 'Fit_DC_6', 'Fit_DC_7',
                                   'Fit_DC_8', 'Fit_DC_9', 'Fit_DC_10', 'Fit_DC_11', 'Fit_DC_12',
                                   'Fit_DC_13', 'Fit_DC_14', 'Fit_DC_15', 'Fit_DC_16', 'Fit_DC_17',
                                   'Fit_DC_18', 'Fit_DC_19', 'Fit_DC_20', 'DC_rmsd', 'weight', 'K_std'], axis = 1).rename(columns = {'An':'a','Bn':'b'})
    return (final_data_dropped)