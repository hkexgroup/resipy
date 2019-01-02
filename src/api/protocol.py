#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 16:05:59 2018

@author: pmclachlan
"""
import numpy as np
import pandas as pd

#dipole dipole
def dpdp1(elec_num, a, n):
    '''
    genetrates measurement matrix for dipole dipole survey
    elec_num is number of electrodes
    a is electode spacing between C and V pairs (a = 1 is the same as skip 0)
    a should be single integer or list
    n is the quadrupole seperation
    n should be single integer or list
    length of a should match n
    '''
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        n_ = np.array(range(0,n))+1
        for j in np.array(range(0,len(n_))):
            A = elec_id
            B = A + a
            M = B + n_[j] * a
            N = M + a
            if (j) == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if (j) > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
    
    else:
        for i in np.array(range(0,len(a))):
            n_ = np.array(range(0,n[i]))+1
            for j in np.array(range(0,len(n_))):
                A = elec_id
                B = A + a[i]
                M = B + n_[j] * a[i]
                N = M + a[i]
                if (i + j) == 0:
                    proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
                if (i + j) > 0:
                    new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                    proto_mtrx = proto_mtrx.append(new_rows)
    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,3] <= elec_num]
    return proto_mtrx


def dpdp2(elec_num, a, n):
   #genetrates measurement matrix for dipole dipole survey
    #elec_num is number of electrodes
    #a is electode spacing between C and V pairs (a = 1 is the same as skip 0)
    #a should be single integer or list
    #n is the quadrupole seperation
    #n should be single integer or list
    #length of a should match n
    
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        n_ = np.array(range(0,n))+1
        for j in np.array(range(0,len(n_))):
            A = elec_id
            B = A + a
            M = B + n_[j]
            N = M + a
            if (j) == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if (j) > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
    
    else:
        for i in np.array(range(0,len(a))):
            n_ = np.array(range(0,n[i]))+1
            for j in np.array(range(0,len(n_))):
                A = elec_id
                B = A + a[i]
                M = B + n_[j]
                N = M + a[i]
                if (i + j) == 0:
                    proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
                if (i + j) > 0:
                    new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                    proto_mtrx = proto_mtrx.append(new_rows)
    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,3] <= elec_num]
    return proto_mtrx

def wenner_alpha(elec_num, a):
    #genetrates measurement matrix for dipole dipole survey
    #elec_num is number of electrodes
    #a is electode spacing
    #a can be list or int
    
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        A = elec_id
        M = A + a
        N = M + a
        B = N + a
        proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
    
    else:
        for i in np.array(range(0,len(a))):
            A = elec_id
            M = A + a[i]
            N = M + a[i]
            B = N + a[i]
            if (i) == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if (i) > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
                    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] <= elec_num]
    return(proto_mtrx)
    
    
def wenner_beta(elec_num, a):
    '''
    genetrates measurement matrix for dipole dipole survey
    elec_num is number of electrodes
    a is electode spacing
    a can be list or int
    '''
    
    elec_id = np.arange(elec_num)+1
    
    if type(a) == int:
        B = elec_id
        A = B + a
        M = A + a
        N = M + a
        proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
    
    if type(a) == list:
        for i in np.array(range(0,len(a))):
            B = elec_id
            A = B + a[i]
            M = A + a[i]
            N = M + a[i]
            if (i) == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if (i) > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
                    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] <= elec_num]
    return(proto_mtrx)
    
def wenner_gamma(elec_num, a):
    '''
    genetrates measurement matrix for dipole dipole survey
    elec_num is number of electrodes
    a is electode spacing
    a can be list or int
    '''
    
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        A = elec_id
        M = A + a
        B = M + a
        N = B + a
        proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
    
    else:
        for i in np.array(range(0,len(a))):
            A = elec_id
            M = A + a[i]
            B = M + a[i]
            N = B + a[i]
            if (i) == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if (i) > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
                    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] <= elec_num]
    return(proto_mtrx)
    
def schlum1(elec_num, a, n):
    #genetrates measurement matrix for dipole dipole survey
    #elec_num is number of electrodes
    #a is electode spacing
    
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        n_ = np.array(range(0,n))+1
        for j in np.array(range(0,len(n_))):
            A = elec_id
            M = A + n_[j] * a
            N = M + a
            B = N + n_[j] * a
            if j == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if j > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
                
    else:
        for i in np.array(range(0,len(a))):
            n_ = np.array(range(0,n[i]))+1
            for j in np.array(range(0,len(n_))):
                A = elec_id
                M = A + n_[j] * a[i]
                N = M + a[i]
                B = N + n_[j] * a[i]
                if (i + j) == 0:
                    proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
                if (i + j) > 0:
                    new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                    proto_mtrx = proto_mtrx.append(new_rows)
    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] <= elec_num]
    return(proto_mtrx)
    
    
def schlum2(elec_num, a, n):
    #genetrates measurement matrix for dipole dipole survey
    #elec_num is number of electrodes
    #a is electode spacing
    
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        n_ = np.array(range(0,n))+1
        for j in np.array(range(0,len(n_))):
            A = elec_id
            M = A + n_[j] * a
            N = M + a
            B = N + n_[j] * a
            if j == 0:
                proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
            if j > 0:
                new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                proto_mtrx = proto_mtrx.append(new_rows)
                
    else:
        for i in np.array(range(0,len(a))):
            n_ = np.array(range(0,n[i]))+1
            for j in np.array(range(0,len(n_))):
                A = elec_id
                M = A + n_[j] 
                N = M + a[i]
                B = N + n_[j] 
                if (i + j) == 0:
                    proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
                if (i + j) > 0:
                    new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                    proto_mtrx = proto_mtrx.append(new_rows)
    
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] <= elec_num]
    return proto_mtrx
    
    
def multigrad(elec_num, a, n, s):
    #genetrates measurement matrix for multigradient array Torleif Dahlin
    #a is spacing between potential electrodes
    #n is multiplier for a to determine spacing from A to M
    #s is seperation factor for current electrodes, should be the intermediate numbers 
    
    elec_id = np.arange(elec_num)+1
    
    if isinstance(a, list) is False:
        n_ = np.array(range(0, n))+1
        s_ = np.array(range(0, s))+1
        for j in np.array(range(0, len(n_))):
            for k in np.array(range(0, len(s_))):
                A = elec_id
                B = A + s_[k] + 2
                M = A + n_[j] * a
                N = M + a
                if (j + k) == 0:
                    proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
                if (j + k) > 0:
                    new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                    proto_mtrx = proto_mtrx.append(new_rows)
                    
    else:
        for i in np.array(range(0,len(a))):
            n_ = np.array(range(0, n[i]))+1
            s_ = np.array(range(0, s[i]))+1
            for j in np.array(range(0, len(n_))):
                for k in np.array(range(0, len(s_))):
                    A = elec_id
                    B = A + s_[k] + 2
                    M = A + n_[j] * a[i]
                    N = M + a[i]
                    if (i + j + k) == 0:
                        proto_mtrx = pd.DataFrame(np.column_stack((A, B, M, N)))
                    if (i + j + k) > 0:
                        new_rows = pd.DataFrame(np.column_stack((A, B, M, N)))
                        proto_mtrx = proto_mtrx.append(new_rows)
                
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] > proto_mtrx.iloc[:,3]]          
    proto_mtrx = proto_mtrx[proto_mtrx.iloc[:,1] <= elec_num]
    return proto_mtrx    

# test code
#x1 = dpdp1(24, 2, 8)
#x2 = dpdp2(24, 2, 8)
#x3 = wenner_alpha(24, 1)
#x4 = wenner_beta(24, 1)
#x5 = wenner_gamma(24, 1)
#x6 = schlum1(24, 1, 10)
#x7 = schlum2(24, 1, 10)
#x8 = multigrad(24, 1, 10, 2)



#%% show pseudoSection
#import matplotlib.pyplot as plt    
#
#N = 24
#elecpos = np.linspace(0, 8, N)
#quad = pd.concat([wenner_alpha(N, 1), wenner_alpha(N, 2)]).values
#array = np.sort(quad, axis=1)
#
#cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
#    + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
#pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
#    + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
#xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
#ypos = - np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
#
#
#fig, ax = plt.subplots()
#cax = ax.plot(xpos, ypos, 'o')
#ax.set_title('Pseudo Section')
#fig.show()