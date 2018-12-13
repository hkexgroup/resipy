#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 13:12:37 2018
Topography interpolation schemes for pyr2 
@author: jimmy
"""
import numpy as np
from api.isinpolygon import isinpolygon
from tqdm import tqdm

#%% compute thin plate spline /bilinear models  for irregular grid
# see solution @ https://math.stackexchange.com/questions/828392/spatial-interpolation-for-irregular-grid
def thin_plate_spline_mod(x,y,z):
    """
    Returns the thin plate spline model. 
    """
    if len(x) != len(y) and len(y)!=len(x):
        raise ValueError("Mismatch in the number of elements in either x, y, or z arrays")
    #construct inverse operator G 
    dimn = len(x)
    G = np.matrix(np.concatenate((x**2,x*y,y**2,x,y,np.ones((dimn,1))),axis=1))
    #compute model
    mod = ((G.T * G)**-1) * G.T * z 
    return mod.A

def bilinear_mod(x,y,z):
    """
    Returns the bilinear model. 
    """
    if len(x) != len(y) and len(y)!=len(x):
        raise ValueError("Mismatch in the number of elements in either x, y, or z arrays")
    #construct inverse operator G 
    dimn = len(x)
    G = np.matrix(np.concatenate((x*y,x,y,np.ones((dimn,1))),axis=1)) 
    #compute model
    mod = ((G.T * G)**-1) * G.T * z 
    return mod.A

def compute(mod,xnew,ynew):
    """
    Compute the znew at xnew and ynew given a model  
    """
    if len(mod)==4:
        znew = mod[0]*xnew*ynew + mod[1]*xnew + mod[2]*ynew + mod[3]
    elif len(mod)==6:
        znew = mod[0]*xnew**2 + mod[1]*xnew*ynew + mod[2]*ynew**2 + mod[3]*xnew + mod[4]*ynew + mod[5]
    else:
        raise ValueError("length of model vector is unrecognised")
    return znew

#%% compute distances between points (2D)
def pdist(x1, y1, x2, y2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2)

def cdist(x1, y1, x2, y2):
    """
    compute distances to other points, compare set 1 to set 2
    """
    dist = np.zeros((len(x1), len(x2)))
    for i, (x, y) in enumerate(zip(x1,y1)):
        dist[i,:] = pdist(x, y, x2, y2)
    return dist
        
#    
#def cdist(x_query,y_query,x0,y0):
    
    
#%% start 
def irregular_grid(xnew, ynew, x_grid, y_grid, z_grid, method="bilinear", extrapolate=True):
    """
    Compute z values for an irregular grid
    """
    #get grid shape 
    dimns = np.shape(x_grid)
    #preallocate array for new z coordinates / interpolated values  
    znew = np.zeros_like(xnew)
    znew.fill(np.nan)
    
    #compute new values inside survey
    for i in tqdm(range(dimns[0]-1),ncols=100,desc="Interpolating values"): # note this done inside 2 for loops - yikes!!! 
        for j in range(dimns[1]-1): #### TODO : make this more efficient 
            x = np.array((x_grid[i,j],x_grid[i,j+1],x_grid[i+1,j+1],x_grid[i+1,j]))
            y = np.array((y_grid[i,j],y_grid[i,j+1],y_grid[i+1,j+1],y_grid[i+1,j]))
            z = np.array((z_grid[i,j],z_grid[i,j+1],z_grid[i+1,j+1],z_grid[i+1,j]))
            if len(x.shape):#bug fix to deal with numpy being a finicky twat 
                x.shape += (1,)
                z.shape += (1,)#append 1 dimension to the numpy array shape (otherwise np.concentrate wont work)
                y.shape += (1,)
            if method=="bilinear":
                mod = bilinear_mod(x,y,z)
            elif method=="spline":
                mod = thin_plate_spline_mod(x,y,z)
            else:
                raise Exception("unrecognised method")
            inside = isinpolygon(xnew,ynew,(x,y))
            znew[inside] = compute(mod,xnew[inside],ynew[inside])
    
    idx_nan = np.isnan(znew) # boolian indexes of where nans are
    idx_num = np.where(idx_nan == False)
    #extrapolate nans using nearest nieghbough interpolation
    if extrapolate:
        #combine known and interpolated values 
        known_x = np.append(x_grid.flatten(),xnew[idx_num])
        known_y = np.append(y_grid.flatten(),ynew[idx_num])
        known_z = np.append(z_grid.flatten(),znew[idx_num])
        extrap_x = xnew[idx_nan] # extrapolate using the gridded and interpolated data
        extrap_y = ynew[idx_nan]
        extrap_z = znew[idx_nan]
        
        for i in tqdm(range(len(extrap_x)),ncols=100,desc="Extrapolating values"):#go through each extrapolated point and find the closest known coordinate
            dist = pdist(extrap_x[i],extrap_y[i],known_x,known_y)
            ref = np.argmin(dist)
            extrap_z[i] = known_z[ref]
        
        znew[idx_nan] = extrap_z
        
    return znew # return new interpolated values 
    
    
    