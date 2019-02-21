#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 13:12:37 2018
Topography interpolation schemes for pyr2 
@author: jimmy
"""
import numpy as np
#from tqdm import tqdm # progess bar package - disabled 

#%% compute thin plate spline /bilinear models  for irregular grid
# see solution @ https://math.stackexchange.com/questions/828392/spatial-interpolation-for-irregular-grid
def thin_plate_spline_mod(x,y,z): # apparently unstable 
    """
    Returns the thin plate spline model. Unstable function. 
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
    """
    Vectorised distance computation between 2 sets of points, 1 and 2. 
    """
    return np.sqrt((x1-x2)**2+(y1-y2)**2)

def cdist(x1, y1, x2, y2):
    """
    compute distances to other points, compare set 1 to set 2
    """
    dist = np.zeros((len(x1), len(x2)))
    for i, (x, y) in enumerate(zip(x1,y1)):
        dist[i,:] = pdist(x, y, x2, y2)
    return dist
        

#%% bilinear interpolation - use closest 4 points 
def bilinear(xnew, ynew, xknown, yknown, zknown, extrapolate=True):
    """
    Compute z values for unstructured data using bilinear interpolation. Coordinates
    outside the bounds of interpolation can be extrapolated using nearest neighbour
    algorithm on the interpolated and known z coordinates.  
    
    Bilinear interpolation requires knowledge of 4 points orientated around the 
    xyz coordinate to be estimated. 
    
    Parameters
    ------------
    xnew: array like
        x coordinates for the interpolated values
    ynew: array like 
        y coordinates for the interpolated values
    xknown: array like
        x coordinates for the known values 
    yknown: array like
        y coordinates for the known values 
    zknown: array like
        z coordinates for the known values 
    extrapolate: bool, optional
        Flag for if extrapolation is to be used if new coordinates lie outside 
        the bounds where it is not possible to interpolate a value. 
        
    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew.
        
    """
    #preallocate array for new z coordinates / interpolated values  
    znew = np.zeros_like(xnew)
    znew.fill(np.nan)
    #outside = np.logical_not(inside)
    num_pts = len(xnew)
    fudgex=0.05#add a bit of padding to prevent artefacts?
    fudgey=0.05
    #compute new values inside survey
    for i in range(num_pts):
    #for i in tqdm(range(num_pts),desc='interpolating values',ncols=100):
        #find closest 4 points in each quad
        quad1 = (xknown < xnew[i]-fudgex) & (yknown < ynew[i]-fudgey) # bottom left quad
        quad2 = (xknown < xnew[i]-fudgex) & (yknown > ynew[i]+fudgey) # top left quad
        quad3 = (xknown > xnew[i]+fudgex) & (yknown > ynew[i]+fudgey) # top right quad
        quad4 = (xknown > xnew[i]+fudgex) & (yknown < ynew[i]-fudgey) # bottom right quad
        
        dist1 = pdist(xnew[i], ynew[i], xknown[quad1], yknown[quad1])#distances to each quad 
        dist2 = pdist(xnew[i], ynew[i], xknown[quad2], yknown[quad2])
        dist3 = pdist(xnew[i], ynew[i], xknown[quad3], yknown[quad3])
        dist4 = pdist(xnew[i], ynew[i], xknown[quad4], yknown[quad4])
        
        if len(dist1)!=0 and len(dist2)!=0 and len(dist3)!=0 and len(dist4)!=0:
            #then the conditions need to interpolate in a quad are met
            idx1 = np.argmin(dist1)#find closest index for each quad 
            idx2 = np.argmin(dist2)
            idx3 = np.argmin(dist3)
            idx4 = np.argmin(dist4)
            
            x = np.array((xknown[quad1][idx1],
                          xknown[quad2][idx2],
                          xknown[quad3][idx3],
                          xknown[quad4][idx4]))
            y = np.array((yknown[quad1][idx1],
                          yknown[quad2][idx2],
                          yknown[quad3][idx3],
                          yknown[quad4][idx4]))
            z = np.array((zknown[quad1][idx1],
                          zknown[quad2][idx2],
                          zknown[quad3][idx3],
                          zknown[quad4][idx4]))
            if len(x.shape)==1:#bug fix to deal with numpy being a finicky twat 
                x.shape += (1,)
                z.shape += (1,)#append 1 dimension to the numpy array shape (otherwise np.concentrate wont work)
                y.shape += (1,)
           
            mod = bilinear_mod(x,y,z) # generate model     
            znew[i] = compute(mod,xnew[i],ynew[i])#interpolate point using model 
        #else: znew is nan (already assigned)
        
    idx_nan = np.isnan(znew) # boolian indexes of where nans are
    idx_num = np.where(idx_nan == False)
    #extrapolate nans using nearest nieghbough interpolation
    if extrapolate:
        #combine known and interpolated values 
        known_x = np.append(xknown,xnew[idx_num])
        known_y = np.append(yknown,ynew[idx_num])
        known_z = np.append(zknown,znew[idx_num])
        extrap_x = xnew[idx_nan] # extrapolate using the gridded and interpolated data
        extrap_y = ynew[idx_nan]
        extrap_z = znew[idx_nan]
        
        for i in range(len(extrap_x)):
        #for i in tqdm(range(len(extrap_x)),desc='extrapolating unknowns',ncols=100):#go through each extrapolated point and find the closest known coordinate
            dist = pdist(extrap_x[i],extrap_y[i],known_x,known_y)
            ref = np.argmin(dist)
            extrap_z[i] = known_z[ref]
        
        znew[idx_nan] = extrap_z
        
    return znew # return new interpolated values 
    
#%% inverse weighted distance
def idw(xnew, ynew, xknown, yknown, zknown, power=2, radius = 10000, extrapolate=True):
    """
    Compute z values for unstructured data using inverse distance weighting. Coordinates
    outside the bounds of interpolation can be extrapolated using nearest neighbour
    algorithm on the interpolated and known z coordinates.  
    
    Parameters
    ------------
    xnew: array like
        x coordinates for the interpolated values
    ynew: array like 
        y coordinates for the interpolated values
    xknown: array like
        x coordinates for the known values 
    yknown: array like
        y coordinates for the known values 
    zknown: array like
        z coordinates for the known values 
    power: float, optional
        Power the weighting function is raised to. 
    raduis: float, optional,
        Search raduis in which points will be selected for IDW interpolation. 
        By default all points in a 10000 unit raduis are selected. 
    extrapolate: bool, optional
        Flag for if extrapolation is to be used if new coordinates lie outside 
        the bounds where it is not possible to interpolate a value. 
        
    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew. 
        
    """
    znew = np.zeros_like(xnew)
    znew.fill(np.nan)
    for i,(x,y) in enumerate(zip(xnew, ynew)):
    #for i,(x,y) in tqdm(enumerate(zip(xnew, ynew)),ncols=100,desc="Interpolating topo"):
        dist = pdist(x, y, xknown, yknown)
        search = dist<=radius #get boolian array where dist is smaller than search radius 
        w = (1/dist)**power # exponent to be chosen
        znew[i] = np.sum(zknown[search]*w[search])/np.sum(w[search])
        
    idx_nan = np.isnan(znew) # boolian indexes of where nans are
    idx_num = np.where(idx_nan == False)
    #extrapolate nans using nearest nieghbough interpolation
    if extrapolate:
        #combine known and interpolated values 
        known_x = np.append(xknown,xnew[idx_num])
        known_y = np.append(yknown,ynew[idx_num])
        known_z = np.append(zknown,znew[idx_num])
        extrap_x = xnew[idx_nan] # extrapolate using the gridded and interpolated data
        extrap_y = ynew[idx_nan]
        extrap_z = znew[idx_nan]
        
        for i in range(len(extrap_x)):#,ncols=100,desc="Extrapolating values"):#go through each extrapolated point and find the closest known coordinate
            dist = pdist(extrap_x[i],extrap_y[i],known_x,known_y)
            ref = np.argmin(dist)
            extrap_z[i] = known_z[ref]     
        znew[idx_nan] = extrap_z
        
    return znew

#%% pure nearest neighbour interpolation
def nearest(xnew, ynew, xknown, yknown, zknown, maxDist=None): 
    """
    Compute z values for unstructured data using nearest neighbour extrapolation.
    Suitable where dense known coordinates occur.ie. in the case of a DEM.  
    
    Parameters
    ------------
    xnew: array like
        x coordinates for the interpolated values
    ynew: array like 
        y coordinates for the interpolated values
    xknown: array like
        x coordinates for the known values 
    yknown: array like
        y coordinates for the known values 
    zknown: array like
        z coordinates for the known values 
    maxDist : float, optional
        Maximum distance for nearest neighbour interpolation. If None then this
        argument is ignored. If given a value then the values outiside the 
        maximum distance will be returned as NaN. 

    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew.  
        
    """
    znew = np.zeros_like(xnew)
    znew.fill(np.nan)        
    if maxDist is None:
        check=False
    else:
        if not isinstance(maxDist,float):# or not isinstance(maxDist,int):
            raise ValueError("maxDist argument must be of type int or float, not %s"%type(maxDist))
        check=True
    for i in range(len(xnew)):#,ncols=100,desc="Extrapolating values"):#go through each extrapolated point and find the closest known coordinate
        dist = pdist(xnew[i],ynew[i],xknown,yknown)
        ref = np.argmin(dist)
        znew[i] = zknown[ref] 
        if check:
            if dist[ref]>maxDist:
                znew[i] = float('NaN')
            
    return znew
    