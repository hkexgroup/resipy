#Topography interpolation schemes for ResIPy, python 3. Dated 2019. 
#author: jimmy boyd 
import sys
import numpy as np
import resipy.isinpolygon as iip
import warnings 

#%% compute thin plate spline /bilinear models  for irregular grid
# see solution @ https://math.stackexchange.com/questions/828392/spatial-interpolation-for-irregular-grid
def thin_plate_spline_mod(x,y,z): # apparently unstable 
    """Returns the thin plate spline model for 4 points  
    Needs lagrange multiplier solution. 
    """
    if len(x) != len(y) and len(x)!=len(z):
        raise ValueError("Mismatch in the number of elements in either x, y, or z arrays")
    #following solution posted at https://math.stackexchange.com/questions/828392/spatial-interpolation-for-irregular-grid
    dimn = len(x)
    X = np.matrix(np.concatenate((x**2,x*y,y**2,x,y,np.ones((dimn,1))),axis=1))
    E = np.zeros((6,6))
    for i in range(3):
        E[i,i]=1
    #construct inverse operator G 
    tG = np.concatenate((E,X.T),axis=1)
    bG = np.concatenate((X,np.zeros((4,4))),axis=1)
    G = np.concatenate((tG,bG))  
    #construct data column 
#    print(z.shape)
    D = np.concatenate((np.zeros((6,1)),z))
    #compute model
    mod = np.linalg.inv(G) * D 
    out = mod.A[0:6] # get first 6 parameters of modelled array 
    return out

#compute bilinear model - using system of linear equations 
def bilinear_mod(x,y,z):
    """Returns the bilinear model for 4 points. 
    """
    if len(x) != len(y) and len(x)!=len(z):
        raise ValueError("Mismatch in the number of elements in either x, y, or z arrays")
    #construct inverse operator G 
    dimn = len(x)
    G = np.matrix(np.concatenate((x*y,x,y,np.ones((dimn,1))),axis=1)) 
    #compute model
    mod = ((G.T * G)**-1) * G.T * z 
    return mod.A

def compute(mod,xnew,ynew):
    """ Compute the znew at xnew and ynew given a model from either 
    """
    shp = xnew.copy().shape
    if len(mod)==4:
        znew = mod[0]*xnew*ynew + mod[1]*xnew + mod[2]*ynew + mod[3]
    elif len(mod)==6:
        znew = mod[0]*xnew**2 + mod[1]*xnew*ynew + mod[2]*ynew**2 + mod[3]*xnew + mod[4]*ynew + mod[5]
    elif len(mod)==16: # then its a cubic model , not currently in use
        xnew.shape=(len(xnew),1)
        ynew.shape=(len(ynew),1)
        G = np.matrix(np.ones((len(xnew),16)))#construct operator 
        itr = 0
        for i in range(4):
            for j in range(4):
                G[:,itr] = xnew**i * ynew**j
                itr+=1
        znew = G*np.matrix(mod)
        znew = znew.A
        znew.shape = shp
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

#%% interpolation using 4 known points. 
def interp2d(xnew, ynew, xknown, yknown, zknown, extrapolate=True,method='bilinear'):
    """Compute z values for unstructured data using bilinear or spline interpolation. Coordinates
    outside the bounds of interpolation can be extrapolated using nearest neighbour
    algorithm on the interpolated and known z coordinates.  
    
    Bilinear /cubic /spline interpolation requires knowledge of 4 points orientated around the 
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
    method: 
        kind of interpolation done, options are bilinear and spline (more to come)
    extrapolate: bool, optional
        Flag for if extrapolation is to be used if new coordinates lie outside 
        the bounds where it is not possible to interpolate a value. 
        
    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew.
        
    """
    #check if inputs are arrays
    if len(xknown)==0 or len(xnew)==0:
        raise ValueError('Empty array passed to interp2d!')
    return_list = False 
    if type(xnew) != 'numpy.ndarray':
        xnew = np.array(xnew)
        #return_list = True
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
        
        
    #preallocate array for new z coordinates / interpolated values  
    znew = np.zeros_like(xnew)
    znew.fill(np.nan)
    #outside = np.logical_not(inside)
    num_pts = len(xnew)
    
    #add a bit of padding to prevent artefacts
    x_diff = np.min(np.diff(np.sort(xknown)))
    y_diff = np.min(np.diff(np.sort(yknown)))
    fudgex=0.001*x_diff
    fudgey=0.001*y_diff
    #compute new values inside survey
    for i in range(num_pts):
        if np.isnan(znew[i]): # point is not already known, else there is no need to interpolate
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
                #then the conditions needed to interpolate in a quad are met
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
                
                if len(x.shape)==1:#bug fix to deal with numpy being finicky, arrays need 2 dimensions 
                    x.shape += (1,)
                    z.shape += (1,)#append 1 dimension to the numpy array shape (otherwise np.concentrate wont work)
                    y.shape += (1,)
               
                if method == 'spline':
                    mod = thin_plate_spline_mod(x,y,z)
                else:
                    mod = bilinear_mod(x,y,z) # generate model   
                    
                inside_quad = iip.isinpolygon(xnew,ynew,(x,y))# any points which also fall inside the quad are interpolated
                znew[inside_quad] = compute(mod,xnew[inside_quad],ynew[inside_quad])# interpolate point using model 
        
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
        
    if return_list:
        znew = list(znew)
        
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
    """Nearest neighbour look up for 2D unstructured data.
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

#%% nearest neighbour interpolation in 3D 
def nearest3d(xnew,ynew,znew,xknown, yknown, zknown, iknown, return_idx=False):
    """Nearest neighbour look up for 3D unstructured data. This process can be 
    RAM and CPU demanding. 
    
    Parameters
    ------------
    xnew: array like
        x coordinates for the interpolated values
    ynew: array like 
        y coordinates for the interpolated values
    znew: array like 
        z coordinates for the interpolated values
    xknown: array like
        x coordinates for the known values 
    yknown: array like
        y coordinates for the known values 
    zknown: array like
        z coordinates for the known values 
    iknown: array like
        known values in 3D space. 
    return_idx: bool 
        Also return the look indexes of the iknown array. 

    Returns
    ------------
    inew: numpy array
        extrapolated values at (xnew ynew znew) coordinates 
        
    """
    #error checking 
    if len(xnew) != len(ynew) or len(xnew) != len(znew):
        raise ValueError('Mismatch in interpolated coordinate array lengths')
    if len(xknown) != len(yknown) or len(xknown) != len(zknown) or len(xknown) != len(iknown):
        raise ValueError('Mismatch in known coordinate array lengths')
    #as this process can take some time to compute progress is output to the screen
    sys.stdout.write('Running 3D nearest neighbour interpolation job...\n')

    sys.stdout.write('Computing delta x matrix ... ')
    a,b = np.meshgrid(xknown,xnew)
    dX = (a - b)**2
    sys.stdout.flush()
    sys.stdout.write('\rComputing delta y matrix ... ')
    a,b = np.meshgrid(yknown,ynew)
    dY = (a - b)**2
    sys.stdout.flush()
    sys.stdout.write('\rComputing delta z matrix ... ')
    a,b = np.meshgrid(zknown,znew)
    dZ = (a - b)**2
    #use meshgrid function to make a euclidian matrix 
    
    sys.stdout.flush()
    sys.stdout.write('\rComputing distance matrix .. ')
    dist = np.sqrt(dX+dY+dZ)#compute distance matrix 
    
    idx = np.argmin(dist,axis=1)
    sys.stdout.flush()
    sys.stdout.write('\rDone ...                     \n')   
    
    if return_idx:
        return iknown[idx], idx
    else: 
        return iknown[idx]

