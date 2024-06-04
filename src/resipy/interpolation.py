#Topography interpolation schemes for ResIPy, python 3. Dated 2019. 
#author: jimmy boyd 
import sys, warnings 
import numpy as np
from scipy.spatial import Delaunay, ConvexHull, cKDTree
from scipy.interpolate import LinearNDInterpolator
import matplotlib.path as mpltPath 

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

#%% quadrangle construction
def angles_in_quad(x,y):
    if len(x)!=4 or len(y)!=4:
        raise ValueError('Input coordinates do not make a quad!')
    theta = np.zeros(4)
    idx1 = [1,2,3,0]
    idx2 = [3,0,1,2]
    for i in range(4):
        dx10 = x[idx1[i]]-x[i]
        dx20 = x[idx2[i]]-x[i]
        dy10 = y[idx1[i]]-y[i]
        dy20 = y[idx2[i]]-y[i]
        m01 = np.sqrt( dx10*dx10 + dy10*dy10 )
        m02 = np.sqrt( dx20*dx20 + dy20*dy20 )
        theta[i] = np.arccos( (dx10*dx20 + dy10*dy20) / (m01 * m02))
    return np.rad2deg(theta)

def nerve_centre(x,y):
    """Compute veronoi nerve cell centre from three points 
    """
    x = np.array(x[0:3])
    y = np.array(y[0:3])
    x.shape = (3,1)
    y.shape = (3,1)
    M12 = np.concatenate((x**2 + y**2,y,np.ones((3,1))),axis=1)
    M13 = np.concatenate((x**2 + y**2,x,np.ones((3,1))),axis=1)
    M11 = np.concatenate((x,y,np.ones((3,1))),axis=1)
    x0 = 0.5*np.linalg.det(M12)/np.linalg.det(M11)
    y0 = -0.5*np.linalg.det(M13)/np.linalg.det(M11)
    return x0,y0

#%% triangle centriod 
def tri_cent(p,q,r):
    """Compute the centre coordinates for a 2d triangle given the x,y coordinates 
    of the vertices.
            
    Parameters
    ----------
    p : tuple,list,np array
        Coordinates of triangle vertices in the form (x,y).
    q : tuple,list,np array
        Coordinates of triangle vertices in the form (x,y).
    r : tuple,list,np array
        Coordinates of triangle vertices in the form (x,y).
            
    Returns
    ----------
    coordinates : tuple
        In the format (x,y).    
    """
    Xm=(p[0]+q[0])/2
    Ym=(p[1]+q[1])/2
    k=2/3
    Xc=r[0]+(k*(Xm-r[0]))
    Yc=r[1]+(k*(Ym-r[1]))
    return(Xc,Yc)

#%% ordering of points 
def ccw(p,q,r):#code expects points as p=(x,y) and so on ... 
    """Checks if points in a triangle are ordered counter clockwise. When using R2,
    mesh nodes should be given in a counter clockwise order otherwise you'll 
    get nonsensical results.
    
    Parameters
    ----------
    p - tuple, list,
        The x y coordinates of the point/vertex ie. (x,y) in desired order
    q - " 
    r - " 
    
    Returns
    ----------
    0 if colinear points, 1 if clockwise order, 2 if points are ordered counter clockwise
    """
    val=((q[1]-p[1])*(r[0]-q[0]))-((q[0]-p[0])*(r[1]-q[1]))
    if val == 0:
        return 0 # lines are colinear
    elif val > 0:
        return 1 # points are oreintated clockwise
    elif val < 0:
        return 2 # points are oreintated counter clockwise 
    
    
def ccwv(p, q, r):  # vectorize version
    """Checks if points in a triangle are ordered counter clockwise (vectorised 
    version). When using R2, mesh nodes should be given in a counter clockwise 
    order otherwise you'll get nonsensical results. 
    
    Parameters
    ----------
    p : 2D array for point 1 of shape N x 2 (first column X, second column Y).
    q : idem
    r : idem
    
    Returns
    -------
    array with :
        0 if colinear points
        1 if clockwise order,
        2 if points are ordered counter clockwise
    """
    val=((q[:,1]-p[:,1])*(r[:,0]-q[:,0]))-((q[:,0]-p[:,0])*(r[:,1]-q[:,1]))
    out = np.zeros(len(val)) # assume all colinear
    out[val > 0] = 1 # oriented clockwise
    out[val < 0] = 2 # oriented counter clockwise
    return out
    
    
def check_tetra(x,y,z):
    """Check if 4 points in a tetrahedra are ordered counter clocwise
    Parameters
    -----------
    x : array like 
        x coordinates of tetrahedral cell, 4 by 1 array. 
    y : array like 
        y coordinates of tetrahedral cell, 4 by 1 array. 
    z : array like 
        z coordinates of tetrahedral cell, 4 by 1 array.
    
    Returns
    ----------
    0 if coplanar points, 1 if clockwise order, 2 if points are counter clockwise
    """

    p = np.array([x,y,z]).T
    
    S = np.cross(p[1] - p[0], p[2] - p[0])
    N = np.dot(S, p[0] - p[3])
    
    if N>0:
        return 1 # points are clockwise
    elif N<0:
        return 2 # points are counter clockwise
    else:
        return 0 # something dogdey has happened as all points are on the same plane 
    

#check order of 4 points  
def order_quad(x,y,raycast=100000):
    """Order 4 points in a quad 
    Parameters
    -------------
    x: np array
        x coordinates 
    y: np array
        y coordinates 
    Return
    -----------
    order: np array of ints 
        indexes of x y points ordered in counter clockwise fashion
    """
    #find starting x coordinate 
    start_idx = np.argmin(x)
    counts = np.count_nonzero(x==x[start_idx])
    
    #protect against colinear points on left hand side of quad 
    if counts >1:#take the point with the smallest y value 
        idx = np.argwhere(x==x[start_idx]) # get all the points with that angle and take clostest
        ytmp = y[idx]
        start_idx = idx[np.argmin(ytmp)]
    
    #order angle 
    ry = y[start_idx]-raycast # ray casted coordinate below min y point
    rx = x[start_idx]
    dx10 = rx-x[start_idx]
    dx20 = x-x[start_idx]
    dy10 = ry-y[start_idx]
    dy20 = y-y[start_idx]
    m01 = np.sqrt( dx10*dx10 + dy10*dy10 )
    m02 = np.sqrt( dx20*dx20 + dy20*dy20 )
    
    m02[start_idx] = np.nan # to avoid / by 0 warning 
    
    theta = np.arccos((dx10*dx20 + dy10*dy20) / (m01 * m02)) # find array angle in radians 
    theta[np.isnan(theta)] = 0
    order = np.argsort(theta) #order 
    return order

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

#%% geometric median calculation 
def fdist(u,x,y,z):
    dx = (x-u[0])**2
    dy = (y-u[1])**2
    dz = (z-u[2])**2
    return np.sqrt(dx+dy+dz)

def geometricMedian(x,y,z, max_iter=30): 
    """
    Find the Geometric Median of a set of points using Weiszfeld's algorithm. 
    See here: https://en.wikipedia.org/wiki/Geometric_median#cite_note-vz-13
    

    Parameters
    ----------
    x : array like 
        DESCRIPTION.
    y : array like 
        DESCRIPTION.
    z : array like 
        DESCRIPTION.

    Returns
    -------
    gm: nd array 
        3 by 1 vector with the xyz coordinate of the point.

    """
    # make arrays in correct format 
    x = np.asarray(x,dtype=float)
    y = np.asarray(y,dtype=float)
    z = np.asarray(z,dtype=float)
    
    # values cannot be absolutely 0 
    x[x==0] = 1e-16
    y[y==0] = 1e-16
    z[z==0] = 1e-16 
    
    u = np.array([np.mean(x),np.mean(y),np.mean(z)]) # start point 
    r = fdist(u,x,y,z)
    if any(r==0):
        # shift start point 
        u = np.array([x[0], y[1], z[2]]) 
    
    pu = np.zeros(3,dtype=float) # new proposed point 
    delta = np.ones(3,dtype=float) # change between proposed and last model 
    c = 0 # counter for while loop 
    while sum(delta)>1e-16: 
        sqdist = ((x - u[0])**2)  + ((y - u[1])**2) + ((z - u[2])**2)
        dist = np.sqrt(sqdist)
        
        denom = sum(1/dist)
        pu[0] = sum(x/dist)/denom 
        pu[1] = sum(y/dist)/denom 
        pu[2] = sum(z/dist)/denom 
        
        delta = np.abs(pu - u) # compute change 
        u = np.array([pu[i] for i in range(3)])
        c += 1 
        if c > max_iter:
            warnings.warn('could not converge when computing geometric median')
            break
        
    return u 

#%% interpolation using 4 known points over an irregular grid 
#### New interpolation scheme coming in 2020 ####
def interp2d(xnew, ynew, xknown, yknown, zknown, extrapolate=True,method='spline'):
    """Compute z values for unstructured data using bilinear or spline interpolation. Coordinates
    outside the bounds of interpolation can be extrapolated using nearest neighbour
    algorithm on the interpolated and known z coordinates. Code firstly generates 
    a system of quadrangles (formed from deluany triangulation) and interpolates 
    over an irregular grid. 
    
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
    if type(xnew) != 'numpy.ndarray':xnew = np.array(xnew)
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
        
    # construct quadrangles / irregular quads 
    points = np.array([xknown,yknown]).T
    chull = ConvexHull(points) # convex hull 
    tri = Delaunay(points) # triangulation

    cindex = chull.vertices # get chull indices 
    
    tindex = tri.simplices.copy()
    nieghbours =  tri.neighbors
    
    vorx = np.zeros(len(tindex),dtype=float)
    vory = np.zeros(len(tindex),dtype=float)
    
    for i in range(len(tindex)):
        vorx[i],vory[i] = nerve_centre(xknown[tindex[i,:]],yknown[tindex[i,:]])
    
    path = mpltPath.Path(np.array([xknown[cindex],yknown[cindex]]).T)
    #inside = isinpolygon(vorx,vory,xknown[cindex],yknown[cindex])# ignore triangles with nerve centres outside the convex hull  
    inside = path.contains_points(np.array([vorx,vory]).T) # now using matplotlib function
    
    indexed_pairs = np.array([False]*len(vorx),dtype=bool)
    skip = 0
    vert = [] # quadrangle matrix >>> vertices of polygons, referenced later 
    
    #go through triangles and pair them up to form quads (an irregular grid)
    for i in range(len(tindex)):
        if inside[i]:
            nvorx = vorx[nieghbours[i]]
            nvory = vory[nieghbours[i]]
            dx = vorx[i] - nvorx
            dy = vory[i] - nvory
            dist = np.sqrt(dx**2 + dy**2)
            best_match = nieghbours[i,np.argmin(dist)]
            #check best match doesnt have a better canidate
            ncheckx = vorx[nieghbours[best_match]]
            nchecky = vory[nieghbours[best_match]]
            cdx = vorx[best_match] - ncheckx
            cdy = vory[best_match] - nchecky
            cndist = np.sqrt(cdx**2 + cdy**2) 
            rcheck = nieghbours[best_match,np.argmin(cndist)]          
            idx = np.append(tindex[i],tindex[best_match])
            if rcheck != i: # then use the next best neighbour ... 
                ri = np.argsort(dist)
                best_match = nieghbours[i,ri[1]]
                idx = np.append(tindex[i],tindex[best_match])
            
            if indexed_pairs[i] or indexed_pairs[best_match]:
                #then the triangles have already been paired skip ahead 
                skip += 1
    
            else:
                indexed_pairs[i]=True;indexed_pairs[best_match]=True
                unix = np.unique(idx)
            
                xu = xknown[unix];yu = yknown[unix] # unique x values , unique y values
                zu = zknown[unix]
                
                order = order_quad(xu,yu) # order vertices counter clockwise
        
                xuf = xu[order] # reorder the vertices counter clockwise
                yuf = yu[order] # reorder the vertices counter clockwise
                zuf = zu[order]
         
                vert.append(list(zip(xuf, yuf, zuf)))
                    
    #preallocate array for new z coordinates / interpolated values  
    znew = np.zeros_like(xnew)
    znew.fill(np.nan)
    
    #go through each quad in the quadrangle matrix (vert)
    #print('interpolating over irregular grid')
    for i in range(len(vert)):
        poly = np.array(vert[i])
        path = mpltPath.Path(np.array([poly[:,0],poly[:,1]]).T) 
        selection = path.contains_points(np.array([xnew,ynew]).T) # now using matplotlib function
        #selection = isinpolygon(xnew,ynew,poly[:,0],poly[:,1])# get interpolated values inside quad 
        
        x = poly[:,0] # quad vertices 
        y = poly[:,1]
        z = poly[:,2]
        
        #if len(x.shape)==1:#bug fix to deal with numpy being finicky, arrays need 2 dimensions 
        x.shape += (1,)
        z.shape += (1,)#append 1 dimension to the numpy array shape (otherwise np.concentrate wont work)
        y.shape += (1,)
        
        # generate model    
        if method == 'spline':
            mod = thin_plate_spline_mod(x,y,z) 
        else:
            mod = bilinear_mod(x,y,z) 
        #compute interpolated points     
        znew[selection] = compute(mod,xnew[selection],ynew[selection])

    znew = np.array(znew,dtype='float64').flatten()#need to specify data type on unix 
    
    idx_nan = np.isnan(znew).flatten() # boolian indexes of where nans are
    idx_num = np.invert(idx_nan).flatten()
    
    #extrapolate nans using nearest nieghbough interpolation(cKDTree)
    if extrapolate:
        #combine known and interpolated values 
        xknown = np.append(xknown,xnew[idx_num])
        yknown  = np.append(yknown,ynew[idx_num])
        zknown  = np.append(zknown,znew[idx_num])
        xextrap = xnew[idx_nan] 
        yextrap = ynew[idx_nan]
        
        # extrapolate using the gridded and interpolated data
        pknown = np.array([xknown,yknown]).T # known points 
        pextrap = np.array([xextrap,yextrap]).T # extrapolated points 
        tree = cKDTree(pknown)#tree object 
        dist,idx = tree.query(pextrap)# map known points to new points 
        
        znew[idx_nan] = zknown[idx]
        
    if return_list:
        znew = list(znew)
        
    return znew # return new interpolated values 

#%% interpolation using 4 known points - legacy. 
def interp2d_old(xnew, ynew, xknown, yknown, zknown, extrapolate=True,method='bilinear'):
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
#    znew = np.zeros_like(xnew)
#    znew.fill(np.nan)
    #outside = np.logical_not(inside)
    num_pts = len(xnew)
    
    #add a bit of padding to prevent artefacts
    x_diff = np.min(np.diff(np.sort(xknown)))
    y_diff = np.min(np.diff(np.sort(yknown)))
    fudgex=0.001*x_diff
    fudgey=0.001*y_diff
    
    def pnt_interp(xnew, ynew, xknown, yknown, zknown,fudgex,fudgey,method):
        quad1 = (xknown < xnew-fudgex) & (yknown < ynew-fudgey) # bottom left quad
        quad2 = (xknown < xnew-fudgex) & (yknown > ynew+fudgey) # top left quad
        quad3 = (xknown > xnew+fudgex) & (yknown > ynew+fudgey) # top right quad
        quad4 = (xknown > xnew+fudgex) & (yknown < ynew-fudgey) # bottom right quad
        
        dist1 = pdist(xnew, ynew, xknown[quad1], yknown[quad1])#distances to each quad 
        dist2 = pdist(xnew, ynew, xknown[quad2], yknown[quad2])
        dist3 = pdist(xnew, ynew, xknown[quad3], yknown[quad3])
        dist4 = pdist(xnew, ynew, xknown[quad4], yknown[quad4])
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
            znew = compute(mod,xnew,ynew)
        else:
            znew = np.nan

        return znew
    
    #compute new values inside survey
    znew = [pnt_interp(xnew[i], ynew[i], xknown, yknown, zknown,fudgex,fudgey,method) for i in range(num_pts)]

    znew = np.array(znew,dtype='float64').flatten()#need to specify data type on unix 
    # its worth parallising the function in the future over 2 cores for more speed 
    
    idx_nan = np.isnan(znew).flatten() # boolian indexes of where nans are
    idx_num = np.invert(idx_nan).flatten()
    
    #extrapolate nans using nearest nieghbough interpolation(cKDtree)
    if extrapolate:
        #combine known and interpolated values 
        xknown = np.append(xknown,xnew[idx_num])
        yknown  = np.append(yknown,ynew[idx_num])
        zknown  = np.append(zknown,znew[idx_num])
        xextrap = xnew[idx_nan] 
        yextrap = ynew[idx_nan]
#        zextrap = znew[idx_nan]
        
        # extrapolate using the gridded and interpolated data
        pknown = np.array([xknown,yknown]).T # known points 
        pextrap = np.array([xextrap,yextrap]).T # extrapolated points 
        tree = cKDTree(pknown)#tree object 
        dist,idx = tree.query(pextrap)# map known points to new points 
#        zextrap = zknown[idx]
        
        znew[idx_nan] = zknown[idx]
        
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
    
    if type(xnew) != 'numpy.ndarray':xnew = np.array(xnew)
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
    
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

    #extrapolate nans using nearest nieghbough interpolation(cKDtree)
    if extrapolate:
        #combine known and interpolated values 
        xknown = np.append(xknown,xnew[idx_num])
        yknown  = np.append(yknown,ynew[idx_num])
        zknown  = np.append(zknown,znew[idx_num])
        xextrap = xnew[idx_nan] 
        yextrap = ynew[idx_nan]
#        zextrap = znew[idx_nan]
        
        # extrapolate using the gridded and interpolated data
        pknown = np.array([xknown,yknown]).T # known points 
        pextrap = np.array([xextrap,yextrap]).T # extrapolated points 
        tree = cKDTree(pknown)#tree object 
        dist,idx = tree.query(pextrap)# map known points to new points 
#        zextrap = zknown[idx]
        
        znew[idx_nan] = zknown[idx]
        
    return znew

#%% interpolate using triangulation (can look odd in some situations)
def triangulate(xnew, ynew, xknown, yknown, zknown, extrapolate=True):
    """Tri-linear interpolation acheived using a triangulation scheme from Qhull,
    uses the scipy.interpolate.LinearNDinterpolator function. See here: 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html#scipy.interpolate.LinearNDInterpolator   
    
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

    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew.  
    """
    
    if type(xnew) != 'numpy.ndarray':xnew = np.array(xnew)
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
    
    #error checking 
    if len(xnew) != len(ynew):
        raise ValueError('Mismatch in interpolated coordinate array lengths')
    if len(xknown) != len(yknown) or len(xknown) != len(zknown):
        raise ValueError('Mismatch in known coordinate array lengths')
        
    pknown = np.array([xknown,yknown]).T # known points 
    pnew = np.array([xnew,ynew]).T # new points
    
    interpolator = LinearNDInterpolator(pknown,np.array(zknown))
    znew = interpolator(pnew)
    
    idx_nan = np.isnan(znew).flatten() # boolian indexes of where nans are
    idx_num = np.invert(idx_nan).flatten()
    
    #extrapolate nans using nearest nieghbough interpolation(cKDtree)
    if extrapolate:
        #combine known and interpolated values 
        xknown = np.append(xknown,xnew[idx_num])
        yknown  = np.append(yknown,ynew[idx_num])
        zknown  = np.append(zknown,znew[idx_num])
        xextrap = xnew[idx_nan] 
        yextrap = ynew[idx_nan]

        # extrapolate using the gridded and interpolated data
        pknown = np.array([xknown,yknown]).T # known points 
        pextrap = np.array([xextrap,yextrap]).T # extrapolated points 
        tree = cKDTree(pknown)#tree object 
        dist,idx = tree.query(pextrap)# map known points to new points 
        
        znew[idx_nan] = zknown[idx]
    
    return znew        

#%% interpolate using triangulation (can look odd in some situations)
def triangulate3d(xnew, ynew, znew, xknown, yknown, zknown, iknown, extrapolate=True):
    """Tri-linear interpolation acheived using a triangulation scheme from Qhull,
    uses the scipy.interpolate.LinearNDinterpolator function. See here: 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html#scipy.interpolate.LinearNDInterpolator   
    
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
    extropolate:
        use nearest nieghbour lookup to find nan values 

    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew.  
    """
    
    if type(xnew) != 'numpy.ndarray':xnew = np.array(xnew)
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(znew) != 'numpy.ndarray':znew = np.array(znew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
    if type(iknown) != 'numpy.ndarray':iknown = np.array(iknown)
    
    #error checking 
    if len(xnew) != len(ynew):
        raise ValueError('Mismatch in interpolated coordinate array lengths')
    if len(xknown) != len(yknown) or len(xknown) != len(zknown):
        raise ValueError('Mismatch in known coordinate array lengths')
        
    pknown = np.array([xknown,yknown,zknown]).T # known points 
    pnew = np.array([xnew,ynew,znew]).T # new points
    
    interpolator = LinearNDInterpolator(pknown,iknown)
    inew = interpolator(pnew)
    
    idx_nan = np.isnan(inew).flatten() # boolian indexes of where nans are
    idx_num = np.invert(idx_nan).flatten()
    
    #extrapolate nans using nearest nieghbough interpolation(cKDtree)
    if extrapolate:
        #combine known and interpolated values 
        xknown = np.append(xknown,xnew[idx_num])
        yknown  = np.append(yknown,ynew[idx_num])
        zknown  = np.append(zknown,znew[idx_num])
        iknown  = np.append(iknown,inew[idx_num])
        xextrap = xnew[idx_nan] 
        yextrap = ynew[idx_nan]
        zextrap = znew[idx_nan]

        # extrapolate using the gridded and interpolated data
        pknown = np.array([xknown,yknown,zknown]).T # known points 
        pextrap = np.array([xextrap,yextrap,zextrap]).T # extrapolated points 
        tree = cKDTree(pknown)#tree object 
        dist,idx = tree.query(pextrap)# map known points to new points 
        
        inew[idx_nan] = iknown[idx]
    
    return inew        
    

#%% pure nearest neighbour interpolation
def nearest(xnew, ynew, xknown, yknown, zknown, return_idx=False,
            num_threads=2): 
    """Nearest neighbour lookup using scipy.spatial.cKDTree
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
    return_idx: bool 
        Also return the look up indexes of the zknown array. 

    Returns
    ------------
    znew: numpy array
        z coordinates at xnew and ynew.  
        
    """
    
    if type(xnew) != 'numpy.ndarray':xnew = np.array(xnew)
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
    
    #error checking 
    if len(xnew) != len(ynew):
        raise ValueError('Mismatch in interpolated coordinate array lengths')
    if len(xknown) != len(yknown) or len(xknown) != len(zknown):
        raise ValueError('Mismatch in known coordinate array lengths')
    # >>> map the known points to the new points using cKDTree
    pknown = np.array([xknown,yknown]).T # known points 
    pnew = np.array([xnew,ynew]).T # new points 

    tree = cKDTree(pknown)#tree object 
    dist,idx = tree.query(pnew)# map known points to new points 
    
    if return_idx:
        return zknown[idx], idx
    else: 
        return zknown[idx]

#%% nearest neighbour interpolation in 3D 
def nearest3d(xnew,ynew,znew,xknown, yknown, zknown, iknown, return_idx=False,
              num_threads=2):
    """Nearest neighbour lookup using scipy.spatial.cKDTree
    
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
        Also return the look indexes up of the iknown array. 

    Returns
    ------------
    inew: numpy array
        extrapolated values at (xnew ynew znew) coordinates 
        
    """
    
    if type(xnew) != 'numpy.ndarray':xnew = np.array(xnew)
    if type(ynew) != 'numpy.ndarray':ynew = np.array(ynew)
    if type(xknown) != 'numpy.ndarray':xknown = np.array(xknown)
    if type(yknown) != 'numpy.ndarray':yknown = np.array(yknown)
    if type(zknown) != 'numpy.ndarray':zknown = np.array(zknown)
    
    #error checking 
    if len(xnew) != len(ynew) or len(xnew) != len(znew):
        raise ValueError('Mismatch in interpolated coordinate array lengths')
    if len(xknown) != len(yknown) or len(xknown) != len(zknown) or len(xknown) != len(iknown):
        raise ValueError('Mismatch in known coordinate array lengths')
    # >>> map the known points to the new points using cKDTree
    pknown = np.array([xknown,yknown,zknown]).T # known points 
    pnew = np.array([xnew,ynew,znew]).T # new points 
    tree = cKDTree(pknown)#tree object 
    dist,idx = tree.query(pnew)# map known points to new points 
    
    if return_idx:
        return iknown[idx], idx
    else: 
        return iknown[idx]

#%% National grid handling (needed for accurate interpolation)
def rotGridData(ngx, ngy, x0=None, y0=None, rotAngle=0):
    """ Rotate National grid coordinates into local coordinates. 
    
    Parameters
    ------------
    ngx: array like 
        Easting in National grid
    ngy: array like
        Northing in National grid
    x0: int, float, optional
        Initial easting
    y0: int, float, optional
        Initial northing
    rotAngle:
        Rotation angle (in degrees)
    
    Returns 
    ------------   
    rotx: array like
        Local x coordinate 
    roty: array like 
        Local Y coordinate 
    """
    if len(ngx) != len(ngy):
        raise ValueError('Got arrays of 2 different lengths in rotGridData')
        
    rotx = np.zeros(len(ngx))
    roty = np.zeros(len(ngx))
    
    if x0 is None:
        x0 = min(ngx)
    if y0 is None:
        y0 = min(ngy)
        
    rotAngle = np.deg2rad(rotAngle) # convert to radians for calculations 

    for i in range(len(ngx)):
        trans_x = ngx[i] - x0
        trans_y = ngy[i] - y0
        rotx[i] = trans_x*np.cos(rotAngle) - trans_y*np.sin(rotAngle)
        roty[i] = trans_x*np.sin(rotAngle) + trans_y*np.cos(rotAngle) 
        
    return rotx,roty

def invRotGridData(local_x, local_y, x0, y0, rotAngle=0):
    """ Rotate local coordinates back into National Grid
    Parameters
    ------------
    local_x: array like 
        Local X array 
    local_y: array like
        Local Y array 
    x0: int, float, 
        Initial easting. Same as 
    y0: int, float, 
        Initial northing
    rotAngle:
        Rotation angle (in radians)
    
    Returns 
    ------------   
    rotx: array like
        Rotated coordinate in National grid 
    roty: array like 
        Rotated coordinate in National grid 
    """
    if len(local_x) != len(local_y):
        raise ValueError('Got arrays of 2 different lengths in invRotGridData')
    rotx = np.zeros(len(local_x))
    roty = np.zeros(len(local_y))
    
    rotAngle = np.deg2rad(rotAngle) # convert to radians for calculations 
    for i in range(len(local_x)):
        col1 = np.array((np.cos(rotAngle),np.sin(rotAngle)))
        col2 = np.array((-np.sin(rotAngle),np.cos(rotAngle)))
        col1.shape=(2,1)
        col2.shape=(2,1)
        G = np.matrix(np.concatenate((col1,col2),axis=1)) 
        d = np.array([local_x[i], local_y[i]])
        d.shape = (2,1)
        mod = ((G.T * G)**-1) * G.T * d 
        t = mod.A
        rotx[i] = t[0][0]+ x0
        roty[i] = t[1][0]+ y0
        
    return rotx,roty
