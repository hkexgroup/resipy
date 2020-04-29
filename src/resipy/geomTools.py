# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:11:01 2018, python 3.6.6
vectorised isinpolygon/isinvolume algorithm using numpy
@author: jamyd91
"""
import numpy as np
import warnings

def ccw(p0,p1,q0,q1,r0,r1):#code expects points as p0=x, p1=y, etc... 
    #Check the ordering of points 
    val=((q1-p1)*(r0-q0)) - ((q0-p0)*(r1-q1))
    order = np.zeros_like(val, dtype=int)
    #order[val==0] = 0 # lines are colinear
    order[val>0] = 1 # points are oreintated clockwise 
    order[val<0] = 2 # points are counter clockwise
    return order
    
def isinpolygon(x,y,poly_data,ray_cast=10000):
    """
     "Is point inside region?" code
    following solution posted at:
    https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
    
    Parameters
    ----------
    x: array like, float
        x coordinate of query point
    y: array like, float
        y coordinate of query point 
    poly_data: tuple 
        contains column of x and y coordinates, in the form (polyx, polyy)
    ray_cast: float, optional
        determines how the far to the right a point is ray casted 
    Returns
    ----------
    inside: boolian, numpy array 
        true indexes where point is inside volume
    Notes
    ----------
    polyx: array like
        polygon x coordinates 
    polyy: array like
        polygon y coordinates
    """
    warnings.warn('Depreciation Warning: Inispolygon is obsolete, please use matplotlib.path module instead')
    polyx = poly_data[0]
    polyy = poly_data[1]
    if len(polyx)!=len(polyy):
        raise ValueError('polygon vertices, xy coordinate, arrays are not the same length, check polygon generation scheme!')
    if len(polyx)<3:
        raise ValueError('A polygon cannot have less than 3 coordinates')
    if len(x) != len(y):
        raise ValueError('The number of x coordinates doesnt match the numer of y coordinates')

    Bx=np.array(x)+ray_cast
    By=np.array(y)
    #line AB is the coordinate and raycasted coordinate (B)
    Ax=np.array(x)
    Ay=np.array(y)
    #line CD will represent a segment from the polygon
    count=np.zeros_like(Ax, dtype = int)#number of times an intersection has been detected 
    tmp=np.zeros_like(Ax, dtype = int)
    for i in range(len(polyx)):
        if i==len(polyx)-1:#its the last iteration
            Cx=polyx[-1]
            Cy=polyy[-1]#take last and first polygon point
            Dx=polyx[0]
            Dy=polyy[0]
        else:
            Cx=polyx[i]
            Cy=polyy[i]
            Dx=polyx[i+1]
            Dy=polyy[i+1]
            
        idx = (ccw(Ax,Ay,Bx,By,Cx,Cy) != ccw(Ax,Ay,Bx,By,Dx,Dy)) & (ccw(Cx,Cy,Dx,Dy,Ax,Ay) != ccw(Cx,Cy,Dx,Dy,Bx,By))
        #tmp[idx] = ccw(Cx,Cy,Dx,Dy,Ax,Ay)
        count[idx]+=1 # add one to the indexes where true 
        #0 means the point is to the right of the polygon 
        #2 means the point is to the left of the polygon
        #1 indicates the point is inside the polygon 
        
        #NB: if a you have complex geometry the number of times the point is counted inside a polygon will be more, 
        # say a star. Rule is that the number of counts should be odd.
        
    return count%2 != 0#returns boolian indexes

#%% determine if line intersects plane 
def planeintersect(Ax, Ay, Az, 
                   Bx, By, Bz,
                   plane_x,
                   plane_y,
                   plane_z):#pass all variables as arguments so there is no need for unpacking 
    """
    Does lines AB intersect a plane in 3D coordinates 
    """
    count=np.zeros_like(Ax, dtype = int)#number of times an intersection has been detected 
    intersect = np.zeros_like(Ax, dtype=int)
    for i in range(len(plane_x)):
        if i==len(plane_x)-1:#its the last iteration
            Cx=plane_x[-1]
            Cy=plane_y[-1]#take last and first planegon point
            Cz=plane_z[-1]
            Dx=plane_x[0]
            Dy=plane_y[0]
            Dz=plane_z[0]
        else:
            Cx=plane_x[i]
            Cy=plane_y[i]
            Cz=plane_z[i]
            Dx=plane_x[i+1]
            Dy=plane_y[i+1]
            Dz=plane_z[i+1]
        #line CD will represent a segment from the planegon
        idxY = (ccw(Ax,Ay,Bx,By,Cx,Cy) != ccw(Ax,Ay,Bx,By,Dx,Dy)) & (ccw(Cx,Cy,Dx,Dy,Ax,Ay) != ccw(Cx,Cy,Dx,Dy,Bx,By))
        idxZ = (ccw(Ax,Az,Bx,Bz,Cx,Cz) != ccw(Ax,Az,Bx,Bz,Dx,Dz)) & (ccw(Cx,Cz,Dx,Dz,Ax,Az) != ccw(Cx,Cz,Dx,Dz,Bx,Bz))
        count[idxY]+=1 # add one to the indexes where true
        count[idxZ]+=1 # add one to the indexes where true
        
    hit = np.argwhere(count>1)
    intersect[hit]=1 
    
    return intersect # return binary 

#%% is inside a volume
def isinvolume(x,y,z,volume_data,ray_cast=10000):
    """
    Determine if a point lies inside a bounding volume (The code will not check 
    if supplied planes are bound).
    
    Parameters
    ----------
    x: array like, float
        x coordinate of query point
    y: array like, float
        y coordinate of query point 
    z: array like, float
        z coordinate of query point
    volume_data: list 
        contains column of poly_data, in the form (polyx, polyy, polyz). Each 
        poly_data argument describes the side of a volume. 
    ray_cast: float, optional
        determines how the far in the x axis a point is ray casted 
    Returns
    ----------
    inside: boolian, numpy array 
        true indexes where point is inside volume
    Notes
    ----------
    polyx: array like
        polygon x coordinates 
    polyy: array like
        polygon y coordinates
    polyz: array like
        polygon z coordinates
    """
    if not isinstance(volume_data,list):
        raise TypeError("Expected list type argument for 'volume_data'")
    Ax=np.array(x) # make coordinates into numpy arrays 
    Ay=np.array(y)
    Az=np.array(z)
    num_faces=len(volume_data)
    hit_array=np.zeros((num_faces,len(x)))
    #sort lines AB
    Bx = Ax + ray_cast
    for i in range(num_faces):
        polyx = volume_data[i][0]#retrive face data
        polyy = volume_data[i][1]
        polyz = volume_data[i][2]
        hits = planeintersect(Ax,Ay,Az,Bx,Ay,Az,polyx,polyy,polyz) # find 0 1 array of hits 
        hit_array[i,:]=hits#append into result array
    
    count = np.sum(hit_array,0)#sum results along the horizontal axis
    
    #the number of intersects with a plane should be  odd   
    return count%2 != 0 #returns boolian indexes


#%% determine if points are inside cubiod  
def in_box(x,y,z,xmax,xmin,ymax,ymin,zmax,zmin):
    """
    Determine if a point lies inside a bounding volume 
    Parameters
    ----------
    x: array like, float
        x coordinate of query point
    y: array like, float
        y coordinate of query point 
    z: array like, float
        z coordinate of query point
    volume_data: list 
        contains column of poly_data, in the form (polyx, polyy, polyz)
    ray_cast: float, optional
        determines how the far in the x axis a point is ray casted 
    Returns
    ----------
    inside: boolian, numpy array 
        true indexes where point is inside volume
    """
    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    #first check x coordinates 
    idx_x_in = (x>xmin) & (x<xmax)
    #check y coordinates 
    idx_y_in = (y>ymin) & (y<ymax)
    #check Z coordinates 
    idx_z_in = (z>zmin) & (z<zmax)
    #finally
    idx = (idx_x_in==True) & (idx_y_in==True) & (idx_z_in==True)
    return idx