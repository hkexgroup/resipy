# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:11:01 2018, python 3.6.6
vectorised isinpolygon algorithm using numpy, ask jamyd91 for a pyx version, it may give a slight performance boost. 
@author: jamyd91
"""
import numpy as np

def ccw(p0,p1,q0,q1,r0,r1):#code expects points as p0=x, p1=y, etc... 
    #Check the ordering of points 
    val=((q1-p1)*(r0-q0)) - ((q0-p0)*(r1-q1))
    order = np.zeros_like(val, dtype=int)
    #order[val==0] = 0 # lines are colinear
    order[val>0] = 1 # points are oreintated counter clockwise 
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
    ans: boolian, array like
        true indexes where point is inside polygon 
    Notes
    ----------
    polyx: array like
        polygon x coordinates 
    polyy: array like
        polygon y coordinates
    """
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
        count[idx]+=1 # add one to the indexes where true 
        #0 means the point is to the right of the polygon 
        #2 means the point is to the left of the polygon
        #1 indicates the point is inside the polygon 
        
        #NB: if a you have complex geometry the number of times the point is counted inside a polygon will be more, 
        # say a star. Rule of thumb is that the number of counts should be odd.
        
    return count%2 != 0 #returns boolian indexes

def planeintersect(x,y,z,poly_x,poly_y,poly_z,ray_cast=float('inf')):
    pass
    
