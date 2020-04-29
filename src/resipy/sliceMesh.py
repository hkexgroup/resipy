#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 12:50:36 2018

@author: jkl
"""

import matplotlib.tri as tri
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider



#%% multidimensional slicing !

def sliceMesh(nodes, elms, values, label='', dim=2 ,vmin=None, vmax=None, ax=None):
    """ Produces an interactive graph of a 3D mesh with slice according to `axis`.
    
    Parameters
    ----------
    nodes : numpy.array
        Matrix with the nodes position (3 columns for x,y,z).
    elms : numpy.array
        Connection matrix to form elements.
    values : numpy.array
        Value for each element.
    label : str, optional
        String of the label for the colorbar.
    dim : int, optional
        Dimension to slice, 0:x, 1,y, 2:z
    vmin : float, optional
        Minimum value for colorbar.
    vmax : float, optional
        Maximum value for colorbar.
    ax : matplotlib.Axes, optional
        If provided, the plot will be drawn on this axis if not a new figure
        is created.
        
    Returns
    -------
    slider : matplotlib.widget.Slider
    
    Notes
    -----
    This function requires matplotlib interactive mode `plt.ion()`.
    """
    # compute centroids of those elements (or get mesh.elm_centre)
    centroids = np.array([np.mean(nodes[elm,:], axis=0) for elm in elms])

    def getSlice(elms, nodes, orig=None, normal=None):
        # let's define a plan in the Z direction
        if orig is None:
            orig = np.mean(nodes, axis=0)
        if normal is None:
            normal = np.array([0,0,1])
        
        # check on which side are those nodes
        iside = np.array([np.dot(nod-orig, normal) < 0 for nod in nodes])
        
        # element with nodes on both side are intersected by the plane
        elmCount = np.sum(iside[elms], axis=1)
        ielm = (elmCount > 0) & (elmCount < elms.shape[1])
    #    print(np.sum(ielm), '/', ielm.shape[0], 'intersected')
        
        return ielm # boolean array, True if element intersect the plane
    
    
    def getNormal(i):
        x = np.zeros(3)
        x[i] = 1
        return x
    
    def getOrig(i, val):
        x = np.zeros(3)
        x[i] = val
        return x
    
    def otherDim(i):
        dims = np.array([0,1,2])
        otherDims = np.delete(dims, i)
        return otherDims
    
    
    if dim == 0:
        title = 'X slice'
    elif dim == 1:
        title = 'Y slice'
    elif dim == 2:
        title = 'Z slice'
        
    ix, iy = otherDim(dim)
    normal = getNormal(dim)
    
    if vmin is None:
#        vmin = np.percentile(values,1)
        vmin = np.nanmin(values)
    if vmax is None:
#        vmax = np.percentile(values,99)
        vmax = np.nanmax(values)
    levels = np.linspace(vmin, vmax, 15)
    levels = None # the level with change for each slice
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    ax.set_title(title)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ielm = getSlice(elms, nodes, normal=normal)
    triang = tri.Triangulation(centroids[ielm,ix],centroids[ielm,iy])
    cax = ax.tricontourf(triang, values[ielm], levels=levels)
    #ax.triplot(triang, 'k-')
    #ax.scatter(centroids[:,1], centroids[:,2], s=5, c=values[ielm])          
    fig.colorbar(cax, ax=ax, label=label)
    
    def callback(i):
        ielm = getSlice(elms, nodes, orig=getOrig(dim,i), normal=normal)
        if np.sum(ielm) > 3:
            ax.clear()
            triang = tri.Triangulation(centroids[ielm,ix],centroids[ielm,iy])
            ax.tricontourf(triang, values[ielm], levels=levels)
            fig.canvas.draw()
    
    ax0 = plt.axes([0.25, 0.1, 0.65, 0.03])
    
    global meshSliceSlider
    meshSliceSlider = Slider(ax0, title, np.min(centroids[:,dim]),
                             np.max(centroids[:,dim]), 
                             valinit=np.mean(centroids[:,dim]))
    
    meshSliceSlider.on_changed(callback)
    
