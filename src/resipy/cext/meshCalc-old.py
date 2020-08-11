cimport cython  #import relevant modules 
import math as ma
import numpy as np
cimport numpy as np
from cython.parallel import prange, parallel 
cimport openmp

cdef extern from "math.h" nogil:
    cpdef double acos(double x)

DINT = np.intc
# mesh calculations - need to be in cython otherwise they will take too long 
@cython.boundscheck(False)#speeds up indexing, however it can be dangerous if the indexes are not in the range of the actual array, 
@cython.wraparound(False)        
def bisection_search(list arr, long long int var):
    """Efficent search algorithm for sorted lists of ints 
    Parameters
    -------------
    arr: list 
        sorted list 
    var: int
        item to be searched / indexed. If not found False is returned 
    """
    cdef int L = 0
    cdef int n = len(arr)
    cdef int R = n-1
    cdef int m
    while L <= R:
        m = ma.floor((L+R)/2)
        if arr[m]<var:
            L = m+1
        elif arr[m]>var:
            R = m-1
        else:
            return m
    return -1

@cython.boundscheck(False)
def unique(list arr): # unique with out numpy for an array of ints 
    """Find the unique values inside a list of ints (ONLY!!).Relies on the 
    efficeint bisection search as above. 
    
    Parameters 
    -------------
    arr: list 
        an unordered list of ints 
    
    Returns
    ------------
    uni: list
        list of unique ints in arr
    idx: list
        list of the corresponding indexes for the unique values in arr
    """
    cdef int i 
    cdef list temp_idx
    cdef list arrs
    cdef list uni, idx
    cdef int search 
    
    temp_idx = [iterat[0] for iterat in sorted(enumerate(arr), key=lambda x:x[1])] # sorted indexes  (pure python code)
    arrs = [arr[i] for i in temp_idx] #sorted array 
    
    #first element 
    uni = [arrs[0]]
    idx = [temp_idx[0]]
    
    for i in range(1,len(arrs)):
        search = bisection_search(uni,arrs[i])
        if search == -1:
            uni.append(arrs[i])
            idx.append(temp_idx[i])
    return uni,idx

def int2str(int a, int pad):
    cdef str s
    s = '{:0>%id}'%pad
    return s.format(a)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double fmin(double[:] arr) nogil: # import note: You need to pass the memory view
    # to a nogil function like this otherwise you get instability 
    cdef int n = arr.shape[0]
    cdef double tmp = arr[0]
    cdef int i 
    for i in range(1,n):
        if tmp>arr[i]:
            tmp = arr[i]
    return tmp

@cython.boundscheck(False)
@cython.wraparound(False)
cdef long tetra_signp(double[:] a, double[:] b, double [:] c, double[:] d,
                      double[:,:] vv, double[:] v1v, double[:] v2v,
                      double[:] v3v, double[:] Sv) nogil:
    #parrall nogil tetra sign code, but needs memory views passed to it 
    cdef double N
    cdef int k 

    for k in range(3):
        vv[k,0] = b[k] - a[k]
        vv[k,1] = c[k] - a[k]
        vv[k,2] = d[k] - a[k]
    
    for k in range(3): #extract delta vectors 
        v1v[k] = vv[k,0]
        v2v[k] = vv[k,1]
        v3v[k] = vv[k,2]
    
    #compute cross product
    Sv[0] = v1v[1]*v2v[2] - v1v[2]*v2v[1]
    Sv[1] = v1v[2]*v2v[0] - v1v[0]*v2v[2]
    Sv[2] = v1v[0]*v2v[1] - v1v[1]*v2v[0]
    #compute dot product (gives orientation)
    N = Sv[0]*v3v[0] + Sv[1]*v3v[1] + Sv[2]*v3v[2]
        
    if N>0:
        return 1 # points are clockwise
    elif N<0:
        return 2 # points are counter clockwise
    else:
        return 0 # something dogdey has happened as all points are on the same plane    
    
@cython.boundscheck(False)                 
def neigh3d(long[:,:] connection, int return_tri_combo):
    """Compute neighbours of each element within a 3D tetrahedral mesh 
    Parameters
    -------------
    connection: np.array 
        N by 4 array, describes how mesh nodes map to elements 
    return_tri_combo: int
        Binary, must be zero or 1. 1 to return the node combination matrix 
    
    Returns
    --------------
    neigh: np.array 
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    """
    
    cdef int i 
    cdef str idx1,idx2,idx3,idx4
    cdef int numel = connection.shape[0]
    cdef list face1s, face2s, face3s, face4s
    cdef str face1t, face2t, face3t, face4t
    cdef long long int face1, face2, face3, face4
    cdef tuple tri_combo = ([0]*numel,[0]*numel,[0]*numel,[0]*numel) # allocate space for tri_combo 
    cdef int num_node = max([max(connection[0]),max(connection[1]),max(connection[2]),max(connection[3])])
    cdef int pad = len(str(num_node))
    
    for i in range(numel):
        idx1 = int2str(connection[i,0],pad)#extract indexes 
        idx2 = int2str(connection[i,1],pad)
        idx3 = int2str(connection[i,2],pad)
        idx4 = int2str(connection[i,3],pad)
        
        # assign each face an organised and unique code
        #sort face indexes 
        face1s = sorted((idx2,idx3,idx4))
        face2s = sorted((idx1,idx4,idx3))
        face3s = sorted((idx1,idx2,idx4))
        face4s = sorted((idx1,idx2,idx3))
        face1t = str(face1s[0])+str(face1s[1])+str(face1s[2]) # add one to deal with 0 node number bug 
        face2t = str(face2s[0])+str(face2s[1])+str(face2s[2])
        face3t = str(face3s[0])+str(face3s[1])+str(face3s[2])
        face4t = str(face4s[0])+str(face4s[1])+str(face4s[2])
        face1 = int(face1t)
        face2 = int(face2t)
        face3 = int(face3t)
        face4 = int(face4t)
    
        tri_combo[0][i] = face1#face 1 
        tri_combo[1][i] = face2#face 2 
        tri_combo[2][i] = face3#face 3 
        tri_combo[3][i] = face4#face 4 
        
    cdef np.ndarray[long, ndim=2] neigh = np.zeros((numel,4),dtype=int) # allocate space for neighbour matrix        
    cdef long[:,:] neighv = neigh  
    
    #using binary search and sorted lists for efficient index lookup 
    cdef list tri_list = tri_combo[0] + tri_combo[1] + tri_combo[2] + tri_combo[3] #all faces together 
    cdef list temp_idx 
    cdef list tri_sort
	
    temp_idx = [iterat[0] for iterat in sorted(enumerate(tri_list), key=lambda x:x[1])] # sorted indexes  (pure python code)
    tri_sort= [tri_list[i] for i in temp_idx] # sorted faces (forward direction)
    cdef int maxo = len(tri_sort)-1
    cdef long long int o
    
    cdef list lookup 
    cdef list lookup_idx 
    lookup = [i for i in range(numel)]*4 # index lookup 
    lookup_idx = [lookup[i] for i in temp_idx] # sorted index lookup
    
    for i in range(numel):
        for j in range(4):
            o = bisection_search(tri_sort,tri_combo[j][i]) # only works on a sorted list 
            #find the reference index
            if lookup_idx[o]==i:# then there are 2 options 
                #1 ) the face in question is unique or
                #2 ) the index is o+1 or o-1
                if o!=maxo and tri_sort[o+1] == tri_combo[j][i]:
                    o+=1
                elif o!=0 and tri_sort[o-1] == tri_combo[j][i]:
                    o-=1
                else: # unique face on edge of mesh 
                    o = -1   
            
            if o==-1: 
                idx = -1  
            else:
                idx = lookup_idx[o]
            neighv[i,j] = idx
        
    if return_tri_combo==1:
        return neigh,tri_combo
    else:
        return neigh
    
@cython.boundscheck(False)
def faces3d(long[:,:] connection, long[:,:] neigh):
    """Return external faces of a 3D tetrahedral mesh 
    Parameters
    ----------
    connection: np.array (int)
        N by 4 array, describes how mesh nodes map to elements 

    neigh: np.array (int)
        N by 4 array, describes indices of neighbour elements for each cell 
        (output of neigh3d)

    Returns
    -------
    fconnection : np.array
        N by 3 matrix with the with the nodes of external faces 
    idxa: np.array
        N by 1 array. Indices of original mesh elements 

    """
    cdef int i,j #loop variables
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    cdef int c = 0
    cdef np.ndarray[long, ndim=2] nmap = np.array([[1, 2, 3], [0, 3, 2], 
                                                   [0, 1, 3], [0, 1, 2]])
    
    #first thing is to find all the cases where neigh == 1 
    cdef list idx = [] # this is the index where a outside edge is 
    cdef list fnum = [] # this is the face number ( 1 2 3 or 4)
    for i in range(numel):
        for j in range(npere):
            if neigh[i,j] == -1:
                idx.append(i)
                fnum.append(j)
                c+=1
    
    cdef np.ndarray[long, ndim=1] idxa = np.array(idx)
    cdef np.ndarray[long, ndim=1] fnuma = np.array(fnum)
    cdef long[:] fnumav = fnuma
    cdef long[:] idxav = idxa
    
    cdef int nfaces = len(idx)
    cdef int fidx, fnumi 
    
    cdef np.ndarray[long, ndim=2] fconnection = np.zeros((nfaces,3),dtype=int) # face connection matrix 
    cdef long[:,:] fconnectionv = fconnection
    
    for i in range(nfaces):
        fidx = idxav[i]
        fnumi = fnumav[i]
        for j in range(3):
            fconnectionv[i,j] = connection[fidx,nmap[fnumi,j]]
            #find missing node in future? 
        
    return fconnection, idxa

@cython.boundscheck(False)
def neigh2d(long[:,:] connection, int return_tri_combo=0):
    """Compute neighbours of each element within a 2D triangular mesh 
    Parameters
    -------------
    connection: np.array 
        N by 3 array, describes how mesh nodes map to elements 
    return_tri_combo: int
        Binary, must be zero or 1. 1 to return the node combination matrix 
    
    Returns
    --------------
    neigh: np.array 
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    """
    
    cdef int i
    cdef str idx1,idx2,idx3,idx4
    cdef int numel = connection.shape[0]
    cdef list face1s, face2s, face3s
    cdef str face1t, face2t, face3t
    cdef long long int face1, face2, face3
    cdef tuple tri_combo = ([0]*numel,[0]*numel,[0]*numel) # allocate space for tri_combo 
    cdef int num_node = max([max(connection[0]),max(connection[1]),max(connection[2])])
    cdef int pad = len(str(num_node))
    
    for i in range(numel):
        idx1 = int2str(connection[i,0],pad) # extract indexes 
        idx2 = int2str(connection[i,1],pad) 
        idx3 = int2str(connection[i,2],pad)
        
        # assign each face an organised and unique code
        #sort face indexes 
        face1s = sorted([idx1,idx2])
        face2s = sorted([idx2,idx3])
        face3s = sorted([idx3,idx1])

        face1t = str(face1s[0])+str(face1s[1]) 
        face2t = str(face2s[0])+str(face2s[1])
        face3t = str(face3s[0])+str(face3s[1])

        face1 = int(face1t)
        face2 = int(face2t)
        face3 = int(face3t)
    
        tri_combo[0][i] = face1#face 1 
        tri_combo[1][i] = face2#face 2 
        tri_combo[2][i] = face3#face 3 
        
    cdef np.ndarray[long, ndim=2] neigh = np.zeros((numel,3),dtype=int) # allocate space for neighbour matrix        
    cdef long[:,:] neighv = neigh             
    
    #using binary search and sorted lists for efficient index lookup 
    cdef list tri_list = tri_combo[0] + tri_combo[1] + tri_combo[2] #all faces together 
    cdef list temp_idx 
    cdef list tri_sort
	
    temp_idx = [iterat[0] for iterat in sorted(enumerate(tri_list), key=lambda x:x[1])] # sorted indexes  (pure python code)
    tri_sort= [tri_list[i] for i in temp_idx] # sorted faces (forward direction)
    cdef int maxo = len(tri_sort)-1
    cdef long long int o
    
    cdef list lookup 
    cdef list lookup_idx 
    lookup = [i for i in range(numel)]*4 # index lookup 
    lookup_idx = [lookup[i] for i in temp_idx] # sorted index lookup
    
    for i in range(numel):
        for j in range(3):
            o = bisection_search(tri_sort,tri_combo[j][i]) # only works on a sorted list 
            #find the reference index
            if lookup_idx[o]==i:# then there are 2 options 
                #1 ) the face in question is unique or
                #2 ) the index is o+1 or o-1
                if o!=maxo and tri_sort[o+1] == tri_combo[j][i]:
                    o+=1
                elif o!=0 and tri_sort[o-1] == tri_combo[j][i]:
                    o-=1
                else: # unique face on edge of mesh 
                    o = -1   
            
            if o==-1: 
                idx = -1  
            else:
                idx = lookup_idx[o]
            neigh[i,j] = idx
        
    if return_tri_combo==1:
        return neigh,tri_combo
    else:
        return neigh
    
@cython.boundscheck(False)
def split_tri(list connection, list node_x, list node_y, list node_z):
    """Split triangle elements down into 4 smaller triangles 
    Parameters
    -----------
    connection: list
        3 by N lists of the node numbers forming the triangle mesh
    node_x: list 
        Node x coordinates 
    node_y: list 
        Node x coordinates 
    node_z: list 
        Node z coordinates 
        
    Returns
    -----------
    new_connection: list
    outx: list
    outy: list
    outz: list
    num_elms: int
    num_nodes: int
    """
    #define variables 
    cdef i=0
    cdef list x = [0]*3
    cdef list y = [0]*3
    cdef list z = [0]*3
    cdef list n = [0]*6 # new node numbers 
    cdef float mx, my, mz
    cdef list nodes
    cdef unsigned long long int mn
    cdef int j, nno, tmpi #loop variables 

    cdef int num_nodes = len(node_x) # number of nodes 
    cdef int num_elms = len(connection[0])
    cdef list nnodex 
    cdef list nnodey # new node lists 
    cdef list nnodez 
    cdef list new_connection = [[0]*num_elms*4,[0]*num_elms*4,[0]*num_elms*4] # new connection matrix 
    cdef list new_map = [[0]*num_elms*3,[0]*num_elms*3,[0]*num_elms*3,[0]*num_elms*3,[0]*num_elms*3] # used to reference already computed node centres
    #the columns are as follows: nodex,nodey,nodez,node_config, original element no
    cdef int pad = len(str(num_nodes))+1
    
    cdef list a = [0,1,2]
    cdef list b = [1,2,0]
    
    #remap the element nodes with this matrix 
    cdef list remap = [[0, 5, 3],
                       [1, 4, 3],
                       [2, 5, 4],
                       [3, 4, 5]]
    
    cdef list unicl, idx, node_id 

    
    ### loop through the elements ### all possible node configurations 
    for i in range(0,num_elms):
        for j in range(3):
            nno = connection[j][i] # node number
            x[j] = node_x[nno]
            y[j] = node_y[nno]
            z[j] = node_z[nno]
            n[j] = nno 
            
        tmpi = 3*i
        for j in range(3):
            mx = (x[a[j]]+x[b[j]])/2
            my = (y[a[j]]+y[b[j]])/2
            mz = (z[a[j]]+z[b[j]])/2
            nodes = sorted([int2str(n[a[j]],pad), int2str(n[b[j]],pad)])
            mn = int(nodes[0]+nodes[1])
            
            new_map[0][tmpi+j] = mx
            new_map[1][tmpi+j] = my
            new_map[2][tmpi+j] = mz
            new_map[3][tmpi+j] = mn
            new_map[4][tmpi+j] = i
    
    ### find unique node configs # ###         
    unicl,idx = unique(new_map[3]) # unique and ordered node configs 
    node_id = [num_nodes + i for i in range(len(unicl))]
    
    ### map back to elements #### 
    for i in range(0,num_elms):
         
        for j in range(3):
            nno = connection[j][i] # node number
            n[j] = nno 
        
        for j in range(3):
            nodes = sorted([int2str(n[a[j]],pad), int2str(n[b[j]],pad)])
            mn = int(nodes[0]+nodes[1])
            search = bisection_search(unicl,mn)
            
            n[j+3] = node_id[search]#reference index
    
        tmpi = i*4 # temporary index for indexing new connection matrix
        for j in range(4):
            new_connection[0][tmpi+j] = n[remap[j][0]]
            new_connection[1][tmpi+j] = n[remap[j][1]]
            new_connection[2][tmpi+j] = n[remap[j][2]]
            
    nnodex = [new_map[0][i] for i in idx]
    nnodey = [new_map[1][i] for i in idx]
    nnodez = [new_map[2][i] for i in idx]
                
    cdef list outx = node_x+nnodex
    cdef list outy = node_y+nnodey
    cdef list outz = node_z+nnodez
                
    return new_connection, outx, outy, outz, len(new_connection[0]), len(node_x)+len(nnodex)

@cython.boundscheck(False)
def split_tetra(list connection, list node_x, list node_y, list node_z):
    """Split tetrahedral elements down into 8 smaller tetrahedra 
    Parameters
    -----------
    connection: list
        3 by N lists of the node numbers forming the triangle mesh
    node_x: list 
        Node x coordinates 
    node_y: list 
        Node x coordinates 
    node_z: list 
        Node z coordinates 
        
    Returns
    -----------
    new_connection: list
    outx: list
    outy: list
    outz: list
    num_elms: int
    num_nodes: int
    """
    cdef i=0
    cdef list x = [0]*4
    cdef list y = [0]*4
    cdef list z = [0]*4
    cdef list n = [0]*10 # new node numbers 
    cdef float mx, my, mz
    cdef list nodes
    cdef long long int mn
    cdef int j, nno, tmpi #loop variables 

    cdef int num_nodes = len(node_x)
    cdef int num_elms = len(connection[0])

    cdef list nnodex = [] 
    cdef list nnodey = [] # new node lists 
    cdef list nnodez = [] 

    cdef list new_connection = [[0]*num_elms*8,[0]*num_elms*8,[0]*num_elms*8,[0]*num_elms*8] # new connection matrix 
    cdef list new_map = [[0]*num_elms*6,[0]*num_elms*6,[0]*num_elms*6,[0]*num_elms*6,[0]*num_elms*6] # used to reference already computed node centres
    #the columns are as follows: nodex,nodey,nodez,node_config, original element no
    cdef int pad = len(str(num_nodes))
    
    cdef list a = [0,1,2,0,1,2]
    cdef list b = [1,2,0,3,3,3]
    
    #remap the element nodes with this matrix 
    cdef list remap = [[4, 7, 6, 0],
                       [5, 9, 6, 2],
                       [8, 9, 7, 3],
                       [8, 9, 7, 6],
                       [8, 4, 7, 6],
                       [8, 5, 9, 6],
                       [4, 5, 6, 8],
                       [4, 1, 5, 8]]
    
    cdef list unicl, idx, node_id 
    
    ### loop through the elements ### all possible node configurations 
    for i in range(0,num_elms):
        x = [0]*4
        y = [0]*4
        z = [0]*4
        n = [0]*10 # new node numbers 
        for j in range(4):
            nno = connection[j][i] # node number
            x[j] = node_x[nno]
            y[j] = node_y[nno]
            z[j] = node_z[nno]
            n[j] = nno 
            
        tmpi = 6*i
        for j in range(6):
            mx = (x[a[j]]+x[b[j]])/2
            my = (y[a[j]]+y[b[j]])/2
            mz = (z[a[j]]+z[b[j]])/2
            nodes = sorted([int2str(n[a[j]],pad), int2str(n[b[j]],pad)])
            mn = int(nodes[0]+nodes[1])
            
            new_map[0][tmpi+j] = mx
            new_map[1][tmpi+j] = my
            new_map[2][tmpi+j] = mz
            new_map[3][tmpi+j] = mn
            new_map[4][tmpi+j] = i
    
    ### find unique node configs # ###         
    unicl,idx = unique(new_map[3]) # unique and ordered node configs 
    node_id = [num_nodes + i for i in range(len(unicl))]
    
    ### map back to elements #### 
    for i in range(0,num_elms):
        n = [0]*10 # new node numbers 
        
        for j in range(4):
            nno = connection[j][i] # node number
            n[j] = nno 
        
        for j in range(6):
            nodes = sorted([int2str(n[a[j]],pad), int2str(n[b[j]],pad)])
            mn = int(nodes[0]+nodes[1])
            search = bisection_search(unicl,mn)
            
            n[j+4] = node_id[search]#reference index
    
        tmpi = i*8 # temporary index for indexing new connection matrix
        for j in range(8):
            new_connection[0][tmpi+j] = n[remap[j][0]]
            new_connection[1][tmpi+j] = n[remap[j][1]]
            new_connection[2][tmpi+j] = n[remap[j][2]]
            new_connection[3][tmpi+j] = n[remap[j][3]]
            
    nnodex = [new_map[0][i] for i in idx]
    nnodey = [new_map[1][i] for i in idx]
    nnodez = [new_map[2][i] for i in idx]
                
    cdef list outx = node_x+nnodex
    cdef list outy = node_y+nnodey
    cdef list outz = node_z+nnodez
                    
    return new_connection, outx, outy, outz, len(new_connection[0]), len(node_x)+len(nnodex)

@cython.boundscheck(False)
@cython.wraparound(False)   
def orderTetra(long[:,:] connection, double[:,:] node, int num_threads=2):
    """ Organise tetrahedral element nodes into a clockwise order
    
    following solution at: 
    https://stackoverflow.com/questions/10612829/tetrahedron-orientation-for-triangle-meshes
    
    Parameters
    ----------
    connection: np array (long)
        Mesh connection matrix, N by 4 array. 
    node: np array (float64)
        N by 3 array describing the mesh nodes 
    num_threads: int 
        Number of threads used during the operation. Default is 2. 
    Returns
    -------
    con : np array (long)
        Mesh connection matrix which has been 'corrected' such that all elements
        are counterclockwise. 
    count : int
        Number of elements with switched nodes (ie been 'corrected')
    ccw : np array(long)
        clockwise orientation vector
    """
    #get the number of nodes and elements etc 
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    #cdef int count = 0 #rolling total for the number of corrected elements 
    
    con = np.zeros((numel,npere), dtype=DINT)# new connection matrix 
    cdef int[:,:] conv = con # connection memeory view
    
    cdef np.ndarray[double, ndim=1] node_x = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef np.ndarray[double, ndim=1] node_y = np.asarray(node[:,1],dtype=float)
    cdef np.ndarray[double, ndim=1] node_z = np.asarray(node[:,2],dtype=float)
    cdef double[:] node_xv = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef double[:] node_yv = np.asarray(node[:,1],dtype=float)
    cdef double[:] node_zv = np.asarray(node[:,2],dtype=float)
    
    # cdef int num_threads = 2 # number of threads to use (using 2 for now)

    #looping variables 
    cdef np.ndarray[double,ndim=3] v = np.zeros((3,3,num_threads),dtype=float) # matrix
    cdef np.ndarray[double,ndim=2] v1 = np.zeros((3,num_threads)) # comparison vector 
    cdef np.ndarray[double,ndim=2] v2 = np.zeros((3,num_threads))
    cdef np.ndarray[double,ndim=2] v3 = np.zeros((3,num_threads))
    cdef np.ndarray[double,ndim=2] S = np.zeros((3,num_threads))
    #memory views of the above
    cdef double[:,:,:] vv = v
    cdef double[:,:] v1v = v1 
    cdef double[:,:] v2v = v2 
    cdef double[:,:] v3v = v3
    cdef double[:,:] Sv = S 
    
    cdef double[:] N = np.zeros(num_threads,dtype=float)# orientation indicator 
    cdef np.ndarray[long, ndim=1] ccw = np.zeros(numel,dtype=int) # clockwise array 
    cdef long[:] ccwv = ccw #clockwise view
    cdef Py_ssize_t i 
    cdef int k, ei, tid # loop integers 
    cdef np.ndarray[long, ndim=1] count = np.zeros(numel,dtype=int)
    cdef long[:] countv = count
    cdef int count_out

    with nogil, parallel(num_threads=num_threads):
        tid = openmp.omp_get_thread_num()
        for i in prange(numel):
            for k in range(3): # work out delta matrix 
                vv[0,k,tid] = node_xv[connection[i,k+1]] - node_xv[connection[i,0]] 
                vv[1,k,tid] = node_yv[connection[i,k+1]] - node_yv[connection[i,0]] 
                vv[2,k,tid] = node_zv[connection[i,0]] - node_zv[connection[i,k+1]] 
                
            for k in range(3): #extract delta vectors 
                v1v[k,tid] = vv[k,0,tid]
                v2v[k,tid] = vv[k,1,tid]
                v3v[k,tid] = vv[k,2,tid]
            
            #compute cross product
            Sv[0,tid] = v1v[1,tid]*v2v[2,tid] - v1v[2,tid]*v2v[1,tid]
            Sv[1,tid] = v1v[2,tid]*v2v[0,tid] - v1v[0,tid]*v2v[2,tid]
            Sv[2,tid] = v1v[0,tid]*v2v[1,tid] - v1v[1,tid]*v2v[0,tid]
            #compute dot product (gives orientation)
            N[tid] = Sv[0,tid]*v3v[0,tid] + Sv[1,tid]*v3v[1,tid] + Sv[2,tid]*v3v[2,tid]
            
            if N[tid]>0: # then tetrahedra is clockwise 
                countv[i] = 1
                conv[i,1] = connection[i,0]
                conv[i,0] = connection[i,1]
                ccwv[i] = 1
            elif N[tid]<0: # then its counter clockwise 
                conv[i,0] = connection[i,0]
                conv[i,1] = connection[i,1]
                ccwv[i] = 2
            else: # shouldn't be possible 
                ccwv[i] = 0
                ei = i #problem index 
                
            conv[i,2] = connection[i,2]
            conv[i,3] = connection[i,3]
        
    for i in range(numel):
        if ccw[i]==0 and ei>-1:
            raise ValueError('Element %i has all colinear nodes, thus is poorly formed and can not be ordered'%ei)
    count_out = sum(count)
            
    return con, count_out, ccw

@cython.boundscheck(False)
@cython.wraparound(False)
def orderQuad(long[:,:] connection, double[:,:] node):
    """Order quaderlateral elements counterclockwise 
    
    Parameters
    ----------
    connection: np array (long)
        Mesh connection matrix, N by 4 array. 
    node: np array (float64)
        N by 3 array describing the mesh nodes 

    Returns
    -------
    con : np array (long)
        Mesh connection matrix which has been 'corrected' such that all elements
        are counterclockwise. 
    count : int
        Number of elements with switched nodes (ie been 'corrected')
    """
    #get the number of nodes and elements etc 
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    
    if npere!=4: 
        raise ValueError('Got the wrong number of nodes per element when ordering mesh nodes')
    
    con = np.zeros((numel,npere), dtype=DINT)# new connection matrix 
    cdef int[:,:] conv = con # connection view
    
    cdef np.ndarray[double, ndim=1] node_x = np.asarray(node[:,0],dtype=float) 
    cdef np.ndarray[double, ndim=1] node_y = np.asarray(node[:,1],dtype=float)
    cdef np.ndarray[double, ndim=1] node_z = np.asarray(node[:,2],dtype=float)
    
    #loop varaibles 
    cdef np.ndarray[double, ndim=1] x = np.zeros(4,dtype=float) # x array 
    cdef np.ndarray[double, ndim=1] z = np.zeros(4,dtype=float) # z array 
    cdef np.ndarray[double, ndim=1] xtmp = np.zeros(4,dtype=float)
    cdef np.ndarray[double, ndim=1] ztmp = np.zeros(4,dtype=float)
    cdef double[:] xv = x # xview
    cdef double[:] zv = z # z view 
    cdef double[:] xtmpv = xtmp 
    cdef double[:] ztmpv = ztmp 

    cdef np.ndarray[long, ndim=1] count = np.zeros(numel,dtype=int)
    cdef long[:] countv = count
    cdef int count_out 
    
    cdef np.ndarray[double, ndim=1] theta = np.zeros(4,dtype=float) # theta array
    cdef double[:] thetav = theta # theta view 
    
    cdef np.ndarray[long, ndim=1] order = np.zeros(4,dtype=int) #ordering array
    cdef long[:] orderv = order #ordering view
    
    cdef int i, k, c, j
    cdef int eflag = 0 # error flag
    cdef int ei # error number 
    cdef int num_threads = 1
    cdef double minx, minz
    cdef double xinf = 10*max(node_x) # infinite like x 
    cdef double zinf = 10*max(node_z) # infinite like z 
    cdef double tinf = 314  #infinite like theta 
    cdef double mtheta
    cdef double rz, rx # ray casted coordinate
    cdef double dx10, dx20, dz10, dz20, m01, m02

    #aim is to work out angles from point with minimum x and z coordinate 
    #min angle >>> bottom left most point 
    #max angle >>> upper left most point
    
    for i in range(numel):
        #could make loop parallel in future? 
        for k in range(4):
            xv[k] = node_x[connection[i,k]]
            zv[k] = node_z[connection[i,k]]
        
        #find starting x coordinate 
        minx = fmin(xv)
        c = 0 # rolling count 
        #protect against colinear points on left hand side of quad 
        for k in range(4):
            if x[k] == minx:
                ztmpv[k] = zv[k]
                xtmpv[k] = xv[k] # put in array where x == min(x)
                c = c + 1
            else:
                ztmpv[k] = zinf
                xtmpv[k] = xinf
                
        minz = fmin(ztmpv)        
        if c>2: # then there is more than one min x coordinate (presumably 2 at most)
            eflag = 1 
            ei = i 
        
        #create order angle 
        for k in range(4):
            if xv[k] == minx and zv[k] == minz:
                thetav[k] = 0 # angle is zero becuase we are comparing on left bottom most point 
            else:
                rz = minz - zinf # ray casted coordinate below min z point
                rx = minx
                dx10 = rx - minx
                dx20 = x[k] - minx
                dz10 = rz - minz
                dz20 = z[k] - minz
                m01 = (dx10*dx10 + dz10*dz10)**0.5
                m02 = (dx20*dx20 + dz20*dz20)**0.5
                thetav[k] = acos((dx10*dx20 + dz10*dz20) / (m01 * m02))
       
        mtheta = fmin(thetav) # min theta 
        for j in range(4):
            for k in range(4):
                if thetav[k] == mtheta:
                    thetav[k] = tinf
                    orderv[j] = k
                    mtheta = fmin(thetav)
                    break  
        
        for k in range(4):# flag if order changes >>> count as ordered element 
            if orderv[k] != k:#order has been changed 
                countv[i] = 1 
                break
        
        for k in range(4):
            conv[i,k] = connection[i,order[k]]
                                
    if eflag == 1: 
        raise ValueError('Element %i has more than 2 colinear points and therefore not a quad'%ei)
        
    count_out = sum(count)
    
    return con, count_out 

@cython.boundscheck(False)
@cython.wraparound(False)
def surfaceCall(long[:,:] fconnection, double[:,:] node, double[:,:] cellcentres,
                int num_threads=2):
    """Call to determine if outward facing mesh faces are on the top surface 
    of a tetrahedral mesh. Roughly following solution posted at: 
    
    https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d

    Parameters
    ----------
    fconnection: np.array (int)
        Connection matrix which describes the faces of elements on the edge 
        of the mesh. N by 3. 
    node: np.array (float)
        Coordinates of mesh nodes (M by 3)
    cellcentres: np.array (float)           
        Coordinates of the centre point of the respective elements which are 
        on the outside of the mesh. 
    int num_threads, optional
        Number of threads to use during the operation. The default is 1.

    Returns
    -------
    ocheck : np.array(int)
        An array of 1 and 0. 1 means the face is likley on the top surface. 

    """
    #print('number of threads = %i'%num_threads)
    #pull out important arrays / info 
    cdef int nfaces = fconnection.shape[0]
    cdef double[:] node_xv = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef double[:] node_yv = np.asarray(node[:,1],dtype=float)
    cdef double[:] node_zv = np.asarray(node[:,2],dtype=float)
    cdef np.ndarray[double,ndim=1] nodez =  np.asarray(node[:,2],dtype=float)
    cdef double tqz = (max(nodez) - min(nodez)) + max(nodez) # query point in z axis
    
    cdef np.ndarray[long,ndim=1] ocheck = np.zeros(nfaces,dtype=int)
    cdef long[:] ocheckv = ocheck

    #looping variables 
    cdef int i,j, tid
    cdef double[:] bqz = np.zeros(num_threads)
    cdef double[:,:] x = np.zeros((3,num_threads))
    cdef double[:,:] y = np.zeros((3,num_threads))
    cdef double[:,:] z = np.zeros((3,num_threads))
    cdef double[:] xm = np.zeros(num_threads)
    cdef double[:] ym = np.zeros(num_threads)
    
    cdef double[:,:] q1 = np.zeros((3,num_threads))
    cdef double[:,:] q2 = np.zeros((3,num_threads))
    cdef double[:,:] p1 = np.zeros((3,num_threads))
    cdef double[:,:] p2 = np.zeros((3,num_threads))
    cdef double[:,:] p3 = np.zeros((3,num_threads))
    
    cdef long[:] s1 = np.zeros(num_threads,dtype=int)
    cdef long[:] s2 = np.zeros(num_threads,dtype=int)
    cdef long[:] s3 = np.zeros(num_threads,dtype=int)
    cdef long[:] s4 = np.zeros(num_threads,dtype=int)
    cdef long[:] s5 = np.zeros(num_threads,dtype=int)

    #memory views for parallel sign calculation 
    cdef double[:,:,:] vv = np.zeros((3,3,num_threads),dtype=float) # matrix
    cdef double[:,:] v1v = np.zeros((3,num_threads)) # comparison vector 
    cdef double[:,:] v2v = np.zeros((3,num_threads)) # comparison vector 
    cdef double[:,:] v3v = np.zeros((3,num_threads)) # comparison vector 
    cdef double[:,:] Sv = np.zeros((3,num_threads)) # cross product 
    

    # extract surface cells 
    with nogil, parallel(num_threads=num_threads):
        tid = openmp.omp_get_thread_num()
        for i in prange(nfaces):
            for j in range(3):
                x[j,tid] = node_xv[fconnection[i,j]]
                y[j,tid] = node_yv[fconnection[i,j]]
                z[j,tid] = node_zv[fconnection[i,j]]
            
            xm[tid] = (x[0,tid] + x[1,tid] + x[2,tid])/3
            ym[tid] = (y[0,tid] + y[1,tid] + y[2,tid])/3
            
            bqz[tid] = cellcentres[i,2]
            q1[0,tid] = xm[tid]; q1[1,tid] = ym[tid]; q1[2,tid] = tqz # top query point 
            q2[0,tid] = xm[tid]; q2[1,tid] = ym[tid]; q2[2,tid] = bqz[tid] #bottom query point 
            
            p1[0,tid] = x[0,tid]; p1[1,tid] = y[0,tid];  p1[2,tid] = z[0,tid]
            p2[0,tid] = x[1,tid]; p2[1,tid] = y[1,tid];  p2[2,tid] = z[1,tid]
            p3[0,tid] = x[2,tid]; p3[1,tid] = y[2,tid];  p3[2,tid] = z[2,tid]
            
            s1[tid] = tetra_signp(q1[:,tid],p1[:,tid],p2[:,tid],p3[:,tid],
                                  vv[:,:,tid], v1v[:,tid], v2v[:,tid],
                                  v3v[:,tid], Sv[:,tid])
            s2[tid] = tetra_signp(q2[:,tid],p1[:,tid],p2[:,tid],p3[:,tid],
                                  vv[:,:,tid], v1v[:,tid], v2v[:,tid],
                                  v3v[:,tid], Sv[:,tid])
            
            if s1[tid] != s2[tid]:
                s3[tid] = tetra_signp(q1[:,tid],q2[:,tid],p1[:,tid],p2[:,tid],
                                      vv[:,:,tid], v1v[:,tid], v2v[:,tid],
                                      v3v[:,tid], Sv[:,tid])
                s4[tid] = tetra_signp(q1[:,tid],q2[:,tid],p2[:,tid],p3[:,tid],
                                      vv[:,:,tid], v1v[:,tid], v2v[:,tid],
                                      v3v[:,tid], Sv[:,tid])
                s5[tid] = tetra_signp(q1[:,tid],q2[:,tid],p3[:,tid],p1[:,tid],
                                      vv[:,:,tid], v1v[:,tid], v2v[:,tid],
                                      v3v[:,tid], Sv[:,tid])
                if s3[tid] == s4[tid] and s4[tid] == s5[tid]:
                    ocheckv[i] = 1 
                
    return ocheck