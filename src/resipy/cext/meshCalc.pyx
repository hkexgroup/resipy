cimport cython  #import relevant modules 
import math as ma
import numpy as np
cimport numpy as np
from cython.parallel import prange

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
def neigh3d(tuple con_mat, int return_tri_combo):
    """Compute neighbours of each element within a 3D tetrahedral mesh 
    Parameters
    -------------
    con_mat: tuple of length 4
        Each entry in tuple is a list column on the connection matrix in the 
        3D mesh. 
    return_tri_combo: int
        Binary, must be zero or 1. 1 to return the node combination matrix 
    
    Returns
    --------------
    neigh: tuple
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    """
    
    cdef int i 
    cdef str idx1,idx2,idx3,idx4
    cdef int num_elem = len(con_mat[0])
    cdef list face1s, face2s, face3s, face4s
    cdef str face1t, face2t, face3t, face4t
    cdef long long int face1, face2, face3, face4
    cdef tuple tri_combo = ([0]*num_elem,[0]*num_elem,[0]*num_elem,[0]*num_elem) # allocate space for tri_combo 
    cdef int num_node = max([max(con_mat[0]),max(con_mat[1]),max(con_mat[2]),max(con_mat[3])])
    cdef int pad = len(str(num_node))
    
    for i in range(num_elem):
        idx1 = int2str(con_mat[0][i],pad)#extract indexes 
        idx2 = int2str(con_mat[1][i],pad)
        idx3 = int2str(con_mat[2][i],pad)
        idx4 = int2str(con_mat[3][i],pad)
        
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
        
    cdef tuple neigh = ([0]*num_elem,[0]*num_elem,[0]*num_elem,[0]*num_elem) # allocate space for neighbour matrix              
    
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
    lookup = [i for i in range(num_elem)]*4 # index lookup 
    lookup_idx = [lookup[i] for i in temp_idx] # sorted index lookup
    
    for i in range(num_elem):
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
            neigh[j][i] = idx
        
    if return_tri_combo==1:
        return neigh,tri_combo
    else:
        return neigh

@cython.boundscheck(False)
def neigh2d(tuple con_mat, int return_tri_combo=0):
    """Compute neighbours of each element within a 2D triangular mesh 
    Parameters
    -------------
    con_mat: tuple of length 3
        Each entry in tuple is a list column on the connection matrix in the 
        3D mesh.  
    return_tri_combo: int
        Binary, must be zero or 1. 1 to return the node combination matrix 
    
    Returns
    --------------
    neigh: tuple
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    """
    
    cdef int i
    cdef str idx1,idx2,idx3,idx4
    cdef int num_elem = len(con_mat[0])
    cdef list face1s, face2s, face3s
    cdef str face1t, face2t, face3t
    cdef long long int face1, face2, face3
    cdef tuple tri_combo = ([0]*num_elem,[0]*num_elem,[0]*num_elem) # allocate space for tri_combo 
    cdef int num_node = max([max(con_mat[0]),max(con_mat[1]),max(con_mat[2])])
    cdef int pad = len(str(num_node))
    
    for i in range(num_elem):
        idx1 = int2str(con_mat[0][i],pad) # extract indexes 
        idx2 = int2str(con_mat[1][i],pad) 
        idx3 = int2str(con_mat[2][i],pad)
        
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
        
    cdef tuple neigh = ([0]*num_elem,[0]*num_elem,[0]*num_elem) # allocate space for neighbour matrix              
    
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
    lookup = [i for i in range(num_elem)]*4 # index lookup 
    lookup_idx = [lookup[i] for i in temp_idx] # sorted index lookup
    
    for i in range(num_elem):
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
            neigh[j][i] = idx
        
    if return_tri_combo==1:
        return neigh,tri_combo
    else:
        return neigh
    
@cython.boundscheck(False)
def split_tri(list con_mat, list node_x, list node_y, list node_z):
    """Split triangle elements down into 4 smaller triangles 
    Parameters
    -----------
    con_mat: list
        3 by N lists of the node numbers forming the triangle mesh
    node_x: list 
        Node x coordinates 
    node_y: list 
        Node x coordinates 
    node_z: list 
        Node z coordinates 
        
    Returns
    -----------
    new_con_mat: list
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
    cdef int num_elms = len(con_mat[0])
    cdef list nnodex 
    cdef list nnodey # new node lists 
    cdef list nnodez 
    cdef list new_con_mat = [[0]*num_elms*4,[0]*num_elms*4,[0]*num_elms*4] # new connection matrix 
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
            nno = con_mat[j][i] # node number
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
            nno = con_mat[j][i] # node number
            n[j] = nno 
        
        for j in range(3):
            nodes = sorted([int2str(n[a[j]],pad), int2str(n[b[j]],pad)])
            mn = int(nodes[0]+nodes[1])
            search = bisection_search(unicl,mn)
            
            n[j+3] = node_id[search]#reference index
    
        tmpi = i*4 # temporary index for indexing new connection matrix
        for j in range(4):
            new_con_mat[0][tmpi+j] = n[remap[j][0]]
            new_con_mat[1][tmpi+j] = n[remap[j][1]]
            new_con_mat[2][tmpi+j] = n[remap[j][2]]
            
    nnodex = [new_map[0][i] for i in idx]
    nnodey = [new_map[1][i] for i in idx]
    nnodez = [new_map[2][i] for i in idx]
                
    cdef list outx = node_x+nnodex
    cdef list outy = node_y+nnodey
    cdef list outz = node_z+nnodez
                
    return new_con_mat, outx, outy, outz, len(new_con_mat[0]), len(node_x)+len(nnodex)

@cython.boundscheck(False)
def split_tetra(list con_mat, list node_x, list node_y, list node_z):
    """Split tetrahedral elements down into 8 smaller tetrahedra 
    Parameters
    -----------
    con_mat: list
        3 by N lists of the node numbers forming the triangle mesh
    node_x: list 
        Node x coordinates 
    node_y: list 
        Node x coordinates 
    node_z: list 
        Node z coordinates 
        
    Returns
    -----------
    new_con_mat: list
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
    cdef int num_elms = len(con_mat[0])

    cdef list nnodex = [] 
    cdef list nnodey = [] # new node lists 
    cdef list nnodez = [] 

    cdef list new_con_mat = [[0]*num_elms*8,[0]*num_elms*8,[0]*num_elms*8,[0]*num_elms*8] # new connection matrix 
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
            nno = con_mat[j][i] # node number
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
            nno = con_mat[j][i] # node number
            n[j] = nno 
        
        for j in range(6):
            nodes = sorted([int2str(n[a[j]],pad), int2str(n[b[j]],pad)])
            mn = int(nodes[0]+nodes[1])
            search = bisection_search(unicl,mn)
            
            n[j+4] = node_id[search]#reference index
    
        tmpi = i*8 # temporary index for indexing new connection matrix
        for j in range(8):
            new_con_mat[0][tmpi+j] = n[remap[j][0]]
            new_con_mat[1][tmpi+j] = n[remap[j][1]]
            new_con_mat[2][tmpi+j] = n[remap[j][2]]
            new_con_mat[3][tmpi+j] = n[remap[j][3]]
            
    nnodex = [new_map[0][i] for i in idx]
    nnodey = [new_map[1][i] for i in idx]
    nnodez = [new_map[2][i] for i in idx]
                
    cdef list outx = node_x+nnodex
    cdef list outy = node_y+nnodey
    cdef list outz = node_z+nnodez
                    
    return new_con_mat, outx, outy, outz, len(new_con_mat[0]), len(node_x)+len(nnodex)

@cython.boundscheck(False)
@cython.wraparound(False)   
def orderTetra(long[:,:] connection, double[:,:] node):
    """ Organise tetrahedral element nodes into a clockwise order
    
    following solution at: 
    https://stackoverflow.com/questions/10612829/tetrahedron-orientation-for-triangle-meshes
    
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
    count : TYPE
        Number of elements with switched nodes (ie been 'corrected')
    """
    #get the number of nodes and elements etc 
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    cdef int count = 0 #rolling total for the number of corrected elements 
    
    con = np.zeros((numel,npere), dtype=DINT)# new connection matrix 
    cdef int[:,:] conv = con # connection memeory view
    
    cdef np.ndarray[double, ndim=1] node_x = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef np.ndarray[double, ndim=1] node_y = np.asarray(node[:,1],dtype=float)
    cdef np.ndarray[double, ndim=1] node_z = np.asarray(node[:,2],dtype=float)

    #looping variables 
    cdef np.ndarray[double,ndim=2] v = np.zeros((3,3),dtype=float) # matrix
    cdef np.ndarray[double,ndim=1] v1 = np.zeros(3) # comparison vector 
    cdef np.ndarray[double,ndim=1] v2 = np.zeros(3)
    cdef np.ndarray[double,ndim=1] v3 = np.zeros(3)
    cdef np.ndarray[double,ndim=1] S = np.zeros(3)
    cdef double N # orientation indicator 
    cdef np.ndarray[long, ndim=1] ccw = np.zeros(numel,dtype=int) # clockwise array 
    cdef int i, k, ei # loop integers 
    cdef int num_threads = 2 # number of threads to use (using 2 for now)

    for i in prange(numel, nogil=True, num_threads=num_threads,schedule='static'):
        for k in range(3): # work out delta matrix 
            v[0,k] = node_x[connection[i,k+1]] - node_x[connection[i,0]] 
            v[1,k] = node_y[connection[i,k+1]] - node_y[connection[i,0]] 
            v[2,k] = node_z[connection[i,0]] - node_z[connection[i,k+1]] 
            
        for k in range(3): #extract delta vectors 
            v1[k] = v[k,0]
            v2[k] = v[k,1]
            v3[k] = v[k,2]
        
        #compute cross product
        S[0] = v1[1]*v2[2] - v1[2]*v2[1]
        S[1] = v1[2]*v2[0] - v1[0]*v2[2]
        S[2] = v1[0]*v2[1] - v1[1]*v2[0]
        #compute dot product (gives orientation)
        N = S[0]*v3[0] + S[1]*v3[1] + S[2]*v3[2]
        
        if N>0: # then tetrahedra is clockwise 
            count += 1
            conv[i,1] = connection[i,0]
            conv[i,0] = connection[i,1]
            ccw[i] = 1
        elif N<0: # then its counter clockwise 
            conv[i,0] = connection[i,0]
            conv[i,1] = connection[i,1]
            ccw[i] = 2
        else: # shouldn't be possible 
            ccw[i] = 0
            ei = i #problem index 
            
        conv[i,2] = connection[i,2]
        conv[i,3] = connection[i,3]
    
    for i in range(numel):
        if ccw[i]==0 and ei>-1:
            raise ValueError('Element %i has all colinear nodes, thus is poorly formed and can not be ordered'%ei)
            
    return con, count, ccw

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
    cdef int count = 0 # number of corrected elements 
    
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
    
    cdef np.ndarray[double, ndim=1] theta = np.zeros(4,dtype=float) # theta array
    cdef double[:] thetav = theta # theta view 
    
    cdef np.ndarray[long, ndim=1] order = np.zeros(4,dtype=int) #ordering array
    cdef long[:] orderv = order #ordering view
    
    cdef int i, k, c, j
    cdef int eflag = 0 # error flag
    cdef int ei # error number 
    cdef int num_threads = 2
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
    #for i in prange(numel, nogil=True, num_threads=num_threads):
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
                count += 1 
                break
        
        for k in range(4):
            conv[i,k] = connection[i,order[k]]
            
        #nb: looked at parallising this loop but it yeilded unstable results 
            
        
    if eflag == 1: 
        raise ValueError('Element %i has more than 2 colinear points and therefore not a quad'%ei)
        
    return con, count