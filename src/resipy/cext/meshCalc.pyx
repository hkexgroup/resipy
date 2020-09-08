cimport cython  #import relevant modules 
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
cdef int bisection_search(long[:] arr, long long int var) nogil:
    """Efficent search algorithm for sorted array of postive ints 
    (memory veiw parallel capable version)
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
        m = int((L+R)/2)
        if arr[m]<var:
            L = m+1
        elif arr[m]>var:
            R = m-1
        else:
            return m
    return -1

@cython.boundscheck(False)
@cython.wraparound(False)  
def bisection_searchL(list arr, long long int var):
    """Efficent search algorithm for sorted list of postive ints 
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
        m = int((L+R)/2)
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
        search = bisection_searchL(uni,arrs[i])
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
cdef void sortInt(long[:] arr, int n) nogil: 
    #adapted from https://www.geeksforgeeks.org/insertion-sort/
    #sorts array in place
    cdef int i, j, key
    for i in range(1, n): 
        key = arr[i] 
        j = i-1
        while j >= 0 and key < arr[j] : 
                arr[j + 1] = arr[j] 
                j -= 1
        arr[j + 1] = key 
    #return arr

cdef long long int mergeInts(int a, int b, int c, int pad) nogil: # merge 3 ints 
    cdef long long int x = a*10**pad + b # merge a and b
    return x*10**pad + c # merge c 

cdef long long int mergeInt(int a, int b, int pad) nogil: #merge 2 ints 
    return a*10**pad + b # merge a and b

@cython.boundscheck(False)    
@cython.wraparound(False)             
def neigh3d(long[:,:] connection, int return_tri_combo, int num_threads=2):
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
    
    cdef int i, j, tid #thread indexing 
    cdef int numel = connection.shape[0]
    #face arrays 
    cdef long[:] face1s = np.zeros(3,dtype=int)
    cdef long[:] face2s = np.zeros(3,dtype=int)
    cdef long[:] face3s = np.zeros(3,dtype=int)
    cdef long[:] face4s = np.zeros(3,dtype=int)
    #cdef long[:,:] face np.zeros((4,num_threads),dtype=int)
    #combination array 
    cdef np.ndarray[long, ndim=2] tri_combo = np.zeros((numel,4),dtype=int,order='C') # allocate space for tri_combo 
    cdef long[:,:] tri_combov = tri_combo
    #maximum number of node digits 
    cdef int num_node = max([max(connection[:,0]),max(connection[:,1]),
                             max(connection[:,2]),max(connection[:,3])])
    cdef int pad = len(str(num_node))
    
    #assign unique node combinations to each element edge, repeats pick out neighbours
    for i in range(numel):
        face1s[0] = connection[i,1]; face1s[1] = connection[i,2]; face1s[2] = connection[i,3]
        face2s[0] = connection[i,0]; face2s[1] = connection[i,3]; face2s[2] = connection[i,2]
        face3s[0] = connection[i,0]; face3s[1] = connection[i,1]; face3s[2] = connection[i,3]
        face4s[0] = connection[i,0]; face4s[1] = connection[i,1]; face4s[2] = connection[i,2]
        #sort face indexes 
        sortInt(face1s,3)
        sortInt(face2s,3)
        sortInt(face3s,3)
        sortInt(face4s,3)
        # assign each face an organised and unique code
        tri_combov[i,0] = mergeInts(face1s[0],face1s[1],face1s[2],pad) #face 1 
        tri_combov[i,1] = mergeInts(face2s[0],face2s[1],face2s[2],pad) #face 2 
        tri_combov[i,2] = mergeInts(face3s[0],face3s[1],face3s[2],pad) #face 3 
        tri_combov[i,3] = mergeInts(face4s[0],face4s[1],face4s[2],pad) #face 4

    cdef np.ndarray[long, ndim=2] neigh = np.zeros((numel,4),dtype=int) # allocate space for neighbour matrix        
    cdef long[:,:] neighv = neigh  
    
    #using binary search and sorted lists for efficient index lookup 
    cdef np.ndarray[long, ndim=1] tri_flatten = tri_combo.T.flatten() #all faces in one array 
	
    cdef np.ndarray[long, ndim=1] temp_idx = np.argsort(tri_flatten).astype(int,order='C') # sorted indexes  (pure python code)
    cdef np.ndarray[long, ndim=1] tri_sort = tri_flatten[temp_idx] # sorted faces (forward direction)
    cdef long[:] tri_sortv = tri_sort
    cdef int maxo = len(tri_sort)-1

    cdef long long o, idx
    #create lookup array 
    cdef np.ndarray[long, ndim=1] lookup = np.zeros((numel*4),dtype=int)
    for i in range(numel):
        lookup[i] = i
        lookup[i+numel] = i
        lookup[i+(2*numel)] = i
        lookup[i+(3*numel)] = i

    cdef np.ndarray[long, ndim=1] lookup_idx = np.asarray(lookup[temp_idx],order='C')
    cdef long[:] lookup_idxv = lookup_idx
  
    #loop is parallel becuase it can be intense
    for i in prange(numel,nogil=True,num_threads=num_threads,schedule='dynamic'):
        for j in range(4):
            o = bisection_search(tri_sortv,tri_combov[i,j]) # only works on a sorted array
            #find the reference index
            if lookup_idxv[o]==i:# then there are 2 options 
                #1 ) the face in question is unique or
                #2 ) the index is o+1 or o-1
                if o != maxo and tri_sortv[o+1] == tri_combov[i,j]:
                    o = o+1
                elif o !=0 and tri_sortv[o-1] == tri_combov[i,j]:
                    o = o-1
                else: # unique face on edge of mesh 
                    o = -1   
            
            if o == -1: 
                idx = -1  
            else:
                idx = lookup_idxv[o]
            neighv[i,j] = idx
    
    if return_tri_combo==1:
        return neigh,tri_combo
    else:
        return neigh
    
@cython.boundscheck(False)
@cython.wraparound(False)
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
        
    return fconnection, idx

@cython.boundscheck(False)
@cython.wraparound(False)
def neigh2d(long[:,:] connection, int return_tri_combo=1):
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
    cdef int i,j
    cdef str idx1,idx2,idx3,idx4
    cdef int numel = connection.shape[0]
    #face arrays 
    cdef long[:] face1s = np.zeros(2,dtype=int)
    cdef long[:] face2s = np.zeros(2,dtype=int)
    cdef long[:] face3s = np.zeros(2,dtype=int)
    #combination array 
    cdef np.ndarray[long, ndim=2] tri_combo = np.zeros((numel,3),dtype=int) # allocate space for tri_combo 
    cdef long[:,:] tri_combov = tri_combo
    cdef int num_node = max([max(connection[:,0]),max(connection[:,1]),
                             max(connection[:,2])])
    cdef int pad = len(str(num_node))
    
    for i in range(numel):
        face1s[0] = connection[i,0]; face1s[1] = connection[i,1]; 
        face2s[0] = connection[i,1]; face2s[1] = connection[i,2]; 
        face3s[0] = connection[i,2]; face3s[1] = connection[i,0]; 
        #sort face indexes 
        sortInt(face1s,2)
        sortInt(face2s,2)
        sortInt(face3s,2)
        # assign each face an organised and unique code
        tri_combov[i,0] = mergeInt(face1s[0],face1s[1],pad) #face 1 
        tri_combov[i,1] = mergeInt(face2s[0],face2s[1],pad) #face 2 
        tri_combov[i,2] = mergeInt(face3s[0],face3s[1],pad) #face 3 
        
    cdef np.ndarray[long, ndim=2] neigh = np.zeros((numel,3),dtype=int) # allocate space for neighbour matrix        
    cdef long[:,:] neighv = neigh             
    
    #using binary search and sorted arrays for efficient index lookup 
    #cdef list tri_list = tri_combo[0] + tri_combo[1] + tri_combo[2] + tri_combo[3] #all faces together 
    cdef np.ndarray[long, ndim=1] tri_flatten = tri_combo.T.flatten()
	
    cdef np.ndarray[long, ndim=1] temp_idx = np.argsort(tri_flatten).astype(int,order='C') # sorted indexes  (pure python code)
    cdef np.ndarray[long, ndim=1] tri_sort = tri_flatten[temp_idx] # sorted faces (forward direction)
    cdef int maxo = len(tri_sort)-1
    cdef long long o
    cdef int idx
    
    cdef np.ndarray[long, ndim=1] lookup = np.zeros((numel*4),dtype=int)
    for i in range(numel):
        lookup[i] = i
        lookup[i+numel] = i
        lookup[i+(2*numel)] = i
        lookup[i+(3*numel)] = i

    cdef np.ndarray[long, ndim=1] lookup_idx = lookup[temp_idx]
    
    for i in range(numel):
        for j in range(3):
            o = bisection_search(tri_sort,tri_combo[i,j]) # only works on a sorted list 
            #find the reference index
            if lookup_idx[o]==i:# then there are 2 options 
                #1 ) the face in question is unique or
                #2 ) the index is o+1 or o-1
                if o!=maxo and tri_sort[o+1] == tri_combo[i,j]:
                    o+=1
                elif o!=0 and tri_sort[o-1] == tri_combo[i,j]:
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
@cython.wraparound(False)
def split_tri(long[:,:] connection, double[:,:] node):
    """Split triangle elements down into 4 smaller triangles 
    Parameters
    -----------
    connection: np array 
        3 by N array of int describing indices of element vertices (nodes)
    node: np array 
        3 by N array of node coordinates 
        
    Returns
    -----------
    new_connection: np array 
    node: np array 
    numel: int
    num_nodes: int
    """
    #define variables 
    cdef double[:] node_x = np.asarray(node[:,0], dtype=float)
    cdef double[:] node_y = np.asarray(node[:,1], dtype=float)
    cdef double[:] node_z = np.asarray(node[:,2], dtype=float)
    
    #loop variables 
    cdef double[:] x = np.zeros(3,dtype=float)
    cdef double[:] y = np.zeros(3,dtype=float)
    cdef double[:] z = np.zeros(3,dtype=float)
    cdef long[:] n = np.zeros(6,dtype=int) # new node numbers 
    cdef float mx, my, mz
    cdef long[:] nodes = np.zeros(2,dtype=int)
    cdef unsigned long long int mn
    cdef int i, j, nno, tmpi 

    #matrices 
    cdef int num_nodes = len(node_x) # number of nodes 
    cdef int numel = connection.shape[0] # number of elements 
    
    cdef np.ndarray[long, ndim=2] new_connection = np.zeros((numel*4,3),dtype = int) # new connection matrix 
    cdef long[:,:] new_connectionv = new_connection # memory view 
    #the columns are as follows: nodex,nodey,nodez,node_config, original element no
    cdef np.ndarray[double, ndim=2] new_node = np.zeros((numel*3,3),dtype=float)
    #cdef np.ndarray[long, ndim=1] new_node_idx = np.zeros(numel*3,dtype=int)
    cdef list new_node_idx = [0]*numel*3
    cdef np.ndarray[long, ndim=1] og_el_id = np.zeros(numel*3,dtype=int)
    cdef double[:,:] new_nodev = new_node
    # cdef long[:] new_node_idxv = new_node_idx
    cdef long[:] og_el_idv = og_el_id 
    # cdef list new_map = [[0]*numel*3,[0]*numel*3,[0]*numel*3,[0]*numel*3,[0]*numel*3] # used to reference already computed node centres
    cdef int pad = len(str(num_nodes))+1
    
    cdef long[:] a = np.array([0,1,2],dtype=int)
    cdef long[:] b = np.array([1,2,0],dtype=int)
    
    #remap the element nodes with this matrix 
    cdef list remap = [[0, 5, 3],
                       [1, 4, 3],
                       [2, 5, 4],
                       [3, 4, 5]]
    cdef long[:,:] remapv = np.asarray(remap,dtype=int)
    
    cdef list unicl, idx, node_id 
    
    ### loop through the elements ### all possible node configurations 
    for i in range(numel):
        for j in range(3):
            nno = connection[i,j] # node number
            x[j] = node_x[nno]
            y[j] = node_y[nno]
            z[j] = node_z[nno]
            n[j] = nno 
            
        tmpi = 3*i
        for j in range(3):
            mx = (x[a[j]]+x[b[j]])/2
            my = (y[a[j]]+y[b[j]])/2
            mz = (z[a[j]]+z[b[j]])/2

            nodes[0] = n[a[j]]; nodes[1] = n[b[j]]
            sortInt(nodes,2)
            mn = mergeInt(nodes[0],nodes[1],pad)
            
            new_nodev[tmpi+j,0] = mx # mid point of 2 respective nodes 
            new_nodev[tmpi+j,1] = my
            new_nodev[tmpi+j,2] = mz
            new_node_idx[tmpi+j] = mn # node combination
            og_el_idv[tmpi+j] = i
    
    ### find unique node configs # ###       
    unicl,idx = unique(new_node_idx) # unique and ordered node configs 
    node_id = [num_nodes + i for i in range(len(unicl))]
    cdef np.ndarray[long, ndim=1] unicla = np.asarray(unicl,dtype=int)
    cdef np.ndarray[long, ndim=1] idxa = np.asarray(idx,dtype=int)
    
    ### map back to elements #### 
    for i in range(numel):
        for j in range(3):
            nno = connection[i,j] # node number
            n[j] = nno 
        
        for j in range(3):
            nodes[0] = n[a[j]]; nodes[1] = n[b[j]]
            sortInt(nodes,2)
            mn = mergeInt(nodes[0],nodes[1],pad)
            search = bisection_searchL(unicl,mn)
            
            n[j+3] = node_id[search]#reference index
    
        tmpi = i*4 # temporary index for indexing new connection matrix
        for j in range(4):
            new_connectionv[tmpi+j,0] = n[remap[j,0]]
            new_connectionv[tmpi+j,1] = n[remap[j,1]]
            new_connectionv[tmpi+j,2] = n[remap[j,2]]
            
    ### make new node matrix ### 
    cdef int added_nodes = len(idx)
    cdef np.ndarray[double, ndim=2] node_out = np.zeros((num_nodes+added_nodes,3),dtype=float)
    cdef double[:,:] node_outv = node_out
    for i in range(num_nodes):
        node_outv[i,0] = node_x[i]
        node_outv[i,1] = node_y[i]
        node_outv[i,2] = node_z[i]
    for i in range(added_nodes):
        j = i + num_nodes
        node_outv[j,0] = new_nodev[idxa[i],0]
        node_outv[j,1] = new_nodev[idxa[i],1]
        node_outv[j,2] = new_nodev[idxa[i],2]
                
    return new_connection, node_out, numel*4, num_nodes + added_nodes

@cython.boundscheck(False)
@cython.wraparound(False)
def split_tetra(long[:,:] connection, double[:,:] node):
    """Split tetrahedral elements down into 8 smaller tetrahedra 
    Parameters
    -----------
    connection: np array 
        4 by N array of int describing indices of element vertices (nodes)
    node: np array 
        3 by N array of node coordinates 
        
    Returns
    -----------
    new_connection: np array 
    node: np array 
    numel: int
    num_nodes: int
    """
    cdef double[:] node_x = np.asarray(node[:,0], dtype=float)
    cdef double[:] node_y = np.asarray(node[:,1], dtype=float)
    cdef double[:] node_z = np.asarray(node[:,2], dtype=float)
    cdef double[:] x = np.zeros(4,dtype=float)
    cdef double[:] y = np.zeros(4,dtype=float)
    cdef double[:] z = np.zeros(4,dtype=float)
    cdef long[:] n = np.zeros(10,dtype=int) # new node numbers 
    cdef double mx, my, mz
    cdef long[:] nodes = np.zeros(2,dtype=int)
    cdef long long int mn
    cdef int i, j, nno, tmpi #loop variables 

    cdef int num_nodes = len(node_x)
    cdef int numel = connection.shape[0]

    cdef np.ndarray[long, ndim=2] new_connection = np.zeros((numel*8,4),dtype=int) # new connection matrix 
    cdef long[:,:] new_connectionv = new_connection
    #the columns are as follows: nodex,nodey,nodez,node_config, original element no
    #cdef list new_map = [[0]*numel*6,[0]*numel*6,[0]*numel*6,[0]*numel*6,[0]*numel*6] # used to reference already computed node centres
    cdef np.ndarray[double, ndim=2] new_node = np.zeros((numel*6,3),dtype=float)
    # cdef np.ndarray[long, ndim=1] new_node_idx = np.zeros(numel*6,dtype=int)
    cdef list new_node_idx = [0]*numel*6
    cdef np.ndarray[long, ndim=1] og_el_id = np.zeros(numel*6,dtype=int)
    cdef double[:,:] new_nodev = new_node
    # cdef long[:] new_node_idxv = new_node_idx
    cdef long[:] og_el_idv = og_el_id 
    
    cdef int pad = len(str(num_nodes))
    
    cdef long[:] a = np.array([0,1,2,0,1,2],dtype=int)
    cdef long[:] b = np.array([1,2,0,3,3,3],dtype=int)
    
    #remap the element nodes with this matrix 
    cdef list remap = [[4, 7, 6, 0],
                       [5, 9, 6, 2],
                       [8, 9, 7, 3],
                       [8, 9, 7, 6],
                       [8, 4, 7, 6],
                       [8, 5, 9, 6],
                       [4, 5, 6, 8],
                       [4, 1, 5, 8]]
    cdef long[:,:] remapv = np.asarray(remap,dtype=int)
    
    cdef list unicl, idx, node_id 
    
    ### loop through the elements ### all possible node configurations 
    for i in range(numel):
        for j in range(4):
            nno = connection[i,j] # node number
            x[j] = node_x[nno]
            y[j] = node_y[nno]
            z[j] = node_z[nno]
            n[j] = nno 
            
        tmpi = 6*i
        for j in range(6):
            mx = (x[a[j]]+x[b[j]])/2
            my = (y[a[j]]+y[b[j]])/2
            mz = (z[a[j]]+z[b[j]])/2
            nodes[0] = n[a[j]]; nodes[1] = n[b[j]]
            sortInt(nodes,2)
            mn = mergeInt(nodes[0],nodes[1],pad)
            
            new_nodev[tmpi+j,0] = mx
            new_nodev[tmpi+j,1] = my
            new_nodev[tmpi+j,2] = mz
            new_node_idx[tmpi+j] = mn
            og_el_idv[tmpi+j] = i
    
    ### find unique node configs # ###       
    unicl,idx = unique(new_node_idx) # unique and ordered node configs 
    cdef int added_nodes = len(idx) # number of new nodes added to mesh 
    node_id = [num_nodes + i for i in range(len(unicl))]
    cdef np.ndarray[long, ndim=1] unicla = np.asarray(unicl,dtype=int)
    cdef np.ndarray[long, ndim=1] idxa = np.asarray(idx,dtype=int)
    
    ### map back to elements #### 
    for i in range(numel):
        for j in range(4):
            nno = connection[i,j] # node number
            n[j] = nno 
        
        for j in range(6):
            nodes[0] = n[a[j]]; nodes[1] = n[b[j]]
            sortInt(nodes,2)
            mn = mergeInt(nodes[0],nodes[1],pad)
            search = bisection_searchL(unicl,mn)
            n[j+4] = node_id[search]#reference index
    
        tmpi = i*8 # temporary index for indexing new connection matrix
        for j in range(8):
            new_connectionv[tmpi+j,0] = n[remap[j,0]
            new_connectionv[tmpi+j,1] = n[remap[j,1]
            new_connectionv[tmpi+j,2] = n[remap[j,2]
            new_connectionv[tmpi+j,3] = n[remap[j,3]
            
    ### make new node matrix ### 
    cdef np.ndarray[double, ndim=2] node_out = np.zeros((num_nodes+added_nodes,3),dtype=float)
    cdef double[:,:] node_outv = node_out
    for i in range(num_nodes):
        node_outv[i,0] = node_x[i]
        node_outv[i,1] = node_y[i]
        node_outv[i,2] = node_z[i]
    for i in range(added_nodes):
        j = i + num_nodes
        node_outv[j,0] = new_nodev[idxa[i],0]
        node_outv[j,1] = new_nodev[idxa[i],1]
        node_outv[j,2] = new_nodev[idxa[i],2]
                    
    return new_connection, node_out, numel*8, num_nodes + added_nodes

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
    
    cdef np.ndarray[long,ndim=2] con = np.zeros((numel,npere), dtype=int)# new connection matrix 
    cdef long[:,:] conv = con # connection memeory view
    
    cdef np.ndarray[double, ndim=1] node_x = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef np.ndarray[double, ndim=1] node_y = np.asarray(node[:,1],dtype=float)
    cdef np.ndarray[double, ndim=1] node_z = np.asarray(node[:,2],dtype=float)
    cdef double[:] nodex = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef double[:] nodey = np.asarray(node[:,1],dtype=float)
    cdef double[:] nodez = np.asarray(node[:,2],dtype=float)

    #looping variables 
    cdef double v00,v01,v02,v10,v11,v12,v20,v21,v22,s0,s1,s2, N 
    
    cdef np.ndarray[long, ndim=1] ccw = np.zeros(numel,dtype=int) # clockwise array 
    cdef long[:] ccwv = ccw #clockwise view
    cdef Py_ssize_t i 
    cdef int k, ei, tid # loop integers 
    cdef np.ndarray[long, ndim=1] count = np.zeros(numel,dtype=int)
    cdef long[:] countv = count
    cdef int count_out #rolling total for the number of corrected elements 

    #with nogil, parallel(num_threads=num_threads):
    for i in prange(numel,nogil=True,num_threads=num_threads,schedule='static'):
        v00 = nodex[connection[i,1]] - nodex[connection[i,0]] 
        v01 = nodex[connection[i,2]] - nodex[connection[i,0]] 
        v02 = nodex[connection[i,3]] - nodex[connection[i,0]] 
        v10 = nodey[connection[i,1]] - nodey[connection[i,0]] 
        v11 = nodey[connection[i,2]] - nodey[connection[i,0]] 
        v12 = nodey[connection[i,3]] - nodey[connection[i,0]] 
        v20 = nodez[connection[i,0]] - nodez[connection[i,1]]
        v21 = nodez[connection[i,0]] - nodez[connection[i,2]]
        v22 = nodez[connection[i,0]] - nodez[connection[i,3]]
        
        #compute cross product
        s0 = v01*v12 - v02*v11
        s1 = v02*v10 - v00*v12
        s2 = v00*v11 - v01*v10
        
        #compute dot product (gives orientation)
        N = s0*v20 + s1*v21 + s2*v22
        
        if N>0: # then tetrahedra is clockwise 
            countv[i] = 1
            conv[i,1] = connection[i,0]
            conv[i,0] = connection[i,1]
            ccwv[i] = 1
        elif N<0: # then its counter clockwise 
            conv[i,0] = connection[i,0]
            conv[i,1] = connection[i,1]
            ccwv[i] = 2
        else: # shouldn't be possible 
            ccwv[i] = 0
            ei = i #problem index 
            
        conv[i,2] = connection[i,2]
        conv[i,3] = connection[i,3]
        
    for i in range(numel):
        if ccwv[i]==0:
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
    #with nogil, parallel(num_threads=num_threads):
    for i in prange(nfaces,nogil=True,num_threads=num_threads,schedule='static'):
        tid = openmp.omp_get_thread_num()
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

def paramLookup(long[:] param, long[:,:] neigh):
    
    cdef int numel = param.shape[0]
    cdef int nneigh = neigh.shape[1] # number of nieghbours 
    cdef np.ndarray[long, ndim=2] out = np.zeros(neigh.shape,dtype=int)
    cdef long[:,:] outv = out # memory view 
    cdef int i,j #looping variables 
    
    for i in range(numel):
        for j in range(nneigh):
            outv[i,j] = param[neigh[i,j]]
    return out

#nsizeA and finite element conductance calculation
def conductanceCall(long[:,:] connection, int numnp, int typ=0,
                    int num_threads=1):
    """Calculate the array size needed for the finite element conductance matrix
    in R2 class codes 

    Parameters
    ----------
    connection: numpy array 
        mesh connection matrix.
    node: numpy array
        mesh node matrix.
    typ: vtk cell type
        DESCRIPTION. The default is 0. Will raise an error if takes an 
        unexpected value. 

    Returns
    -------
    nsizeA : int
        Number of connected nodes + number of nodes

    Note
    -----
    typ is 5 for triangle meshes, 8 or 9 for quads, 10 for tetrahedra, 13 for 
    a prism. 
    """
    cdef int nmax = 60
    cdef int numel = connection.shape[0]
    cdef int nedges 
    cdef int pad = len(str(numnp))
    cdef long[:] a, b
    cdef np.ndarray[long long, ndim=1] idx, counts, # unique counts and indices 
    cdef np.ndarray[long long, ndim=1] uni, combof # uni segment combinations 
    cdef np.ndarray[long, ndim=2] Nconnec = np.zeros((numnp,nmax),dtype=int) - 1
    cdef long[:,:] Nconnecv = Nconnec

    #determine number of edges 
    if typ==5:#then elements are triangles
        nedges = 3
        a = np.asarray([0,1,2], dtype=int)
        b = np.asarray([1,2,0], dtype=int)
    elif typ==8 or typ==9:#elements are quads
        nedges = 4
        a = np.asarray([0,1,2,3], dtype=int)
        b = np.asarray([1,2,3,0], dtype=int)
    elif typ == 10:# elements are tetrahedra 
        nedges = 6
        a = np.asarray([0,1,2,3,3,3], dtype=int)
        b = np.asarray([1,2,0,0,1,2], dtype=int)
    elif typ == 13: # elements are 3d wedges 
        nedges = 9
        a = np.asarray([0,1,2, 3,4,5, 0,1,2], dtype=int)
        b = np.asarray([1,2,0, 4,5,3, 3,4,5], dtype=int)
    else: 
        raise ValueError('Cell type argument does not match vtk cell types,', 
                         'used with meshTools must be one of the following:',
                         '5, 8, 9, 10 or 13')

    #looping variables     
    cdef long long[:,:] combo = np.zeros((numel,nedges),dtype=np.int64)
    cdef long[:] nodes = np.zeros(2,dtype=int)
    cdef int i,j,k
    cdef long long merged 
    cdef long na, nb, nid
    
    #find unique node combinations 
    with nogil, parallel(num_threads=num_threads):
        for i in prange(numel,schedule='dynamic', chunksize=1):
            for j in range(nedges):
                na = connection[i,a[j]]
                nb = connection[i,b[j]] 
                if na < nb:
                    merged = mergeInt(na,nb,pad) # merge
                else:
                    merged = mergeInt(nb,na,pad) # merge
                combo[i,j] = merged       
                
                #create the conductance matrix             
                for k in range(nmax):
                    nid= Nconnecv[na,k]
                    if nid == -1: # then its not been assigned as a pair
                        Nconnecv[na,k] = nb
                        break # break out the loop 
                    elif nid == nb: #its already been assigned 
                        break # so break out 
                    
                for k in range(nmax):#same as above but for the inverse
                    nid = Nconnecv[nb,k]
                    if nid == -1:
                        Nconnecv[nb,k] = na
                        break
                    elif nid == na:
                        break

    #pure python code 
    combof = np.asarray(combo,dtype=np.int64).flatten()
    uni, idx, counts = np.unique(combof, 
                                 return_index=True, 
                                 return_counts=True)
    if max(counts) > 60: 
        raise Exception('Too many connected nodes (>%i) check that mesh is not poorly formed.'%nmax)
        
    cdef int nsizeA = len(uni) + numnp # nsizeA can now be computed 
    
    return nsizeA, Nconnec