cimport cython  #import relevant modules 
import math as ma

# mesh calculations - need to be in cython otherwise they will take too long 
@cython.boundscheck(False)#speeds up indexing, however it can be dangerous if the indexes are not in the range of the actual array, 
@cython.wraparound(False)
@cython.nonecheck(False)                      
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
    
def tri_combo(tuple con_mat):
    """Compute node combination for tetrahedra
    -------------
    con_mat: tuple of length 4
        Each entry in tuple is a list column on the connection matrix in the 
        3D mesh. 
    
    Returns
    --------------
    tri_combo: tuple
        Corresponding node combination indexes for each face of the cells in the 
        connection matrix
    """
    
    cdef int i,idx1,idx2,idx3,idx4
    cdef int num_elem = len(con_mat[0])
    cdef list face1s, face2s, face3s, face4s
    cdef str face1t, face2t, face3t, face4t
    cdef long long int face1, face2, face3, face4
    cdef tuple tri_combo = ([0]*num_elem,[0]*num_elem,[0]*num_elem,[0]*num_elem) # allocate space for tri_combo 
    
    for i in range(num_elem):
        idx1 = con_mat[0][i]#extract indexes 
        idx2 = con_mat[1][i]
        idx3 = con_mat[2][i]
        idx4 = con_mat[3][i]
        
        # assign each face an organised and unique code
        #sort face indexes 
        face1s = sorted((idx2,idx3,idx4))
        face2s = sorted((idx1,idx4,idx3))
        face3s = sorted((idx1,idx2,idx4))
        face4s = sorted((idx1,idx2,idx3))
        face1t = str(face1s[0])+str(face1s[1])+str(face1s[2])
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
    
    return tri_combo 

def group(tuple neigh_matrix, list groupings):
    """Group elements according to thier neighbours 
    """
    cdef int i, pt, a, j, h
    cdef int loops = len(groupings) # number of elements to loop through 
    cdef int dims = len(neigh_matrix) # dimensions of neighbour matrix
    cdef list param = [0]*loops
    
    pt = 0 
    for i in range(loops):
        if param[i]==0: # then the element is ungrouped
            pt += 1
            param[i] = pt
            a = 0 # number of parameters assigned per element
            for j in range(dims):
                h = neigh_matrix[j][i]
                if h!=-1 and param[h]==0 and a<=groupings[i]: # make sure its not already paired or on the edge
                    a += 1
                    param[h] = pt
                    
    return param


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