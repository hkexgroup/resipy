cimport cython 
import math as ma

# mesh calculations - need to be in cython otherwise they will take too long 
@cython.boundscheck(False)#speeds up indexing, however it can be dangerous if the indexes are not in the range of the actual array, 
                        #shouldnt be an issue in this code. 
@cython.wraparound(False)
                        
def binary_search(list arr, long long int var):
    """Efficent search algorithm for sorted lists
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
    return False
                        
def neigh3d(tuple con_mat):
    """Compute neighbours of each element within a 3D tetrahedral mesh 
    Parameters
    -------------
    con_mat: tuple of length 4
        Each entry in tuple is a list column on the connection matrix in the 
        3D mesh. 
    
    Returns
    --------------
    neigh: tuple
        Corresponding neighbour indexes for each face of the cells in the 
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
        
    cdef tuple neigh = ([0]*num_elem,[0]*num_elem,[0]*num_elem,[0]*num_elem) # allocate space for neighbour matrix              
    
    #using binary search and sorted lists for efficient index lookup 
    cdef list tri_list = tri_combo[0] + tri_combo[1] + tri_combo[2] + tri_combo[3] #all faces together 
    cdef list temp_idx 
    cdef list tri_sort
	
    temp_idx = [iterat[0] for iterat in sorted(enumerate(tri_list), key=lambda x:x[1])] # sorted indexes  
    tri_sort= [tri_list[i] for i in temp_idx] # sorted faces (forward direction)
    cdef int maxo = len(tri_sort)-1
    cdef long long int o
    
    cdef list lookup 
    cdef list lookup_idx 
    lookup = [i for i in range(num_elem)]*4 # index lookup 
    lookup_idx = [lookup[i] for i in temp_idx] # sorted index lookup
    
    for i in range(num_elem):
        for j in range(4):
            o = binary_search(tri_sort,tri_combo[j][i]) # only works on a sorted list 
            #find the reference index
            if lookup_idx[o]==i:# then there are 2 options 
                #1 ) the face in question is unique or
                #2 ) the index is o+1 or o-1
                if tri_sort[o+1] == tri_combo[j][i] and o!=maxo:
                    o+=1
                elif tri_sort[o-1] == tri_combo[j][i] and o!=0:
                    o-=1
                else: # unique face on edge of mesh 
                    o = -1   
            
            if o==-1: 
                idx = -1  
            else:
                idx = lookup_idx[o]
            neigh[j][i] = idx
        
    return neigh


