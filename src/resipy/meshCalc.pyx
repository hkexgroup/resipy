
# mesh calculations - need to be in cython otherwise they will take too dlong 
cimport cython  #import relevant modules 
import numpy as np
cimport numpy as np
ctypedef {:} dlong # 64bit dlong data type, plaform dependent, edit as necessary 
npy_long = np.int64 # what numpy considers a dlong data type, note that this not platform dependent 

cdef extern from "math.h" nogil:
    cpdef double acos(double x)
cdef extern from "math.h" nogil:
    cpdef double atan(double x)

@cython.boundscheck(False)#speeds up indexing, however it can be dangerous if the indexes are not in the range of the actual array, 
@cython.wraparound(False)        
cdef Py_ssize_t bisectionSearch(dlong [:] arr, dlong var):
    """Efficent search algorithm for sorted array of postive ints 
    (memory veiw parallel capable version)
    
    Parameters
    -------------
    arr: list 
        sorted list 
    var: int
        item to be searched / indexed. If not found False is returned 
    """
    cdef dlong L = 0
    cdef dlong n = len(arr)
    cdef dlong R = n-1
    cdef Py_ssize_t m
    while L <= R:
        m = int((L+R)/2)
        if arr[m]<var:
            L = m+1
        elif arr[m]>var:
            R = m-1
        else:
            return m
    return -1

# @cython.boundscheck(False) - function depreciated 
# def unique(list arr): # unique with out numpy for an array of ints 
#     """Find the unique values inside a list of ints (ONLY!!).Relies on the 
#     efficeint bisection search as above. 
    
#     Parameters 
#     -------------
#     arr: list 
#         an unordered list of ints 
    
#     Returns
#     ------------
#     uni: list
#         list of unique ints in arr
#     idx: list
#         list of the corresponding indexes for the unique values in arr
#     """
#     cdef int i 
#     cdef list temp_idx
#     cdef list arrs
#     cdef list uni, idx
#     cdef int search 
    
#     temp_idx = [iterat[0] for iterat in sorted(enumerate(arr), key=lambda x:x[1])] # sorted indexes  (pure python code)
#     arrs = [arr[i] for i in temp_idx] #sorted array 
    
#     #first element 
#     uni = np.array([arrs[0]])
#     idx = [temp_idx[0]]
    
#     for i in range(1,len(arrs)):
#         search = bisectionSearch(uni,arrs[i])
#         if search == -1:
#             uni.append(arrs[i])
#             idx.append(temp_idx[i])
#     return uni,idx
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double fmin(double[:] arr): # import note: You need to pass the memory view
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
cdef double fmax(double[:] arr): # import note: You need to pass the memory view
    # to a nogil function like this otherwise you get instability 
    cdef int n = arr.shape[0]
    cdef double tmp = arr[0]
    cdef int i 
    for i in range(1,n):
        if tmp<arr[i]:
            tmp = arr[i]
    return tmp

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double mean(double[:] arr):
    #compute mean average of an array 
    cdef int n = arr.shape[0]
    cdef double tmp = 0
    cdef double a 
    cdef int i 
    for i in range(n):
        tmp = tmp + arr[i]
    a = tmp/n
    return a 

@cython.boundscheck(False)
@cython.wraparound(False)
cdef ccw(p,q,r):#code expects points as p=(x,y) and so on ... 
    # Checks if points in a triangle are ordered counter clockwise. 
    val=((q[1]-p[1])*(r[0]-q[0]))-((q[0]-p[0])*(r[1]-q[1]))
    if val == 0:
        return 0 # lines are colinear
    elif val > 0:
        return 1 # points are oreintated clockwise
    elif val < 0:
        return 2 # points are oreintated counter clockwise 

@cython.boundscheck(False)
@cython.wraparound(False)    
cdef dlong tetrasignp(double[:] a, double[:] b, double [:] c,  double [:] d):
    #parrall nogil tetra sign code, but needs memory views passed to it 
    cdef double v00,v01,v02,v10,v11,v12,v20,v21,v22,s0,s1,s2, N 
    
    #matrix entries 
    v00 = b[0] - a[0] # x row
    v01 = c[0] - a[0]
    v02 = d[0] - a[0] 
    v10 = b[1] - a[1] # y row 
    v11 = c[1] - a[1] 
    v12 = d[1] - a[1]
    v20 = b[2] - a[2] # z row 
    v21 = c[2] - a[2]
    v22 = d[2] - a[2]
    
    #compute cross product
    s0 = v01*v12 - v02*v11
    s1 = v02*v10 - v00*v12
    s2 = v00*v11 - v01*v10
    
    #compute dot product (gives orientation)
    N = s0*v20 + s1*v21 + s2*v22
    
    if N>0:
        return 1 # points are clockwise
    elif N<0:
        return 2 # points are counter clockwise
    else:
        return 0 # something dogdey has happened as all points are on the same plane   

@cython.boundscheck(False)    
cdef void sortInt(dlong[:] arr, int n): 
    #adapted from https://www.geeksforgeeks.org/insertion-sort/
    #sorts array in place
    cdef Py_ssize_t i, j
    cdef dlong key
    for i in range(1, n): 
        key = arr[i] 
        j = i-1
        while j >= 0 and key < arr[j] : 
                arr[j + 1] = arr[j] 
                j -= 1
        arr[j + 1] = key 

cdef dlong mergeInts(int a, int b, int c, int pad): # merge 3 ints 
    cdef dlong x = np.floor(a*10**pad + b) # merge a and b
    return np.floor(x*10**pad + c) # merge c 

cdef dlong mergeInt(int a, int b, int pad): #merge 2 ints 
    return np.floor(a*10**pad + b) # merge a and b


# @cython.boundscheck(False)    
# @cython.wraparound(False)             
def neighSearchSafe(dlong[:,:] connection, int cell_type = 10):
    """Search for neighbours of each element within a mesh. 
    Algorithm updated in June 25 to be more robust against 
    integer overflow errors. 
    
    Parameters
    -------------
    connection: np.array 
        N by 4 array, describes how mesh nodes map to elements 
    cell_type: int
        Cell type indentifier (see notes)
    
    Returns
    --------------
    neigh: np.array 
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    
    Notes 
    -----
    cell types for each type of element supported by meshtools 
    
    """
    cdef Py_ssize_t i, j, ii, jj # variables used in loops 
    cdef int numel = connection.shape[0] # number of elements 
    cdef int npere = 4 # number of nodes per element (defaults for tetra)
    cdef int npfac = 3 # number of nodes per face 
    cdef int nnper = 4 # number of neighbours per element 
    cdef dlong[:,:] nmap 
    cdef dlong[:,:] nmap_tri = np.asarray([[0,1,2],[1,2,0]], dtype=npy_long)
    cdef dlong[:,:] nmap_quad = np.asarray([[0,1,2,3],[1,2,3,0]], dtype=npy_long)
    cdef dlong[:,:] nmap_tetra = np.asarray([[1,0,0,0],[2,3,1,1],[3,2,3,2]], dtype=npy_long)
    
    # define parameters 
    if cell_type==5: #then elements are triangles
        npere = 3 
        npfac = 2 
        nnper = 3 
        nmap = nmap_tri 
    elif cell_type==8 or cell_type==9: #elements are quads
        npere = 4 
        npfac = 2 
        nnper = 4
        nmap = nmap_quad 
    # elif cell_type == 11 or cell_type == 12: # elements are voxels (not tested)
    #     npere = 8 
    #     npfac = 4
    #     nnper = 6
    elif cell_type == 10:# elements are tetrahedra 
        npere = 4 # number of nodes per element (defaults for tetra)
        npfac = 3 # number of nodes per face 
        nnper = 4 # number of neighbours per element 
        nmap = nmap_tetra 
        
    else: 
        raise Exception('Unknown cell_type passed to mashCalc')
    
    #face arrays 
    cdef dlong[:] face = np.zeros(npfac, dtype=npy_long)

    # how does this work? Need to find the combination of nodes which are shared on each face 
    
    # create sorted arrays so we can expliot bisection search 
    cdef np.ndarray[dlong, ndim=2] face_connection_expanded = np.zeros((numel*npere,npfac), dtype=npy_long) # expanded table for storing faces 
    cdef np.ndarray[dlong, ndim=1] face_elementid_reference = np.zeros(numel*npere, dtype=npy_long) # reference for face elements 
    cdef np.ndarray[dlong, ndim=2] face_connection_expanded_sorted = np.zeros((numel*npere,npfac), dtype=npy_long) # sorted connection matrix 
    cdef np.ndarray[dlong, ndim=1] face_elementid_reference_sorted = np.zeros(numel*npere, dtype=npy_long) # 
    cdef np.ndarray[dlong, ndim=1] fconnection_isort = np.zeros(numel*npere, dtype=npy_long) # sorted indices of connection matrix 
    cdef np.ndarray[dlong, ndim=2] neigh = np.zeros((numel,nnper),dtype=npy_long)-1 # allocate space for neighbour matrix    
    cdef dlong iimax = numel*npere 
    cdef int npmatch = 0 # number of node matches per face (should be a max of 3)
    cdef int fcmatch = 0 # number of matching faces per node combination (should be a max of 2)
    cdef np.ndarray[dlong, ndim=1] matching_index = np.zeros(2,dtype=npy_long)
    # allocate memory views for efficient lookup 

    for i in range(numel): 
        for j in range(npere):
            #setup nodes which are part of face
            for ii in range(npfac): 
                face[ii] = connection[i,nmap[ii,j]]
            # cache face connection into expanded matrix -> will use this to lookup if the face is referenced by other elements 
            sortInt(face, npfac)
            for jj in range(npfac):
                face_connection_expanded[i+(j*numel), jj] = face[jj]
            face_elementid_reference[i+(j*numel)] = i 

    # sort by the first face node 
    fconnection_isort = np.argsort(face_connection_expanded[:,0])
    face_elementid_reference_sorted = face_elementid_reference[fconnection_isort]
    face_connection_expanded_sorted = face_connection_expanded[fconnection_isort,:]
 
    # loop through each element and find the matching combinations of nodes associated with each face 
    for i in range(numel): 
        for j in range(npere):
            #setup nodes which are part of face
            for ii in range(npfac): 
                face[ii] = connection[i,nmap[ii,j]]
            sortInt(face, npfac)

            ii = bisectionSearch(face_connection_expanded_sorted[:,0], face[0]) # o is the lookup index 
            # rewind until the the start of the part of the matrix/table with the relevant node is found 
            while face_connection_expanded_sorted[ii-1,0] == face[0]:
                ii -= 1 
                if ii == 0: # fail safe against end of array 
                    break 

            # there should be 2 exact matches in the faces which are found if the face has a neighbour 
            fcmatch = 0 
            for jj in range(2): 
                matching_index[jj] = -1 
            while face_connection_expanded_sorted[ii,0] == face[0] and ii < iimax:
                npmatch = 0 # set the number of matches to zero 
                for jj in range(npfac):
                    if face_connection_expanded_sorted[ii,jj] == face[jj]: 
                        npmatch += 1 
                if npmatch == npfac: # then we have found a matching face 
                    matching_index[fcmatch] = face_elementid_reference_sorted[ii]
                    fcmatch += 1 
                if fcmatch == 2: # pick which ever matching index is not the current element in question 
                    if matching_index[0] == i: 
                        neigh[i,j] = matching_index[1] 
                    else: 
                        neigh[i,j] = matching_index[0] 
                    break 
                ii += 1 
                if ii == iimax: 
                    break 
    return neigh 

@cython.boundscheck(False)    
@cython.wraparound(False)             
def neighSearch(dlong[:,:] connection, int cell_type = 10):
    """Search for neighbours of each element within a mesh. 
    Algorithm updated in June 25 to be more robust against 
    integer overflow errors. 
    
    Parameters
    -------------
    connection: np.array 
        N by 4 array, describes how mesh nodes map to elements 
    cell_type: int
        Cell type indentifier (see notes)
    
    Returns
    --------------
    neigh: np.array 
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    
    Notes 
    -----
    cell types for each type of element supported by meshtools 
    
    """
    cdef Py_ssize_t i, j, ii, jj # variables used in loops 
    cdef int numel = connection.shape[0] # number of elements 
    cdef int npere = 4 # number of nodes per element (defaults for tetra)
    cdef int npfac = 3 # number of nodes per face 
    cdef int nnper = 4 # number of neighbours per element 
    cdef dlong[:,:] nmap 
    cdef dlong[:,:] nmap_tri = np.asarray([[0,1,2],[1,2,0]], dtype=npy_long)
    cdef dlong[:,:] nmap_quad = np.asarray([[0,1,2,3],[1,2,3,0]], dtype=npy_long)
    cdef dlong[:,:] nmap_tetra = np.asarray([[1,0,0,0],[2,3,1,1],[3,2,3,2]], dtype=npy_long)
    
    # define parameters 
    if cell_type==5: #then elements are triangles
        npere = 3 
        npfac = 2 
        nnper = 3 
        nmap = nmap_tri 
    elif cell_type==8 or cell_type==9: #elements are quads
        npere = 4 
        npfac = 2 
        nnper = 4
        nmap = nmap_quad 
    # elif cell_type == 11 or cell_type == 12: # elements are voxels (not tested)
    #     npere = 8 
    #     npfac = 4
    #     nnper = 6
    elif cell_type == 10:# elements are tetrahedra 
        npere = 4 # number of nodes per element (defaults for tetra)
        npfac = 3 # number of nodes per face 
        nnper = 4 # number of neighbours per element 
        nmap = nmap_tetra 
        
    else: 
        raise Exception('Unknown cell_type passed to mashCalc')
    
    #face arrays 
    cdef dlong[:] face = np.zeros(npfac, dtype=npy_long)

    # how does this work? Need to find the combination of nodes which are shared on each face 
    
    # create sorted arrays so we can expliot bisection search 
    cdef np.ndarray[dlong, ndim=2] face_connection_expanded = np.zeros((numel*npere,npfac), dtype=npy_long) # expanded table for storing faces 
    cdef np.ndarray[dlong, ndim=1] face_elementid_reference = np.zeros(numel*npere, dtype=npy_long) # reference for face elements 
    cdef np.ndarray[dlong, ndim=2] face_connection_expanded_sorted = np.zeros((numel*npere,npfac), dtype=npy_long) # sorted connection matrix 
    cdef np.ndarray[dlong, ndim=1] face_elementid_reference_sorted = np.zeros(numel*npere, dtype=npy_long) # 
    cdef np.ndarray[dlong, ndim=1] fconnection_isort = np.zeros(numel*npere, dtype=npy_long) # sorted indices of connection matrix 
    cdef np.ndarray[dlong, ndim=2] neigh = np.zeros((numel,nnper),dtype=npy_long)-1 # allocate space for neighbour matrix    
    cdef dlong iimax = numel*npere 
    cdef int npmatch = 0 # number of node matches per face (should be a max of 3)
    cdef int fcmatch = 0 # number of matching faces per node combination (should be a max of 2)
    cdef dlong[:] matching_index = np.zeros(2,dtype=npy_long)
    # allocate memory views for efficient lookup 
    cdef dlong[:,:] face_connection_expandedv = face_connection_expanded 
    cdef dlong[:] face_elementid_referencev = face_elementid_reference
    cdef dlong[:,:] neighv = neigh 

    for i in range(numel): 
        for j in range(npere):
            #setup nodes which are part of face
            for ii in range(npfac): 
                face[ii] = connection[i,nmap[ii,j]]
            # cache face connection into expanded matrix -> will use this to lookup if the face is referenced by other elements 
            sortInt(face, npfac)
            for jj in range(npfac):
                face_connection_expandedv[i+(j*numel), jj] = face[jj]
            face_elementid_referencev[i+(j*numel)] = i 

    # sort by the first face node 
    fconnection_isort = np.argsort(face_connection_expanded[:,0])
    face_elementid_reference_sorted = face_elementid_reference[fconnection_isort]
    face_connection_expanded_sorted = face_connection_expanded[fconnection_isort,:]
    # for i in range(numel*npere): 
    #     for j in range(npfac): 
    #         face_connection_expanded_sorted[i,j] = face_connection_expanded[fconnection_isort[i],j]

    cdef dlong[:,:] face_connection_expanded_sortedv = face_connection_expanded_sorted
    cdef dlong[:] face_elementid_reference_sortedv = face_elementid_reference_sorted
 
    # loop through each element and find the matching combinations of nodes associated with each face 
    for i in range(numel): 
        for j in range(npere):
            #setup nodes which are part of face
            for ii in range(npfac): 
                face[ii] = connection[i,nmap[ii,j]]
            sortInt(face, npfac)

            ii = bisectionSearch(face_connection_expanded_sortedv[:,0], face[0]) # o is the lookup index 
            # rewind until the the start of the part of the matrix/table with the relevant node is found 
            while face_connection_expanded_sortedv[ii-1,0] == face[0]:
                ii -= 1 
                if ii == 0: # fail safe against end of array 
                    break 

            # there should be 2 exact matches in the faces which are found if the face has a neighbour 
            fcmatch = 0 
            for jj in range(2): 
                matching_index[jj] = -1 
            while face_connection_expanded_sortedv[ii,0] == face[0] and ii < iimax:
                npmatch = 0 # set the number of matches to zero 
                for jj in range(npfac):
                    if face_connection_expanded_sortedv[ii,jj] == face[jj]: 
                        npmatch += 1 
                if npmatch == npfac: # then we have found a matching face 
                    matching_index[fcmatch] = face_elementid_reference_sortedv[ii]
                    fcmatch += 1 
                if fcmatch == 2: # pick which ever matching index is not the current element in question 
                    if matching_index[0] == i: 
                        neighv[i,j] = matching_index[1] 
                    else: 
                        neighv[i,j] = matching_index[0] 
                    break 
                ii += 1 
                if ii == iimax: 
                    break 
    return neigh 
    
@cython.boundscheck(False)
@cython.wraparound(False)
def faces3d(dlong[:,:] connection, dlong[:,:] neigh):
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
    cdef Py_ssize_t i,j #loop variables
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    cdef int c = 0
    cdef np.ndarray[dlong, ndim=2] nmap = np.array([[1, 2, 3], [0, 3, 2], 
                                                   [0, 1, 3], [0, 1, 2]], dtype=npy_long)
    
    #first thing is to find all the cases where neigh == 1 
    cdef list idx = [] # this is the index where a outside edge is 
    cdef list fnum = [] # this is the face number ( 1 2 3 or 4)
    for i in range(numel):
        for j in range(npere):
            if neigh[i,j] == -1:
                idx.append(i)
                fnum.append(j)
                c+=1
    
    cdef dlong[:] idxa = np.array(idx,dtype=npy_long)
    cdef dlong[:] fnuma = np.array(fnum,dtype=npy_long)
    cdef np.ndarray[dlong, ndim=1] idxo = np.array(idx,dtype=npy_long)
    
    cdef int nfaces = len(idx)
    cdef int fidx, fnumi 
    
    cdef np.ndarray[dlong, ndim=2] fconnection = np.zeros((nfaces,3),dtype=npy_long) # face connection matrix 
    cdef dlong[:,:] fconnectionv = fconnection
    
    for i in range(nfaces):
        fidx = idxa[i]
        fnumi = fnuma[i]
        for j in range(3):
            fconnectionv[i,j] = connection[fidx,nmap[fnumi,j]]
            #find missing node in future? 
        
    return fconnection, idxo 
    
@cython.boundscheck(False)    
@cython.wraparound(False)             
def neighPrism(dlong[:,:] connection):
    """Compute neighbours of each element within a prism mesh 
    
    Parameters
    -------------
    connection: np.array 
        N by 4 array, describes how mesh nodes map to elements 
    
    Returns
    --------------
    neigh: np.array 
        Corresponding neighbour indexes for each face of the cells in the 
        connection matrix
    """
    
    cdef Py_ssize_t i, j
    cdef int numel = connection.shape[0]
    cdef int npere = 5
    #face arrays 
    cdef dlong[:] faces = np.zeros(4, dtype=npy_long) # face on side 
    cdef dlong[:] facet = np.zeros(3, dtype=npy_long) # face on top / bottom 
    cdef dlong[:] a = np.asarray([0,3], dtype=npy_long)
    cdef dlong[:] b = np.asarray([1,4], dtype=npy_long)  
    cdef dlong[:] c = np.asarray([2,5], dtype=npy_long)  
    cdef dlong[:] d = np.asarray([0,1,2] , dtype=npy_long)
    cdef dlong[:] e = np.asarray([1,2,0] , dtype=npy_long)
    cdef dlong[:] f = np.asarray([3,4,3] , dtype=npy_long)
    cdef dlong[:] g = np.asarray([4,5,5] , dtype=npy_long)

    #combination array 
    cdef np.ndarray[dlong, ndim=2] tri_combo = np.zeros((numel,5),dtype=npy_long ,order='C') # allocate space for tri_combo 
    cdef dlong[:,:] tri_combov = tri_combo
    #maximum number of node digits 
    cdef int num_node = max([max(connection[:,0]),max(connection[:,1]),
                             max(connection[:,2]),max(connection[:,3])])
    cdef int pad = len(str(num_node))
    
    #assign unique node combinations to each element edge, repeats pick out neighbours
    for i in range(numel):
        #5 faces per prism 
        for j in range(2):
            #setup nodes which are part of face
            facet[0] = connection[i,a[j]]
            facet[1] = connection[i,b[j]] 
            facet[2] = connection[i,c[j]]
            # sort face indexes 
            sortInt(facet,3)
            # assign each face an organised and unique code
            tri_combov[i,j] = mergeInts(facet[0],facet[1],facet[2],pad)
        for j in range(3):
            faces[0] = connection[i,d[j]]
            faces[1] = connection[i,e[j]] 
            faces[2] = connection[i,f[j]]
            faces[3] = connection[i,g[j]]
            sortInt(faces,4)
            # assign each face an organised and unique code
            tri_combov[i,j+2] = mergeInt(mergeInt(faces[0],faces[1],pad),
                                         mergeInt(faces[2],faces[3],pad),
                                         pad)
        

    cdef np.ndarray[dlong, ndim=2] neigh = np.zeros((numel,5),dtype=npy_long) # allocate space for neighbour matrix        
    cdef dlong[:,:] neighv = neigh  
    
    #using binary search and sorted lists for efficient index lookup 
    cdef np.ndarray[dlong, ndim=1] tri_flatten = tri_combo.T.flatten() #all faces in one array 
	
    cdef np.ndarray[dlong, ndim=1] temp_idx = np.argsort(tri_flatten).astype(npy_long,order='C') # sorted indexes  (pure python code)
    cdef np.ndarray[dlong, ndim=1] tri_sort = tri_flatten[temp_idx] # sorted faces (forward direction)
    cdef dlong[:] tri_sortv = tri_sort
    cdef int maxo = len(tri_sort)-1

    cdef dlong idx
    cdef Py_ssize_t o 
    #create lookup array 
    cdef np.ndarray[dlong, ndim=1] lookup = np.zeros((numel*5),dtype=npy_long)
    for i in range(numel):
        lookup[i] = i
        lookup[i+numel] = i
        lookup[i+(2*numel)] = i
        lookup[i+(3*numel)] = i
        lookup[i+(4*numel)] = i

    cdef np.ndarray[dlong, ndim=1] lookup_idx = np.asarray(lookup[temp_idx],order='C')
    cdef dlong[:] lookup_idxv = lookup_idx
  
    #loop through each element 
    for i in range(numel):
        for j in range(npere):
            o = bisectionSearch(tri_sortv,tri_combov[i,j]) # only works on a sorted array
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
    
    return neigh
    
@cython.boundscheck(False)
@cython.wraparound(False)
def facesPrism(dlong[:,:] connection, double[:,:] node, dlong[:,:] neigh):
    """Return external faces of a 3D prism mesh 
    
    Parameters
    ----------
    connection: np.array (int)
        N by 4 array, describes how mesh nodes map to elements 

    neigh: np.array (int)
        N by 4 array, describes indices of neighbour elements for each cell 
        (output of neigh3d)

    Returns
    -------
    fconnection : list
        N by 3 matrix with the with the nodes of external faces 
    idxa: np.array
        N by 1 array. Indices of original mesh elements 

    """
    cdef Py_ssize_t i,j,k #loop variables
    cdef int numel = connection.shape[0]
    cdef int c = 0
    cdef list nmap_list = [
        [0, 1, 2, -1], 
        [3, 4, 5, -1], 
        [0, 1, 4, 3], 
        [1, 2, 5, 4],
        [0, 2, 5, 3]]
    cdef np.ndarray[dlong, ndim=2] nmap = np.array(nmap_list, dtype = npy_long)
    
    #first thing is to find all the cases where neigh == 1 
    cdef list idx = [] # this is the index where a outside edge is 
    cdef list fnum = [] # this is the face number ( 1 2 3 4 or 5)
    for i in range(numel):
        for j in range(5):
            if neigh[i,j] == -1:
                idx.append(i)
                fnum.append(j)
                c+=1
    
    cdef dlong[:] idxa = np.array(idx,dtype=npy_long)
    cdef np.ndarray[dlong, ndim=1] idxo = np.array(idx,dtype=npy_long)
    cdef dlong[:] fnuma = np.array(fnum,dtype=npy_long)
    
    cdef int nfaces = len(idx)
    cdef int fidx, fnumi, search 
    
    cdef list fcoords = [[]]*nfaces 
    cdef tuple vert 
    cdef list vertices = []
    
    for i in range(nfaces):
        fidx = idxa[i]
        fnumi = fnuma[i]
        vertices = []
        if fnumi < 2:
            k=3
        else:
            k=4
        for j in range(k):
            search = connection[fidx, nmap[fnumi,j]]
            vert = (node[search,0], node[search,1], node[search,2])
            vertices.append(vert)
        fcoords[i] = vertices 
        
    return fcoords, idxo
 
@cython.boundscheck(False)
@cython.wraparound(False)
def sortNeigh(dlong[:,:] neigh, dlong[:] zone):
    """Sort neighbour matrix for input into R3t. 
    -----------
    neigh: nd array 
        M by N describing element neighbours 
    zone: nd array 
        N rows which corresponds to the zone number of elements. 
        
    Returns
    -----------
    new_connection: np array 
    neigh: nd array
        Prepared neighbour matrix 
    """
    cdef Py_ssize_t i,j 
    cdef int numel = neigh.shape[0]
    cdef int npere = neigh.shape[1]
    cdef np.ndarray[dlong, ndim=2] new_neigh = np.zeros((numel, npere), dtype=npy_long)
    cdef dlong[:,:] new_neighv= new_neigh; 
    
    #reorder each row (acsending but outside elements at end of row) 
    cdef dlong like_inf = numel*2 #make a big number 
    for i in range(numel):
        for j in range(npere):
            new_neighv[i,j] = neigh[i,j]
            if neigh[i,j] == -1: #check if outside element 
                new_neighv[i,j] = like_inf # if outside assign big number 
            elif zone[i] != zone[neigh[i,j]]:# check if neighbour in different zone 
                new_neighv[i,j] = like_inf 
                
        sortInt(new_neigh[i,:],npere) # sort in that part of the row
        for j in range(npere):
            if new_neighv[i,j] == like_inf: # replace outside values with -1 
                new_neighv[i,j] = -1
    return new_neigh
                
    
@cython.boundscheck(False)
@cython.wraparound(False)
def splitTri(dlong[:,:] connection, double[:,:] node):
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
    cdef dlong[:] n = np.zeros(6,dtype=npy_long) # new node numbers 
    cdef float mx, my, mz
    cdef dlong[:] nodes = np.zeros(2,dtype=npy_long)
    cdef dlong mn, nno 
    cdef Py_ssize_t i, j, tmpi, search

    #node handling 
    cdef int num_nodes = len(node_x) # number of nodes 
    cdef int numel = connection.shape[0] # number of elements 
    
    #matrices 
    cdef np.ndarray[dlong, ndim=2] new_connection = np.zeros((numel*4,3),dtype = npy_long) # new connection matrix 
    cdef dlong[:,:] new_connectionv = new_connection # memory view 

    cdef np.ndarray[double, ndim=2] new_node = np.zeros((numel*3,3), dtype=float)
    cdef np.ndarray[dlong, ndim=1] new_node_idx = np.zeros(numel*3, dtype=npy_long)
    cdef np.ndarray[dlong, ndim=1] og_el_id = np.zeros(numel*3,dtype=npy_long)
    cdef double[:,:] new_nodev = new_node
    cdef dlong[:] og_el_idv = og_el_id 
    cdef int pad = len(str(num_nodes))+1
    
    cdef dlong[:] a = np.array([0,1,2],dtype=npy_long)
    cdef dlong[:] b = np.array([1,2,0],dtype=npy_long)
    
    #remap the element nodes with this matrix 
    cdef list remap = [[0, 5, 3],
                       [1, 4, 3],
                       [2, 5, 4],
                       [3, 4, 5]]
    cdef dlong[:,:] remapv = np.asarray(remap,dtype=npy_long)
    
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
    cdef np.ndarray[dlong, ndim=1] idxa, unicl
    unicl, idxa = np.unique(new_node_idx, return_index=True) # unique and ordered node configs 
    cdef dlong[:] node_id = np.arange(len(unicl), dtype=npy_long) + num_nodes 
    
    ### map back to elements #### 
    for i in range(numel):
        for j in range(3):
            nno = connection[i,j] # node number
            n[j] = nno 
        
        for j in range(3):
            nodes[0] = n[a[j]]; nodes[1] = n[b[j]]
            sortInt(nodes,2)
            mn = mergeInt(nodes[0],nodes[1],pad)
            search = bisectionSearch(unicl,mn)
            
            n[j+3] = node_id[search]#reference index
    
        tmpi = i*4 # temporary index for indexing new connection matrix
        for j in range(4):
            new_connectionv[tmpi+j,0] = n[remapv[j,0]]
            new_connectionv[tmpi+j,1] = n[remapv[j,1]]
            new_connectionv[tmpi+j,2] = n[remapv[j,2]]
            
    ### make new node matrix ### 
    cdef int added_nodes = len(idxa)
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
def splitTetra(dlong[:,:] connection, double[:,:] node):
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
    #define vars 
    cdef double[:] node_x = np.asarray(node[:,0], dtype=float)
    cdef double[:] node_y = np.asarray(node[:,1], dtype=float)
    cdef double[:] node_z = np.asarray(node[:,2], dtype=float)
    cdef double[:] x = np.zeros(4,dtype=float)
    cdef double[:] y = np.zeros(4,dtype=float)
    cdef double[:] z = np.zeros(4,dtype=float)
    cdef dlong[:] n = np.zeros(10,dtype=npy_long) # new node numbers 
    cdef double mx, my, mz
    cdef dlong[:] nodes = np.zeros(2,dtype=npy_long)
    cdef dlong mn, nno 
    cdef Py_ssize_t i, j, tmpi, search #loop variables 

    cdef int num_nodes = len(node_x)
    cdef int numel = connection.shape[0]

    cdef np.ndarray[dlong, ndim=2] new_connection = np.zeros((numel*8,4),dtype=npy_long) # new connection matrix 
    cdef dlong[:,:] new_connectionv = new_connection

    cdef np.ndarray[double, ndim=2] new_node = np.zeros((numel*6,3),dtype=float) # used to reference already computed node centres
    cdef np.ndarray[dlong, ndim=1] new_node_idx = np.zeros(numel*6,dtype=npy_long)
   
    cdef np.ndarray[dlong, ndim=1] og_el_id = np.zeros(numel*6,dtype=npy_long)
    cdef double[:,:] new_nodev = new_node
    cdef dlong[:] og_el_idv = og_el_id 
    
    cdef int pad = len(str(num_nodes))
    
    cdef dlong[:] a = np.array([0,1,2,3,3,3],dtype=npy_long)
    cdef dlong[:] b = np.array([1,2,0,0,1,2],dtype=npy_long)
    
    #remap the element nodes with this matrix 
    cdef list remap = [[4, 7, 6, 0],
                       [5, 9, 6, 2],
                       [8, 9, 7, 3],
                       [8, 9, 7, 6],
                       [8, 4, 7, 6],
                       [8, 5, 9, 6],
                       [4, 5, 6, 8],
                       [4, 1, 5, 8]]
    cdef dlong[:,:] remapv = np.asarray(remap,dtype=npy_long)
    
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
    cdef np.ndarray[dlong, ndim=1] idxa, unicl
    unicl, idxa = np.unique(new_node_idx, return_index=True) # unique and ordered node configs 
    cdef dlong[:] node_id = np.arange(len(unicl), dtype=npy_long) + num_nodes
    
    ### map back to elements #### 
    for i in range(numel):
        for j in range(4):
            nno = connection[i,j] # node number
            n[j] = nno 
        
        for j in range(6):
            nodes[0] = n[a[j]]; nodes[1] = n[b[j]]
            sortInt(nodes,2)
            mn = mergeInt(nodes[0],nodes[1],pad)
            search = bisectionSearch(unicl,mn)
            n[j+4] = node_id[search]#reference index
    
        tmpi = i*8 # temporary index for indexing new connection matrix
        for j in range(8):
            new_connectionv[tmpi+j,0] = n[remapv[j,0]]
            new_connectionv[tmpi+j,1] = n[remapv[j,1]]
            new_connectionv[tmpi+j,2] = n[remapv[j,2]]
            new_connectionv[tmpi+j,3] = n[remapv[j,3]]
            
    ### make new node matrix ### 
    cdef int added_nodes = len(idxa) # number of new nodes added to mesh 
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
def orderTetra(dlong[:,:] connection, double[:,:] node):
    """ Organise tetrahedral element nodes into a clockwise order
    
    following solution at: 
    https://stackoverflow.com/questions/10612829/tetrahedron-orientation-for-triangle-meshes
    
    Parameters
    ----------
    connection: np array (dlong)
        Mesh connection matrix, N by 4 array. 
    node: np array (float64)
        N by 3 array describing the mesh nodes 
    Returns
    -------
    con : np array (dlong)
        Mesh connection matrix which has been 'corrected' such that all elements
        are counterclockwise. 
    count : int
        Number of elements with switched nodes (ie been 'corrected')
    ccw : np array(dlong)
        clockwise orientation vector
    """
    #get the number of nodes and elements etc 
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    
    cdef np.ndarray[dlong,ndim=2] con = np.zeros((numel,npere), dtype=npy_long)# new connection matrix 
    cdef dlong[:,:] conv = con # connection memeory view
    
    cdef np.ndarray[double, ndim=1] node_x = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef np.ndarray[double, ndim=1] node_y = np.asarray(node[:,1],dtype=float)
    cdef np.ndarray[double, ndim=1] node_z = np.asarray(node[:,2],dtype=float)
    cdef double[:] nodex = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef double[:] nodey = np.asarray(node[:,1],dtype=float)
    cdef double[:] nodez = np.asarray(node[:,2],dtype=float)

    #looping variables 
    cdef double v00,v01,v02,v10,v11,v12,v20,v21,v22,s0,s1,s2, N 
    
    cdef np.ndarray[dlong, ndim=1] ccw = np.zeros(numel,dtype=npy_long) # clockwise array 
    cdef dlong[:] ccwv = ccw #clockwise view
    cdef Py_ssize_t i 
    cdef int ei # loop integers 
    cdef np.ndarray[dlong, ndim=1] count = np.zeros(numel,dtype=npy_long)
    cdef dlong[:] countv = count
    cdef int count_out #rolling total for the number of corrected elements 

    for i in range(numel):
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
def orderQuad(dlong[:,:] connection, double[:,:] node):
    """Order quaderlateral elements counterclockwise 
    
    Parameters
    ----------
    connection: np array (dlong)
        Mesh connection matrix, N by 4 array. 
    node: np array (float64)
        N by 3 array describing the mesh nodes 

    Returns
    -------
    con : np array (dlong)
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
    
    cdef np.ndarray[dlong, ndim=2] con = np.zeros((numel,npere), dtype=npy_long)# new connection matrix 
    cdef int[:,:] conv = con # connection view
    
    cdef np.ndarray[double, ndim=1] node_x = np.asarray(node[:,0],dtype=float) 
    cdef np.ndarray[double, ndim=1] node_y = np.asarray(node[:,1],dtype=float)
    cdef np.ndarray[double, ndim=1] node_z = np.asarray(node[:,2],dtype=float)
    
    #loop varaibles 
    cdef double[:] x = np.zeros(4,dtype=float) # x array 
    cdef double[:] z = np.zeros(4,dtype=float) # z array 
    cdef double[:] xtmp = np.zeros(4,dtype=float)
    cdef double[:] ztmp = np.zeros(4,dtype=float)

    cdef np.ndarray[dlong, ndim=1] count = np.zeros(numel,dtype=npy_long)
    cdef dlong[:] countv = count
    cdef int count_out 
    
    cdef double[:] theta = np.zeros(4,dtype=float) # theta array
    
    cdef dlong[:] order = np.zeros(4,dtype=npy_long) #ordering array
    
    cdef double pi = np.pi 
    cdef Py_ssize_t i, j, k
    cdef int eflag = 0 # error flag
    cdef int c, ei # count, error number 
    cdef double minx, minz
    cdef double xinf = 10*max(node_x) # infinite like x 
    cdef double zinf = 10*max(node_z) # infinite like z 
    cdef double tinf = 314  #infinite like theta 
    cdef double mtheta
    cdef double dx, dy, a

    #aim is to work out angles from point with minimum x and z coordinate 
    #min angle >>> bottom left most point 
    #max angle >>> upper left most point
    
    for i in range(numel):
        #could make loop parallel in future? 
        for j in range(4):
            x[j] = node_x[connection[i,j]]
            z[j] = node_z[connection[i,j]]
        
        #find starting x coordinate 
        minx = fmin(x)
        c = 0 # rolling count 
        #protect against colinear points on left hand side of quad 
        for j in range(4):
            if x[j] == minx:
                ztmp[j] = z[j]
                xtmp[j] = x[j] # put in array where x == min(x)
                c = c + 1
            else:
                ztmp[j] = zinf
                xtmp[j] = xinf
        
        #min z coordinate 
        minz = fmin(ztmp)        
        if c>2: # then there is more than one min x coordinate (presumably 2 at most)
            eflag = 1 
            ei = i 
        
        #create order angle 
        for j in range(4):
            if x[j] == minx and z[j] == minz:#baseline point 
                theta[j] = 0
            elif x[j] == minx and z[j] != minz: # colinear point 
                theta[j] = pi
            elif x[j] != minx and z[j] == minz: #colinear on z axis 
                theta[j] = pi/2
            else:
                dx = x[j] - minx
                dz = z[j] - minz
                a = atan(dz/dx)
                theta[j] = (pi/2) + a
       
        mtheta = fmin(theta) # min theta 
        for j in range(4):
            for k in range(4):
                if theta[k] == mtheta:
                    theta[k] = tinf
                    order[j] = k
                    mtheta = fmin(theta)
                    break  
        
        for j in range(4):# flag if order changes >>> count as ordered element 
            if order[j] != j:#order has been changed 
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
def surfaceCall(dlong[:,:] fconnection, double[:,:] node, double[:,:] cellcentres):
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
    
    cdef np.ndarray[dlong,ndim=1] ocheck = np.zeros(nfaces,dtype=npy_long)
    cdef dlong[:] ocheckv = ocheck

    #looping variables 
    cdef Py_ssize_t i,j
    cdef double bqz, xx, yy, xm, ym
    
    cdef double[:,:] q0 = np.zeros((nfaces,3))#query points 
    cdef double[:,:] q1 = np.zeros((nfaces,3))
    cdef double[:,:] p0 = np.zeros((nfaces,3))
    cdef double[:,:] p1 = np.zeros((nfaces,3))
    cdef double[:,:] p2 = np.zeros((nfaces,3))
    
    cdef dlong s1,s2,s3,s4,s5
    
    #construct query arrays for each element and extract surface cells 
    for i in range(nfaces):
        xx = 0
        yy = 0 
        for j in range(3):
            xx = xx + node_xv[fconnection[i,j]]
            yy = yy + node_yv[fconnection[i,j]]
            
        xm = xx/3 # x mid 
        ym = yy/3 # y mid 
        bqz = cellcentres[i,2] # zmid 
        
        q0[i,0] = xm
        q0[i,1] = ym
        q0[i,2] = tqz # top query point 
        
        q1[i,0] = xm
        q1[i,1] = ym
        q1[i,2] = bqz #bottom query point 
        
        p0[i,0] = node_xv[fconnection[i,0]] # corner 1 
        p0[i,1] = node_yv[fconnection[i,0]]
        p0[i,2] = node_zv[fconnection[i,0]]
        
        p1[i,0] = node_xv[fconnection[i,1]] # corner 2
        p1[i,1] = node_yv[fconnection[i,1]]
        p1[i,2] = node_zv[fconnection[i,1]]
        
        p2[i,0] = node_xv[fconnection[i,2]] # corner 3
        p2[i,1] = node_yv[fconnection[i,2]]
        p2[i,2] = node_zv[fconnection[i,2]]
        
        s1 = tetrasignp(q0[i,:],p0[i,:],p1[i,:],p2[i,:])
        s2 = tetrasignp(q1[i,:],p0[i,:],p1[i,:],p2[i,:])
        
        if s1 != s2: # then query point is either side of the triangle so probe
            s3 = tetrasignp(q0[i,:],q1[i,:],p0[i,:],p1[i,:])
            s4 = tetrasignp(q0[i,:],q1[i,:],p1[i,:],p2[i,:])
            s5 = tetrasignp(q0[i,:],q1[i,:],p2[i,:],p0[i,:])

            if s3 == s4 and s4 == s5:
                ocheckv[i] = 1 
                
    return ocheck

#nsizeA and finite element conductance calculation
@cython.boundscheck(False)
@cython.wraparound(False)
def conductanceCall(dlong[:,:] connection, int numnp, int typ=0):
    """Calculate the array size needed for the finite element conductance matrix
    in R2 class codes 

    Parameters
    ----------
    connection: numpy array 
        mesh connection matrix.
    numnp: int
        Number of nodes
    typ: vtk cell type
        The default is 0. Will raise an error if takes an unexpected value. 

    Returns
    -------
    nsizeA : int
        Number of connected nodes + number of nodes
    Nconnec: nd array
        Node connection matrix. 

    Note
    -----
    typ is 5 for triangle meshes, 8 or 9 for quads, 10 for tetrahedra, 13 for 
    a prism. 
    """
    cdef int nmax = 60
    cdef int numel = connection.shape[0]
    cdef int nedges 
    cdef int pad = len(str(numnp))
    cdef dlong[:] a
    cdef dlong[:] b
    #allocate lookup arrays a and b for different cell types 
    cdef a5 = np.asarray([0,1,2], dtype=npy_long)
    cdef b5 = np.asarray([1,2,0], dtype=npy_long)
    cdef a8 = np.asarray([0,1,2,3,0,1], dtype=npy_long)
    cdef b8 = np.asarray([1,2,3,0,2,3], dtype=npy_long)
    cdef a10 = np.asarray([0,1,2,3,3,3], dtype=npy_long)
    cdef b10 = np.asarray([1,2,0,0,1,2], dtype=npy_long)
    cdef a13 = np.asarray([0,1,2, 3,4,5, 0,1,2, 0,0, 1,1, 2,2], dtype=npy_long)
    cdef b13 = np.asarray([1,2,0, 4,5,3, 3,4,5, 5,4, 3,5, 3,4], dtype=npy_long)
    cdef np.ndarray[dlong, ndim=1] idx 
    cdef np.ndarray[dlong, ndim=1] counts # unique counts and indices
    cdef np.ndarray[dlong, ndim=1] uni
    cdef np.ndarray[dlong, ndim=1] combof # uni segment combinations 
    cdef np.ndarray[dlong, ndim=2] Nconnec = np.zeros((numnp,nmax),dtype=npy_long) - 1
    cdef dlong[:,:] Nconnecv = Nconnec

    #determine number of edges 
    if typ==5:#then elements are triangles
        nedges = 3
        a = a5 
        b = b5 
    elif typ==8 or typ==9:#elements are quads
        nedges = 6
        a = a8 
        b = b8
    elif typ == 10:# elements are tetrahedra 
        nedges = 6
        a = a10
        b = b10 
    elif typ == 13: # elements are 3d wedges 
        a = a13
        b = b13
        nedges = len(a)
    else: 
        raise ValueError('Cell type argument does not match vtk cell types,', 
                         'used with meshTools must be one of the following:',
                         '5, 8, 9, 10 or 13')

    cdef dlong[:,:] combo = np.zeros((numel,nedges),dtype=npy_long)
    cdef Py_ssize_t i,j,k #looping variables  
    cdef dlong merged 
    cdef dlong na, nb, nid
    
    #find unique node combinations 
    for i in range(numel):
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
                nid = Nconnecv[na,k]
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
    combof = np.asarray(combo,dtype=npy_long).flatten()
    uni, idx, counts = np.unique(combof, 
                                 return_index=True, 
                                 return_counts=True)
    
    #reorder Nconnec so that each row is ordered 
    for i in range(numnp):
        for j in range(nmax):
            if Nconnecv[i,j] == -1: #find up to point that actually values
                k = j
                break # break once -1 is found 
        sortInt(Nconnec[i,:],k) # sort in that part of the row
            
    if max(counts) > 60: 
        raise Exception('Too many connected nodes (>%i) check that mesh is not poorly formed.'%nmax)
        
    cdef int nsizeA = len(uni) + numnp # nsizeA can now be computed 
    
    return nsizeA, Nconnec


@cython.boundscheck(False)
@cython.wraparound(False)
def externalN(dlong[:,:] connection, double[:,:] node, dlong[:,:] neigh):
    """Get the external nodes of a triangle or tetrahedral mesh. Future plan 
    is to add a flag to determine if nodes are on top of the mesh or on the side. 

    Parameters
    ----------
    connection: np array 
        4 by N array of int describing indices of element vertices (nodes)
    node: np array 
        3 by N array of node coordinates 
    neigh: np.array (int)
        N by 4 array, describes indices of neighbour elements for each cell 
        (output of neigh3d)
        
    Returns
    -------
    enodes: np array 
        1 by N array of indices where nodes are on the edge of the mesh. 

    """
    # generic variables for function 
    cdef int flag3d = 0 
    cdef int i, j, k 
    cdef int numel = connection.shape[0]
    cdef int numnp = node.shape[0]
    cdef int nn = neigh.shape[1] # number of possible neighbours 
    cdef int nvert = connection.shape[1] # number of element vertices 
    cdef double xa, ya, za, xm, ym, zm, xmf, ymf, zmf 
    
    if connection.shape[1] == 3:
        flag3d = 0 
    elif connection.shape[1] == 4: 
        flag3d = 1 
    else:
        raise Exception('Unsupported mesh type for computing external nodes')
        
    cdef dlong[:] a, b, c 
    cdef dlong[:] a2d = np.asarray([0,1,2], dtype=npy_long)
    cdef dlong[:] b2d = np.asarray([1,2,0], dtype=npy_long)
    cdef dlong[:] a3d = np.asarray([1,0,0,0], dtype=npy_long) 
    cdef dlong[:] b3d = np.asarray([2,3,1,1], dtype=npy_long)  
    cdef dlong[:] c3d = np.asarray([3,2,3,2], dtype=npy_long)  
    if flag3d == 0: 
        a = a2d 
        b = b2d 
    else:
        a = a3d
        b = b3d
        c = c3d 
      
    cdef list enodesl = [] # will be used to store all found face nodes 
    cdef list surfaceflagl = [] # will be used to store if node is at the surface 
    cdef int enode0, enode1, enode2 # external nodes 
    
    cdef double[:] node_xv = np.asarray(node[:,0],dtype=float) # extract 1D arrays of node coordinates  
    cdef double[:] node_yv = np.asarray(node[:,1],dtype=float)
    cdef double[:] node_zv = np.asarray(node[:,2],dtype=float)
    # cdef np.ndarray[double,ndim=1] nodez =  np.asarray(node[:,2],dtype=float)
    cdef double tqz = (max(node_zv) - min(node_zv)) + max(node_zv) # query point in z axis
    
    cdef double[:] q0 = np.zeros(3)#query points 
    cdef double[:] q1 = np.zeros(3)
    cdef double[:] p0 = np.zeros(3)
    cdef double[:] p1 = np.zeros(3)
    cdef double[:] p2 = np.zeros(3)
    
    cdef dlong s1,s2,s3,s4,s5
    
    for i in range(numel):
        # skip if not an edge element 
        if min(neigh[i,:]) > -1: 
            continue 
        
        # compute cell centre 
        xa = 0
        ya = 0 
        za = 0         
        for j in range(nvert):
            xa += node_xv[connection[i,j]]
            ya += node_yv[connection[i,j]]
            za += node_zv[connection[i,j]]
        
        xm = xa/nvert 
        ym = ya/nvert 
        zm = za/nvert 
                          
        for j in range(nn):
            if neigh[i,j] == -1:
                enode0 = connection[i][a[j]]
                enode1 = connection[i][b[j]]
                enodesl.append(enode0)
                enodesl.append(enode1)
                if flag3d == 1: 
                    enode2 = connection[i][c[j]]
                    enodesl.append(enode2)
                    
                    p0[0] = node_xv[enode0]# corner 1 
                    p0[1] = node_yv[enode0]
                    p0[2] = node_zv[enode0]
                    
                    p1[0] = node_xv[enode1] # corner 2
                    p1[1] = node_yv[enode1]
                    p1[2] = node_zv[enode1]
                    
                    p2[0] = node_xv[enode2] # corner 3
                    p2[1] = node_yv[enode2]
                    p2[2] = node_zv[enode2]
                    
                    # work out mid point of external face 
                    xmf = (p0[0] + p1[0] + p2[0])/3 
                    ymf = (p0[1] + p1[1] + p2[1])/3 
                    
                    q0[0] = xmf
                    q0[1] = ymf
                    q0[2] = tqz # top query point 
                    
                    q1[0] = xmf
                    q1[1] = ymf
                    q1[2] = zm #bottom query point 
                    
                    s1 = tetrasignp(q0,p0,p1,p2)
                    s2 = tetrasignp(q1,p0,p1,p2)
                    
                    if s1 != s2: # then query point is either side of the triangle so probe
                        s3 = tetrasignp(q0,q1,p0,p1)
                        s4 = tetrasignp(q0,q1,p1,p2)
                        s5 = tetrasignp(q0,q1,p2,p0)
            
                        if s3 == s4 and s4 == s5:
                            for k in range(3):
                                surfaceflagl.append(1)
                        else:
                            for k in range(3):
                                surfaceflagl.append(0)
                    else:
                        for k in range(3):
                            surfaceflagl.append(0)
                else:
                    #now to compute if face boundary condition              
                    dx = abs(node_xv[enode0] - node_xv[enode1]) # work out approx edge dimensions 
                    dz = abs(node_zv[enode0] - node_zv[enode1]) # work out approx edge dimensions 
                    
                    if dx < 1e-16: # element on side of mesh 
                        surfaceflagl.append(0) 
                        surfaceflagl.append(0) 
                    else: 
                        p0[0] = node_xv[enode0]
                        p0[1] = node_zv[enode0]
                        p1[0] = node_xv[enode1]
                        p1[1] = node_zv[enode1]
                        p2[0] = xm
                        p2[1] = zm+dx+dz  
                        o = ccw(p0,p1,p2)  

                        if o == 1:
                            surfaceflagl.append(1) 
                            surfaceflagl.append(1) 
                        else:
                            surfaceflagl.append(0) 
                            surfaceflagl.append(0) 
                            
                            
    cdef np.ndarray[dlong, ndim=1] enodeu, enodeidx 
    enodeu, enodeidx = np.unique(enodesl,return_index=True) # get unique nodes 
    
    cdef np.ndarray[dlong, ndim=1] enodes = np.zeros(len(enodeu),dtype=npy_long)
    cdef np.ndarray[dlong, ndim=1] surfaceflag = np.zeros(len(enodeu),dtype=npy_long)

    for i in range(len(enodeu)):
        enodes[i] = enodesl[enodeidx[i]]
        surfaceflag[i] = surfaceflagl[enodeidx[i]]
    
    return enodes, surfaceflag 
    
    
@cython.boundscheck(False)
@cython.wraparound(False)
def fcrosscheck(dlong[:,:] fconnection1, dlong[:,:] fconnection2):
    """Cross check for overlapping faces in 2 different 2D face meshes... 
    Finds repeats of the first face connection matrix in the second. 
    Note, code may be vunrable to integer overflow for extremely large meshes. 
    
    Parameters
    -------------
    fconnection1: np.array 
        N by 3 array, describes how mesh nodes map to face elements 
    fconnection2: np.array 
        N by 3 array, describes how mesh nodes map to face elements 
    
    Returns
    --------------
    repeats: np.array 
        Boolian array, where index == 1 then the face of fonnection1 is repeated 
        in fconnection2. 
    """
    
    cdef Py_ssize_t i
    cdef int numel1 = fconnection1.shape[0]
    cdef int numel2 = fconnection2.shape[0]
    
    #face arrays 
    cdef dlong[:] face = np.zeros(3,dtype=npy_long)

    #combination arrays 
    cdef np.ndarray[dlong, ndim=1] combo1 = np.zeros((numel1,),dtype=npy_long ,order='C') # allocate space for combination matrix 
    cdef dlong[:] combov1 = combo1
    cdef np.ndarray[dlong, ndim=1] combo2 = np.zeros((numel2,),dtype=npy_long ,order='C') # allocate space for combination matrix 
    cdef dlong[:] combov2 = combo2
    #maximum number of node digits 
    cdef int num_node = max([max(fconnection1[:,0]),max(fconnection1[:,1]),
                             max(fconnection1[:,2]),max(fconnection2[:,0]),
                             max(fconnection2[:,1]),max(fconnection2[:,2])])
    cdef int pad = len(str(num_node))
    
    #assign unique node combinations to each face 
    for i in range(numel1):
        #setup nodes which are part of face
        face[0] = fconnection1[i,0]
        face[1] = fconnection1[i,1] 
        face[2] = fconnection1[i,2]
        # sort face indexes 
        sortInt(face,3)
        # assign each face an organised and unique code
        combov1[i] = mergeInts(face[0],face[1],face[2],pad)
    for i in range(numel2):
        #setup nodes which are part of face
        face[0] = fconnection2[i,0]
        face[1] = fconnection2[i,1] 
        face[2] = fconnection2[i,2]
        # sort face indexes 
        sortInt(face,3)
        # assign each face an organised and unique code
        combov2[i] = mergeInts(face[0],face[1],face[2],pad)
        
        
    #using binary search and sorted lists for efficient lookup 
    cdef np.ndarray[dlong, ndim=1] cflatten = combo2.flatten() #all faces in one array 
    cdef np.ndarray[dlong, ndim=1] tempidx = np.argsort(cflatten).astype(npy_long,order='C') 
    
    cdef dlong[:] csort = cflatten[tempidx] # sorted faces (for second matrix)
    cdef Py_ssize_t o

    cdef np.ndarray[dlong, ndim=1] repeats = np.zeros((numel1,),dtype=npy_long ,order='C') 
    cdef dlong[:] repeatv = repeats    
    
    for i in range(numel1):
        o = bisectionSearch(csort,combov1[i]) # only works on a sorted array  
        #if the item is found then o!=-1 and it must be a repeat 
        if o!=-1:
            repeatv[i] = 1
    
    return repeats 

@cython.boundscheck(False)
@cython.wraparound(False)
def boundcall(dlong[:,:] tconnection, dlong[:,:] fconnection, double[:,:] node):
    """
    Compute cell boundary conditions for E4D input. Works for gentle topography. 
    Assumes sides of mesh are flat. 

    Parameters
    ----------
    tconnection: nd array 
        Truncated connection matrix, should be indexed with the output of faces3d. 
        N by 4. 
    fconnection: nd array 
        Connection matrix of outer faces (same number of entries (N) as 
        tconnection). N by 3 
    node: nd array 
        Node matrix. N by 3. 

    Raises
    ------
    ValueError
        If length of tconnection and fconnection are not the same. 

    Returns
    -------
    node_bd : nd array 
        1 dimensional array describing node boundary flags 
    face_bd : nd array 
        1 dimensional array describing face boundary flags 

    """
    
    #declare vars 
    cdef int numel = tconnection.shape[0]
    cdef int numel_check = fconnection.shape[0]
    cdef int numnp = node.shape[0]
    
    #face cell handling 
    cdef int npere_face = 3
    cdef double[:] x = np.zeros(npere_face ,dtype=float)
    cdef double[:] y = np.zeros(npere_face ,dtype=float)
    cdef double[:] z = np.zeros(npere_face ,dtype=float)
    cdef double dx, dy, dz
    
    #tetrahedral cell handling 
    cdef int npere = 4
    cdef double[:] cx = np.zeros(npere ,dtype=float)
    cdef double[:] cy = np.zeros(npere ,dtype=float)
    cdef double[:] cz = np.zeros(npere ,dtype=float)
    
    cdef double[:] tcentriod = np.zeros(npere_face,dtype=float)
    cdef double[:] fcentriod = np.zeros(npere_face,dtype=float)
    
    cdef dlong bd # boundary flags 
    cdef np.ndarray[dlong, ndim=1] node_bd = np.zeros(numnp,dtype=npy_long)
    cdef np.ndarray[dlong, ndim=1] face_bd = np.zeros(numel,dtype=npy_long)
    cdef dlong[:] node_bdv = node_bd 
    cdef dlong[:] face_bdv = face_bd 
    
    #do checks? 
    if numel != numel_check:
        raise ValueError('Mismatch in connection matrix lengths')
    
    #compute cell centres 
    for i in range(numel):
        for j in range(npere_face):
            x[j] = node[fconnection[i,j],0]
            y[j] = node[fconnection[i,j],1]
            z[j] = node[fconnection[i,j],2]
        
        dx = fmax(x) - fmin(x)
        dy = fmax(y) - fmin(y)
        
        if dy < 1e-16 or dx < 1e-16 : #face is on the side of the mesh 
            bd = 2
        else:
            fcentriod[0] = mean(x)
            fcentriod[1] = mean(y)
            fcentriod[2] = mean(z)
            for j in range(npere):
                cx[j] = node[tconnection[i,j],0]
                cy[j] = node[tconnection[i,j],1]
                cz[j] = node[tconnection[i,j],2]
            tcentriod[0] = mean(cx)
            tcentriod[1] = mean(cy)
            tcentriod[2] = mean(cz)              
            #if outward vector faces downwards then then the element is on base of mesh 
            if fcentriod[2] > tcentriod[2]:
                bd = 1
            else:
                bd = 2
        for j in range(npere_face):
            if node_bdv[fconnection[i,j]] != 2: # check its not already assiged to side of mesh 
                node_bdv[fconnection[i,j]] = bd 
        
        face_bdv[i] = bd    
        
    return node_bd, face_bd

@cython.boundscheck(False)
@cython.wraparound(False)
def rmRepeatNodes(double[:,:] node, dlong [:,:] connection, int n=10):
    """
    Remove colocated nodes, useful for some mesh formats or remapping the ABMN 
    matrix when colocated electrodes are present (becuase of multiple survey lines). 

    Parameters
    ----------
    node: numpy array
        Mesh node matrix or in this case can be the electrode coordinates. Should
        be N by 3. 
    connection: numpy array 
        Mesh connection matrix.
    n : int, optional
        Maximum number of expected colacated points per node (limited to save on memory). 
        The default is 10. Nb: n can be set to the N (size of node matrix) but this will 
        entail a dlonger processing time. 

    Raises
    ------
    Exception
        if n is exceeded during the colacated point lookup. 

    Returns
    -------
    idx : nd array 
        n by numnp array describing the groupings of repeated nodes at each 
        node coordinate, where -1 means no grouping. 
    nnodes : nd array 
        New nodes locations with repeated (colocated coordinates) filtered out.
    nconnec : nd array 
        New connection matrix, same size as connection but nodes indices have 
        been remapped to ignore repeated points. 
    """
    cdef int numnp = node.shape[0] # number of nodes 
    cdef double[:] x = node[:,0] # extract x y z columns 
    cdef double[:] y = node[:,1]
    cdef double[:] z = node[:,2]
    cdef np.ndarray[double,ndim=2] nnodes = np.zeros((numnp,3), dtype=float) # allocate array for new point array 
    cdef double[:,:] nnodev = nnodes #(memory view)

    cdef int ndim = node.shape[1] # should be 3 (always)
    
    cdef int c = 0 # rolling count used in for loops 
    cdef np.ndarray[dlong,ndim=2] idx = np.zeros((numnp,n), dtype=npy_long)-1 # index grouping 
    cdef dlong[:,:] idxv = idx 
    cdef int err = 0 # error flag
    cdef Py_ssize_t i,j # looping variables 
    
    cdef int numel = connection.shape[0]
    cdef int npere = connection.shape[1]
    cdef np.ndarray[dlong,ndim=2] nconnec = np.zeros((numel,npere), dtype=npy_long)
    cdef dlong[:,:] nconnecv = nconnec 

    
    cdef dlong[:] remap = np.zeros(numnp,dtype=npy_long)-1 # array to allocate remapped node values to 
    cdef dlong[:] flag = np.zeros(numnp,dtype=npy_long) # decides if point already added
    
    for i in range(numnp):
        #this bit can be parallel 
        idx[i,0] = i
        c=1
        for j in range(numnp):
            if i!=j and x[i]==x[j] and y[i]==y[j] and z[i]==z[j]:
                idxv[i,c] = j
                c=c+1
                if c>n:
                    err = 1
                    
    if err==1:
        raise Exception("Maximum number of colocated points exceeded")
        
    c = 0 # reset count to zero 

    for i in range(numnp): # needs to be processed sequentially 
        if flag[i]==0:#point not already flagged 
            for j in range(ndim):    
                nnodev[c,j] = node[i,j]# put new point into array 
            for j in range(n):
                if idx[i,j]==-1:
                    break
                else:
                    flag[idx[i,j]]=1 # remove colocated points from selection pool 
                    remap[idx[i,j]]=c # add count to remapple array 
            c+=1
         
    nnodes = nnodes[0:c,:] # truncate matrix down to the unique points 
    
    #sort remapping of connection matrix 
    for i in range(numel):
        for j in range(npere):
            nconnecv[i,j] = remap[connection[i,j]]
                    
    return idx, nnodes, nconnec # index groupings, new points (filtered) and remapped connection matrix 