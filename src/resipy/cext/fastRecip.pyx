# methods using C extension
cimport cython  #import relevant modules
import numpy as np
cimport numpy as np
from cpython cimport array 
# from libc.stdlib cimport malloc, free

cdef extern from "math.h" nogil:
    cpdef double acos(double x)
cdef extern from "math.h" nogil:
    cpdef double atan(double x)
cdef extern from 'math.h' nogil: # get c square root function
    cpdef double sqrt(double x)
    
@cython.boundscheck(False)#speeds up indexing, however it can be dangerous if the indexes are not in the range of the actual array, 
@cython.wraparound(False)        
cdef int bisection_search(long[:] arr, long var) nogil:
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
cdef double fmax(double[:] arr) nogil: # import note: You need to pass the memory view
    # to a nogil function like this otherwise you get instability 
    cdef int n = arr.shape[0]
    cdef double tmp = arr[0]
    cdef int i 
    for i in range(1,n):
        if tmp<arr[i]:
            tmp = arr[i]
    return tmp

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

cdef long long mergeInt(int a, int b, int pad): #merge 2 ints 
    return np.floor(a*10**pad + b) # merge a and b
    

def fastReciprocals(long[:,:] conf, double[:] r, double[:] p): 
    """Compute reciprocal measurements with compiled code.
    
    Parameters
    ----------
    conf: nd array 
        N by 4 Array of ints
    r: nd array 
        N by 1 array of floats (resistance measurements)
    p: nd array 
        N by 1 array of floats (phase measurements)
    forceSign: long
        Boolian operator to force reciprocals to have the same polarity in the
        calculation of the reciprocals. If 1, reciprocal pairs are forced to 
        have the same polarity regarding reciprocal error calculations. By 
        default forceSign == 0 (no forcing applied). 
    Notes
    -----
    The method first sorts the dipole AB and MN. Then efficiently searches
    for reciprocal pairs with a bisection search. 
    """
    # declare function variables 
    cdef int ndata = conf.shape[0] # number of measurements 
    cdef int conf_check = conf.shape[1]
    cdef int max_elec = 0 
    cdef int i,j,k # loop variables 

    if conf_check != 4: # check conf is of expected length 
        raise Exception('conf should be 4 by N matrix')
    if ndata != r.shape[0]:
        raise Exception('Number of resistance measurements doesnt match configuration length')
    if ndata != p.shape[0]:
        raise Exception('Number of phase measurements doesnt match configuration length')

    # setup arrays for managing internal workings of function 
    for i in range(4): 
        if max(conf[:,i]) > max_elec:
            max_elec = max(conf[:,i])
    
    cdef int pad = len(str(max_elec)) # order of magnitude of max electrode number 
    cdef long[:,:] ab = conf[:,0:2].copy() # ab configuration 
    cdef long[:,:] mn = conf[:,2:4].copy() # mn configuration 

    cdef long[:] AB = np.zeros(ndata,dtype=int) # holds AB combination as one integer 
    cdef long[:] MN = np.zeros(ndata,dtype=int) # holds MN combination as one integer 
    cdef long[:] comboF = np.zeros(ndata,dtype=int) # forward combination
    cdef long[:] comboR = np.zeros(ndata,dtype=int) # reverse combination 
    
    cdef np.ndarray[long, ndim=1] ifwd = np.zeros(ndata,dtype=int)
    cdef long[:] ifwdv = ifwd # flag for forward and reverse direction 
    cdef double rj, pj # recipocal resistance and phase 
    cdef double re, rer, rm, pe # recip error, relative error, reciprocal mean, phase error  

    # outputs of function (so need to be in python/numpy type format)
    cdef np.ndarray[long, ndim=1] irecip = np.zeros(ndata,dtype=int) # index of reciprocal measurement  
    cdef np.ndarray[double, ndim=1] reciprocalErr = np.zeros(ndata,dtype=float)*np.nan # raw reciprocal error 
    cdef np.ndarray[double, ndim=1] reciprocalErrRel = np.zeros(ndata,dtype=float)*np.nan # relative error 
    cdef np.ndarray[double, ndim=1] reciprocalMean = np.zeros(ndata,dtype=float)*np.nan # mean of 2 measurements 
    cdef np.ndarray[double, ndim=1] reciprocalPhase = np.zeros(ndata,dtype=float)*np.nan # mean of 2 measurements 
    cdef long[:] irecipv = irecip 
    cdef double[:] reciprocalErrv = reciprocalErr
    cdef double[:] reciprocalErrRelv = reciprocalErrRel
    cdef double[:] reciprocalMeanv = reciprocalMean
    cdef double[:] reciprocalPhasev = reciprocalPhase 
    
    #merge integers to create unique codes 
    for i in range(ndata):
        sortInt(ab[i,:],2) # sort ints 
        sortInt(mn[i,:],2)
        AB[i] = mergeInt(ab[i,0], ab[i,1],pad)
        MN[i] = mergeInt(mn[i,0], mn[i,1],pad)
        comboF[i] = mergeInt(AB[i], MN[i], pad*2)
        comboR[i] = mergeInt(MN[i], AB[i], pad*2)

    # sort arrays for efficient lookup 
    cdef np.ndarray[long, ndim=1] idxcr = np.argsort(comboR).astype(int,order='C') # (pure python code)
    cdef long[:] idxcrv = idxcr # memory view of sorted index  
    cdef long[:] sortR = np.zeros(ndata,dtype=int)
    for i in range(ndata):
        sortR[i] = comboR[idxcrv[i]]
    
    ## main loop ### 
    for i in range(ndata): 
        # skip the loop if the ifwd value has already been assigned 
        if ifwdv[i] != 0:
            continue 

        k = bisection_search(sortR,comboF[i]) # nb: only works on a sorted array
        # returned index is -1 if no match found. 
        if k != -1: # then match has been found 
            j = idxcrv[k] # index of reverse measurement 
            # store index of forward measurement (also store at reciprocal) 
            irecipv[i] = i+1 
            irecipv[j] = -(i+1)
            # assign forward and reverse directions 
            ifwdv[i] = 1
            ifwdv[j] = -1 
            # get the reciprocal phase and resistance measurements 
            rj = r[j]
            pj = p[j]
            
            # need to assert normal abmn >>> reciprocal mnab, else swap polarity 
            # reciprocal of opposite polarity has configuration nmab or mnba 
            
            if conf[i,0] == conf[j,2] and conf[i,2] == conf[j,1]: 
                # Case abmn >>> nmab 
                rj *= -1 
            elif conf[i,0] == conf[j,3] and conf[i,2] == conf[j,0]: 
                # Case abmn >>> mnba 
                rj *= -1 

            # compute reciprocal error stats  
            re = rj - r[i] # reciprocal error 
            rm = (rj + r[i])/2 # mean of the forward and reciprocal measurements 
            pe = pj - p[i] # phase error (in the case of IP)
            if re == 0 or rm == 0: # if statement needed to avoid 0 float division error 
                rer = 0
            else:
                rer = re/rm # relative error 
            
            # assign to arrays 
            reciprocalErrv[i] = re; reciprocalErrv[j] = -re 
            reciprocalErrRelv[i] = rer; reciprocalErrRelv[j] = rer 
            reciprocalMeanv[i] = rm; reciprocalMeanv[j] = rm 
            reciprocalPhasev[i] = pe; reciprocalPhasev[j] = -pe
            
        else: # there is no pairing 
            ifwdv[i] = 1 # then the measurement is foward only 
            reciprocalMeanv[i] = r[i] # make mean reciprocal equivalent to forward measurement 
    
    return irecip, reciprocalErr, reciprocalErrRel, reciprocalMean, reciprocalPhase, ifwd 


def perElecRi(long[:,:] conf, double[:] recip, long max_elec):
    """
    Average electrode reciprocal error (can be in percent)

    Parameters
    ----------
    conf : nd array 
        4 by N array of electrode configurations 
    recip : nd array 
        Reciprocal error column, can be in percent, all that happens is this
        Column is averaged per an electrode basis.

    Returns
    -------
    riavg: nd array 
        Average reciprocal error per electrode.

    """
    # per electrode reciprocal error 
    cdef int ndata = conf.shape[0] # number of measurements 
    cdef int i, j, a, idx # looping variables
    # cdef int max_elec = 0
    for i in range(4): 
        if max(conf[:,i]) > max_elec:
            max_elec = max(conf[:,i])
            
    cdef np.ndarray[long, ndim=1] elecid = np.arange(1,max_elec+1,dtype=int)
    # note that id numbers are not forced to be consectutive 
    cdef int nelec = len(elecid) # number of electrodes 
    cdef double[:] risum = np.zeros(nelec,dtype=float) # empty array to hold summed contact resistances 
    cdef long[:] ricount = np.zeros(nelec,dtype=int) # hold number of times contact resistance measured on a electrode 
    cdef np.ndarray[double, ndim=1] riavg = np.zeros(nelec,dtype=float) # holds the average reciprocal measurement 
    
    cdef long[:] abmn 

    
    for i in range(ndata): # loop through each measurement to get per measurement contact resistance stats 
        abmn = conf[i,:]
        if recip[i] == -1: # cant average if nan  
            continue 
        # find index of electrode id 
        for j in range(4): 
            a = abmn[j]
            idx = bisection_search(elecid,a) # quick search for the electrode id index
            risum[idx] += recip[i] # add ri to sum of reciprocal errors 
            ricount[idx] += 1 # add one to count of reciprocal errors 
    
    # compute average  reciprocal error for the electrode 
    for i in range(nelec):
        if ricount[i] > 0:# needed to avoid zero division error     
            riavg[i] = risum[i]/ricount[i]
        else:
            riavg[i] = -1 # negative one >> not a number 
        
    return riavg # return elecid and computed averages
