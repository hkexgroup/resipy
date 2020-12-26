# compute covariance matrix on Nvidia GPU 
import os
import numpy as np 
import cupy as cp 
#%% functions 
def readSize(fname):
    fh = open(fname)
    line1 = fh.readline().split()
    fh.close()
    return [int(k) for k in line1]

def readJacob(fname):
    array = []
    fh = open(fname)
    header = fh.readline()
    size = [int(k) for k in header.split()]
    lines = fh.readlines()
    fh.close()
    for line in lines:
        array += [float(x) for x in line.split()]
    return tuple(size), array

def readRm(fname, jsize, rsize):
    Rmap = np.genfromtxt(fname.replace('R','Rindex'),
                         skip_header=1,dtype=int)-1
    
    Rvals = np.genfromtxt(fname,
                          skip_header=1)
    
    Rn = np.zeros((jsize[1],jsize[1]))
    for i in range(rsize[0]):
        for j in range(rsize[1]):
            k = Rmap[i,j]
            if k!=-1:
                Rn[i,k] = Rvals[i,j]
    
    R = cp.array(Rn,dtype=np.float32)
    return R 

def getAlpha(fname):
    #return final reported alpha value, file should be the .out file from andy's code
    fh = open(fname,'r')
    lines = fh.readlines()
    fh.close()
    idx = []
    c = 0
    for line in lines:
        if line.find('Alpha:') != -1:
            idx.append(c)
        c+=1
    if len(idx) == 0:
        raise ValueError("Can't find alpha line")
        
    fnl = lines[max(idx)].split()
    alpha = float(fnl[1])
    return alpha


def calc(invdir):
    """Compute Resolution and Covariance matrix for 2D problems using GPU. 

    Parameters
    ----------
    invdir : string 
        Inversion directory used by R2.

    Returns
    -------
    covar : nd array 
        Values along the diagonal of the coviarance matrix.
    remat : nd array 
        Values along the diagonal of the Resolution matrix..

    """
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()
    # read in jacobian
    jsize, jdata = readJacob(os.path.join(invdir,'f001_J.dat'))
    Jn = np.array(jdata,dtype=np.float32).reshape(jsize) # numpy equivalent 
    J = cp.array(Jn,dtype=np.float32)
    
    # read in data Weighting matrix
    protocol = np.genfromtxt(os.path.join(invdir,'f001_err.dat'),
                             skip_header=1)
    
    Wd = cp.array(np.diag(protocol[:,8]),dtype=np.float32)
    
    # read in model roughness matrix 
    rsize = readSize(os.path.join(invdir,'f001_R.dat'))
    
    R = readRm(os.path.join(invdir,'f001_R.dat'), jsize, rsize)
    
    #construct A and b on GPU 
    files = os.listdir(invdir)
    for f in files:
        if f.endswith('.out'):
            alpha = getAlpha(os.path.join(invdir,f))
            break
    S = cp.matmul(cp.matmul(J.T,Wd.T), cp.matmul(Wd,J))
    A = S + alpha*R #Form A (Menke et al, 2015)
    
    #get rid of parameters we dont need anymore to free up memory 
    J = None
    Wd = None
    R = None 
    mempool.free_all_blocks()
    
    Cm = cp.linalg.inv(A) # solve inverse of A to get covariance matrix 
    ResM = Cm*S
    A = None
    mempool.free_all_blocks()
    
    # retrieve outputs as numpy arrays 
    covar = np.diagonal(Cm.get())
    remat = np.diagonal(ResM.get())
    
    #finally clear memory once again 
    Cm = None
    S = None
    ResM = None
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()
    
    return covar, remat
