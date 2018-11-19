# create .txt file for Electre 2 software

import numpy as np
import matplotlib.pyplot as plt

def rmnan(x):
    return x[~np.isnan(x)]


def ddskip(nbElectrode=24,spacing=0.25,skip=1,depth=8):
    """
    nbElectrode = 24
    spacing = 0.25
    skip = 1
    depth = 8
    """

    # electrode position
    electrodes = np.zeros((nbElectrode, 3))
    electrodes[:, 0] = np.arange(0, nbElectrode * spacing, spacing)

    # quadrupoles (skip 1 dipole dipole)
    quadrupoles = []
    for i in range(0, nbElectrode):
        for j in range(0, depth):
            a = i + 1
            b = a + 1 + skip
            m = b + 1 + j
            n = m + 1 + skip
            if (b > nbElectrode)|(m > nbElectrode)|(n > nbElectrode):
                break
            quadrupoles.append([a, b, m, n])
    quadrupoles = np.array(quadrupoles, dtype=int)
    print('ddskip : {0:d} quadrupoles generated.'.format(len(quadrupoles)))
    return electrodes,quadrupoles


def write(filename,electrodes,quadrupoles):
    with open(filename, 'w') as f:
        f.write('# X Y Z\n')
        for i in range(0, len(electrodes)):
            f.write('{0:d} {1:f} {2:f} {3:f}\n'.format(i + 1, electrodes[i, 0], electrodes[i, 1], electrodes[i, 2]))
        f.write('# A B M N\n')
        for i in range(0, len(quadrupoles)):
            f.write('{0:d} {1:d} {2:d} {3:d} {4:d}\n'.format(i + 1, quadrupoles[i, 0], quadrupoles[i, 1], quadrupoles[i, 2],
                                                             quadrupoles[i, 3]))


def addReciprocal(array):
    toadd = np.array([])
    counter = 0
    for i in range(0, len(array)):
        rev = [-2, -1, -4, -3]
        index = (array[:, -4:] == array[i, rev]).all(1)
        if len(index[index]) == 0:
            toadd = np.append(toadd, np.hstack(array[i, rev]))
            counter = counter + 1
            print(str(array[i,:])+' add reciprocal : '+str(array[i, rev]))
    l, ll = array.shape

    toadd = toadd.reshape([len(toadd) / ll, ll]) # figure out why this
    # could not be int -> deprecated in the future
    print(str(counter)+ 'reciprocal measurements added.')
    if len(toadd) > 0:
        array = np.vstack((array, toadd))
    return array.astype(int)


def optimize(arg):
    array = np.copy(arg)
    newarray = np.zeros(array.shape, dtype=int)
    c1 = 0  # counter
    for i in range(0, len(array)):
        index = (array[:, 1:2] == array[i, 1:2]).all(1)
        if len(index[index]) > 0:
            c2 = c1 + len(index[index])
            newarray[c1:c2, :] = array[index, :]
            array[index, :] = np.nan  # avoid taking them again
            c1 = c2
    newarray[:, 0] = np.arange(1, len(newarray) + 1)
    print('optimization done.')
    return newarray.astype(int)


def optimize2(arg):
    index = np.lexsort((arg[:, -1], arg[:, -2], arg[:, -3], arg[:, -4]))
    print('optimization2 done.')
    return arg[index, :].astype(int)

    
def optimize3(arg,dummy=True,nchannel=10, showOpti=False): # optimize potential electrodes too !
    ''' optimize measurements scheme for Syscal Pro 10 channel
        
    Arguments:
        arg = 4 columns arrays containing current (2 first columns) 
            and potential electrodes (2 last columns)
        nchannel = number of channel for potential electrodes
        dummy = True : adding non informational measurements for filling the
            channels or optimizing them
        subarray if length 11 (because consecutive electrode so only 10
            different electrodes in a subarray)
        
    Returns :
        electrodes : the optimized electrodea array
        counter : estimation of group of measurements that could be perfomed on
            the same time (helps to estimage survey time)
    '''
    array = np.copy(arg)
    array=array.astype(float)
    #newarray = np.zeros(array.shape, dtype=int)
    newarray=[]
    counter=0
    dummyCounter = 0
    for i in range(0,len(array)):
        # no one can match because meas. are set to nan after assigned
        if ~np.isnan(array[i,0]):
            index = (array[:, 0:2] == array[i, 0:2]).all(1) # same current elct.
            if len(index[index]) > 0:
                subarray=np.copy(array[index,:])
                subarraySorted=[]
                for j in range(0,len(subarray)):
                    if len(rmnan(subarray[:,3])) > 0:
                        k=np.where(subarray[:,3]==np.min(rmnan(subarray[:,3])))[0]
                        subarraySorted.append(subarray[k,:])
                        for m in range(0,len(subarray)):
                            p2=subarray[k,3]
                            subarray[k,:]=np.nan
                            index1=subarray[:,2]==p2 # there should be only one
                            if len(index1[index1]) > 0:
                                k=np.where(index1)[0]
                                subarraySorted.append(subarray[k,:])
                            else: # can't find consecutive next potential meas.
                                if dummy&(len(rmnan(subarray))>1)&(len(subarraySorted)%(nchannel-1)>0): 
                                    # it is useless to add dummy at the end
                                    dum=np.zeros((1,4))
                                    dum[0,0:2]=array[i,0:2]
                                    dum[0,2]=p2
                                    dum[0,3]=np.min(rmnan(subarray[:,2]))
                                    #print("dummy quadrupole added : " + str(dum))
                                    dummyCounter = dummyCounter + 1
                                    subarraySorted.append(dum)
                                break
                    else:
                        break
                subarraySorted=np.vstack(subarraySorted)
                #print('subarray='+str(subarraySorted))
                #newarray[c1:c2, :] = np.array(subarraySorted)
                newarray.append(subarraySorted)
                counter=counter+int(len(subarraySorted)/(nchannel+1))
                if np.mod(len(subarraySorted),nchannel+1) != 0:
                    counter = counter + 1
                array[index, :] = np.nan  # avoid taking them again
    newarray=np.vstack(newarray)
    #newarray[:, 0] = np.arange(1, len(newarray) + 1)
    counter = nbChannelGroup(newarray, show=showOpti, nchannel=nchannel)
    print('optimization3 done : {0:d} quadrupoles ({1:d} dummies; {2:d} channel group).'.format(len(newarray), dummyCounter, counter))
    return newarray.astype(int)


def nbChannelGroup(x, show=False, nchannel=10):
    # compute number of channel group
    subarrays = []
    subarray = [x[0,:]]
    k = x[0,3]
    counter = 1
    for i in range(1,len(x)):
        if len(subarray) == nchannel:
            counter = counter + 1
            subarray = np.vstack(subarray)
            subarrays.append(subarray)
            if show:
                print('subarray (' + str(len(subarray)) + ') : \n' + str(subarray))
            subarray = [x[i,:]]
        else:
            if x[i,2] == k:
                subarray.append(x[i,:])
            else:
                counter = counter + 1
                subarray = np.vstack(subarray)
                subarrays.append(subarray)
                if show:
                    print('subarray (' + str(len(subarray)) + ') : \n' + str(subarray))
                subarray = [x[i,:]]
        k = x[i,3]
    if show:
        print('{0:d} channel groups'.format(counter))
    return counter
             
            

def zigzagElec(elecSpacing=2,lineSpacing=2,nelecLine=6,nbElec=96, ravel=True):
    """create zigzag electrode array positions"""
    electrodes=np.zeros((nbElec,3))
    xpos=np.arange(0,nelecLine*elecSpacing,elecSpacing)
    xpos2 = np.zeros(len(xpos)+1)*np.nan
    xpos2[1:] = xpos
    nbLines=int(nbElec/nelecLine)
    ypos=np.arange(0,nbLines*lineSpacing,lineSpacing)
    flip = 1
    for i in range(0,nbLines):
        for j in range(0,nelecLine):
            a=i*nelecLine+j
            electrodes[a,:]=np.array([a+1, xpos2[flip*(j+1)], ypos[i]])
        flip = flip * -1
    return electrodes


def paralleldd(nbElectrode=24, nbLine=2, spacing=0.25, lineSpacing=0.40,
               skip=1, skipLine=0, depth=3):
    ''' generate configuration to take 3D measurements between two lines
    
    TODO skipLine not supported yet
    max two lines (not more plugs for the 48 Syscal)
    '''
    # electrode position
    electrodes = np.zeros((nbElectrode*nbLine, 3))
    for i in range(0,nbLine):
        electrodes[i*nbElectrode:(i+1)*nbElectrode, 0] = np.arange(0, nbElectrode * spacing, spacing)
        electrodes[i*nbElectrode:(i+1)*nbElectrode, 1] = i*lineSpacing

    # quadrupoles (skip 1 dipole dipole)
    quadrupoles = []
    for h in range(0, nbLine-1):
        for i in range(0, nbElectrode):
            # face to face
            a = i + 1
            b = a + 1 + skip
            m = (h+1)*nbElectrode + a 
            n = (h+1)*nbElectrode + b
            aa = h*nbElectrode + a
            bb = h*nbElectrode + b
            #print('face to face : ' + str([aa,bb,m,n]))
            if (bb > (h+1)*nbElectrode):
                break
            quadrupoles.append([aa, bb, m, n])
            for j in range(0, depth):
                # going upwards
                m = (h+1)*nbElectrode + a + j + 1
                n = (h+1)*nbElectrode + b + j + 1
                #print('upwards : ' + str([aa, bb, m, n]))
                if (m > (h+2)*nbElectrode)|(n > (h+2)*nbElectrode):
                    break
                quadrupoles.append([aa, bb, m, n])
                
                # going downwards
                m = (h+1)*nbElectrode + a - j - 1
                n = (h+1)*nbElectrode + b - j - 1
                #print('downwards : ' + str([aa, bb, m, n]))
                if (m < (h+1)*nbElectrode)|(n < (h+1)*nbElectrode):
                    break
                quadrupoles.append([aa, bb, m, n])
    print('paralleldd : {0:d} quadrupoles generated.'.format(len(quadrupoles)))
    quadrupoles = np.array(quadrupoles, dtype=int)
    return electrodes,quadrupoles


def elecStar(nbLine, nbElec, elecSpacing, centerOffset):
    """ build a star design"""
    elec = []
    for i in range(0,nbLine):
        # compute trigonometric angle
        angle = np.pi/nbLine*i
        rpos = centerOffset + np.linspace(0, nbElec/2*elecSpacing, nbElec/2)
        rpos = rpos[::-1]
        xpos = rpos*np.cos(angle)
        ypos = rpos*np.sin(angle)
        elec.append(np.vstack([xpos,ypos]).T)
        # now for the other part of the line
        angle2 = np.pi + np.pi/nbLine*i
        rpos2 = centerOffset + np.linspace(0, nbElec/2*elecSpacing, nbElec/2)
        xpos2 = rpos2*np.cos(angle2)
        ypos2 = rpos2*np.sin(angle2)
        elec.append(np.vstack([xpos2,ypos2]).T)
    tmp = np.vstack(elec)
    elecs = np.zeros((len(tmp), 3))
    elecs[:,:2] = tmp
    return elecs
    
'''
z = elecStar(4,24,1,1)
plt.close('all')
plt.ioff()
fig, ax = plt.subplots()
ax.plot(z[:,0], z[:,1], 'bo')
ax.grid()
#fig.show()
'''

def tofwd3d(quadrupoles, filename):
    """convert to protocol.dat format"""
    x = np.ones((quadrupoles.shape[0], 9), dtype=int)
    x[:,[2,4,6,8]] = quadrupoles[:,[0,1,2,3]]
    x[:,0] = np.arange(1, len(x)+1)
    np.savetxt(filename, x, header=str(x.shape[0]),
    fmt='%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f', comments='')
    


def togeo(elec, filename):
    """convert to Point for .geo files"""
    with open(filename, 'w') as f:
        for i in range(0, len(elec)):
            f.write('Point(' + str(i+1) + ') = {' + str(elec[i,0]) + ', ' 
                    + str(elec[i,1]) + ', ' + str(elec[i,2]) + ', cl1};\n')



#%% boxford 2018-03-01

#elec, quad0 = ddskip(48, 0.25, skip=0, depth=10)
#elec, quad1 = ddskip(48, 0.25, skip=1, depth=9)
#
#qs = []
#qs.append(optimize3(quad0))
#qs.append(optimize3(quad1))
#
#qs.append(optimize3(quad0[:,[2,3,0,1]]))
#qs.append(optimize3(quad1[:,[2,3,0,1]]))
#
#quad = np.vstack(qs)
#
#write('./out/dd48s01rn.txt', elec, quad)
#
#
