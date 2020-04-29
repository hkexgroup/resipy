# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 15:57:55 2018
test to deal with matplotlib interactive ginput
@author: jamyd91
"""
import matplotlib.pyplot as plt
import numpy as np 

#%% firt block of code 

max_x = 50#input("maximum x value: ")
min_x = 0
max_y = 20#input("maximum Y value: ")
min_y = 0

fh, ax = plt.subplots()
graph = ax.plot([min_y,max_y],[min_y,max_y],'w')

ax.set_ylim([min_y,max_y])
ax.set_xlim([min_x,max_x])

#
print("right click to add point, left click to remove point and middle click to advance")
print("please pick survey topography")
mouse_input = plt.ginput(1000)
coord_mat = np.array(mouse_input)
x_pts = coord_mat[:,0]
y_pts = coord_mat[:,1]
ax.plot(x_pts,y_pts,'kx-')
plt.draw()
no_bounds = 2

for i in range(no_bounds):
    mouse_input = plt.ginput(1000)
    coord_mat = np.array(mouse_input)
    x_pts = coord_mat[:,0]
    y_pts = coord_mat[:,1]
    ax.plot(x_pts,y_pts,'yx-')
    plt.draw()
    

#%%%
import time

import numpy as np
import matplotlib.pyplot as plt


def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

#Define a triangle by clicking three points

plt.clf()
plt.axis([-1., 1., -1., 1.])
plt.setp(plt.gca(), autoscale_on=False)

tellme('You will define a triangle, click to begin')

plt.waitforbuttonpress()

while True:
    pts = []
    while len(pts) < 3:
        tellme('Select 3 corners with mouse')
        pts = np.asarray(plt.ginput(3, timeout=-1))
        if len(pts) < 3:
            tellme('Too few points, starting over')
            time.sleep(1)  # Wait a second

    ph = plt.fill(pts[:, 0], pts[:, 1], 'r', lw=2)

    tellme('Happy? Key click for yes, mouse click for no')

    if plt.waitforbuttonpress():
        break

    # Get rid of fill
    for p in ph:
        p.remove()
        
pressed = plt.waitforbuttonpress() 

#%% test code
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import secrets


def mouse_coord():
    mouse_input = plt.ginput(1000)
    coord_mat = np.array(mouse_input)
    x = coord_mat[:,0]
    y = coord_mat[:,1]
    return x, y 

def write_geom(topo,regions,title='geom input file'):
    fh = open("geom_input.dat",'w')
    fh.write("pyR2 forward model input file\n")
    fh.write("title: %s\n"%title)
    fh.write("dmns: 2D\n")
    fh.write("format: x | y\n")
    fh.write("topography: %i\n"%len(topo[0]))
    for i in range(len(topo[0])):
        fh.write("%f\t%f\n"%(topo[0][i],topo[1][i]))

    fh.write("\nregions: %i \n"%len(regions[0]))
    for i in range(len(regions[0])):
        x_data = regions[0][i][0]
        y_data = regions[1][i][0]
        fh.write("region %i: %i\n"%(i,len(x_data)))
        for k in range(len(x_data)):
            fh.write("%f\t%f\n"%( x_data[k], y_data[k]))
        fh.write("\n")
    fh.close()    

class Button_functions(object):
    #class variables
    topo_x = []
    topo_y = []
    bound_x = []
    bound_y = []
    color_str = ['y','r','g','b','k','m']
    
    def add_topography(self,event):
        ax.plot(self.topo_x,self.topo_y,'wx-',linewidth=2)#write over the last line
        self.topo_x, self.topo_y = mouse_coord()
        ax.plot(self.topo_x,self.topo_y,'gx-')
        plt.draw() 

    def add_boundary(self,event):
        x_pts, y_pts = mouse_coord()
        if len(x_pts)<3:
            print("cant have a polygon with less than 3 control points")
        else:
            self.bound_x.append([x_pts])
            self.bound_y.append([y_pts])
            ax.fill(x_pts,y_pts,secrets.choice(self.color_str))
            plt.draw()  
        
        
    def reset(self,event):
        ax.cla()
        ax.set_ylim([min_y,max_y])
        ax.set_xlim([min_x,max_x])
        self.topo_x = []
        self.topo_y = []
        self.bound_x = []
        self.bound_y = []
        
    def finish(self,event):
        print(self.topo_x,self.topo_y)
        print(self.bound_x,self.bound_y)
        #close the current figure
        plt.close(fig)
        global topo_x, topo_y
        topo_x = self.topo_x
        topo_y = self.topo_y
        write_geom([self.topo_x,self.topo_y],[self.bound_x,self.bound_y])
         
        
    def __str__(self):
        return "Buttons object to be used in conjuction with a figure to create survey geometry"
        
    
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)  

#prelim dimensions
max_x = 50#input("maximum x value: ")
min_x = 0
max_y = 20#input("maximum Y value: ")
min_y = 0

global max_x,max_y,min_x,min_x 

ax.plot([min_y,max_y],[min_y,max_y],'w')

ax.set_ylim([min_y,max_y])
ax.set_xlim([min_x,max_x])


callback = Button_functions()#set call back function 
#define location and widths of buttons 
axtopo = plt.axes([0.1, 0.05, 0.2, 0.075])
axbound = plt.axes([0.35, 0.05, 0.2, 0.075])
axreset = plt.axes([0.6, 0.05, 0.15, 0.075])
ax_exit = plt.axes([0.8, 0.05, 0.15, 0.075])
#set button functions
btopo = Button(axtopo, 'Add topo')
btopo.on_clicked(callback.add_topography)
bbound = Button(axbound, 'Add bound')
bbound.on_clicked(callback.add_boundary)
breset = Button(axreset, 'reset')
breset.on_clicked(callback.reset)
b_exit = Button(ax_exit, 'exit')
b_exit.on_clicked(callback.finish)
#b_exit.disconnect(callback)

def handle_close(evt):
    print('Closed Figure!')
    
fig.canvas.mpl_connect('close_event', handle_close)

plt.show()
    

#use callback to get at the generated variables 

print("press any key to continue")
# plt.waitforbuttonpress() returns true if keyboard is struck, false if mouse is clicked 
#x= callback.topo_x
#y= callback.topo_y


print('do we get here?')

        