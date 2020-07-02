"""
=====================
Polygon Selector Demo
=====================

Shows how one can select indices of a polygon interactively.

"""
import numpy as np

#from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button
from mpl_toolkits.axes_grid1 import make_axes_locatable

class SelectPoints(object):
    """Select indices from a matplotlib collection using `PolygonSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    pts : array
        Arrays of size Nx2 for N points.
    """

    def __init__(self, ax, pts=None, typ='poly', callback=None):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.pts = pts
        self.vertices = []
        self.path = Path
        self.currentIndex = 1
        if self.pts is not None:
            self.iselect = np.zeros(len(self.pts), dtype=bool)
        self.typ = typ
        self.callback = callback
        self.canvas.mpl_connect('key_press_event', self.keyPress)
        self.cax = self.ax.plot([],[], 'r.-')
        self.caxSelected = self.ax.plot([],[], 'k+')
        self.rect = None
        self.xmin = None
        self.xmax = None
        if self.pts is not None:
            self.xmin = np.nanmin(pts[:,0]) - 10 # padding to get extreme points in the zone
            self.xmax = np.nanmax(pts[:,0]) + 10
        self.drawable = False # can draw on canvas
        
        self.addButtons()
        self.connect()

        
    def addButtons(self):
#        self.polyAxes = self.ax.figure.add_axes([0.6, 0.9, 0.1, 0.08])
#        self.rectAxes = self.ax.figure.add_axes([0.7, 0.9, 0.1, 0.08])
#        self.lineAxes = self.ax.figure.add_axes([0.8, 0.9, 0.1, 0.08])
        divider = make_axes_locatable(self.ax)
        buttons = divider.new_vertical(size='5%', pad=0) # 5% cushion between buttons and plot below it
        axbutton = self.ax.figure.add_axes(buttons)
        axbutton.axis('off')
        self.polyAxes = axbutton.figure.add_axes([0.35, 0.91, 0.1, 0.08]) # just to put the buttons out of legend's way
        self.rectAxes = axbutton.figure.add_axes([0.45, 0.91, 0.1, 0.08])
        self.lineAxes = axbutton.figure.add_axes([0.55, 0.91, 0.1, 0.08])
        self.polyButton = Button(self.polyAxes, 'Poly')
        self.rectButton = Button(self.rectAxes, 'Rect')
        self.lineButton = Button(self.lineAxes, 'Line')
            
        
    def connect(self):
        self.canvas.mpl_connect('button_press_event', self.onPress)
        self.canvas.mpl_connect('motion_notify_event', self.onMove)
#        
#        if self.typ == 'poly':
#            cid = self.canvas.mpl_connect('button_press_event', self.onPressPoly)
#            self.connections = [(cid,)]
#        elif self.typ == 'rect':
#            self.rect = None
#            cid1 = self.canvas.mpl_connect('button_press_event', self.onPressRect),
#            cid2 = self.canvas.mpl_connect('motion_notify_event', self.onMotionRect)
#            self.connections = [cid1, (cid2,)]
#        elif self.typ == 'line':
#            cid1 = self.canvas.mpl_connect('button_press_event', self.onPressLine)
#            cid2 = self.canvas.mpl_connect('motion_notify_event', self.onMotionLine)
#            self.connections = [(cid1,), (cid2,)]
#        else:
#            raise ValueError("Type unknown, choose either 'rect' or 'poly'")
#        
    
    def onPress(self, event):
        if event.inaxes == self.polyAxes:
            self.typ = 'poly'
            self.reset()
            self.drawable = True
        elif event.inaxes == self.rectAxes:
            self.typ = 'rect'
            self.reset()
            self.drawable = True
        elif event.inaxes == self.lineAxes:
            self.typ = 'line'
            self.reset()
            self.drawable = True
        elif event.inaxes == self.ax:
            if self.drawable is True:
                if self.typ == 'poly':
                    self.onPressPoly(event)
                elif self.typ == 'rect':
                    self.onPressRect(event)
                elif self.typ == 'line':
                    self.onPressLine(event)
                
        
    def onMove(self, event):
        if self.drawable is True:
            if self.typ == 'rect':
                self.onMotionRect(event)
            elif self.typ == 'line':
                self.onMotionLine(event)
        
    
    def onPressPoly(self, event):
        if event.button == 1:
            self.vertices.append((event.xdata, event.ydata))
            self.drawPath()
        if event.button == 3:
            self.vertices.append(self.vertices[0])
            self.drawPath()
            self.getPointsInside()
            self.disconnect()
    
    
    def onPressRect(self, event):
        if self.rect is None:
            self.rect = self.ax.plot([],[], 'r-')[0]
            self.vertices.append((event.xdata, event.ydata))
            self.drawPath()
        else:
            self.vertices.append((event.xdata, self.vertices[0][1]))
            self.vertices.append((event.xdata, event.ydata))
            self.vertices.append((self.vertices[0][0], event.ydata))
            self.vertices.append((self.vertices[0][0], self.vertices[0][1]))
            self.drawPath()
            self.getPointsInside()
            self.disconnect()
    
    
    def onMotionRect(self, event):
        if (self.rect is not None) & (event.xdata is not None):
            self.rect.set_xdata([self.vertices[0][0],
                                 event.xdata,
                                 event.xdata,
                                 self.vertices[0][0],
                                 self.vertices[0][0]])
            self.rect.set_ydata([self.vertices[0][1],
                                 self.vertices[0][1],
                                 event.ydata,
                                 event.ydata,
                                 self.vertices[0][1]])
            self.canvas.draw()
    
    
    def onPressLine(self, event):
        if self.rect is None:
            self.rect = self.ax.plot([],[], 'r-')[0]
            self.vertices.append((self.xmin, event.ydata))
            self.vertices.append((self.xmax, event.ydata))
            self.drawPath()
        else:
            self.vertices.append((self.xmax, event.ydata))
            self.vertices.append((self.xmin, event.ydata))
            self.vertices.append((self.vertices[0][0], self.vertices[0][1]))
            self.drawPath()
            self.getPointsInside()
            self.disconnect()
    
    
    def onMotionLine(self, event):
        if (self.rect is not None) & (event.ydata is not None):
            self.rect.set_xdata([self.vertices[0][0],
                                 self.vertices[1][0],
                                 self.xmax,
                                 self.xmin,
                                 self.xmin])
            self.rect.set_ydata([self.vertices[0][1],
                                 self.vertices[1][1],
                                 event.ydata,
                                 event.ydata,
                                 self.vertices[0][1]])
            self.canvas.draw()
    
    
    def drawPath(self):
        xy = np.array(self.vertices)
        self.cax[0].set_xdata(xy[:,0])
        self.cax[0].set_ydata(xy[:,1])
        self.canvas.draw_idle()
    
    
    def setVertices(self, xy):
        self.vertices = xy
        if np.sum(xy[0,:] - xy[-1,:]) != 0:
            self.vertices = np.r_[self.vertices, [self.vertices[0,:]]]
        self.drawPath()
        self.getPointsInside()
        
        
    def getPointsInside(self):
        if self.pts is not None:
            path = Path(self.vertices)
            self.iselect = path.contains_points(self.pts)
            # self.caxSelected[0].set_xdata(self.pts[self.iselect, 0])
            # self.caxSelected[0].set_ydata(self.pts[self.iselect, 1])
            if self.callback is not None:
                self.callback(self.iselect)
            else:
                print(np.sum(self.iselect), 'elements selected')
            self.canvas.draw_idle()
        if self.pts is None and self.callback is not None:
            self.callback()
    
    
    def keyPress(self, event):
        if event.key == 'escape':
            self.disconnect()
        if event.key == 'e':
            self.reset()
            self.connect()
        if event.key == 'l':
            self.typ = 'line'
            self.reset()
            self.connect()
        if event.key == 'r':
            self.typ = 'rect'
            self.reset()
            self.connect()
        if event.key == 't':
            self.typ = 'poly'
            self.reset()
            self.connect()
        
                
        
    def disconnect(self):
        self.drawable = False
        
        
    def reset(self):
        if self.rect is not None:
            self.rect.set_xdata([])
            self.rect.set_ydata([])
            self.rect = None
        self.vertices = []
        if self.pts is not None:
            self.iselect = np.zeros(len(self.pts), dtype=bool)
        self.cax[0].set_xdata([])
        self.cax[0].set_ydata([])
        self.caxSelected[0].set_xdata([])
        self.caxSelected[0].set_ydata([])
        self.canvas.draw()
        


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    grid_x = np.random.randn(10)
    grid_y = np.random.randn(10)
    ax.scatter(grid_x, grid_y)
    pts = np.array([grid_x, grid_y]).T
#    rect = Rectangle([0,0], 1,2, alpha=0.3, color='red')
#    ax.add_artist(rect)
    
    def sayHello(arg):
        print('arg = ', arg)
    selector = SelectPoints(ax, pts, typ='poly', callback=sayHello)
#
#    print("Select points in the figure by enclosing them within a polygon.")
#    print("Press the 'esc' key to start a new polygon.")
#    print("Try holding the 'shift' key to move all of the vertices.")
#    print("Try holding the 'ctrl' key to move a single vertex.")
    plt.show()

#    selector.setVertices(np.array([[-0.5, -0.5],
#                                   [-0.5, 0.5],
#                                   [0.5, 0.5],
#                                   [0.5, -0.5]]))
#    selector.disconnect()

    # After figure is closed print the coordinates of the selected points
    print('\nSelected points:')
    print(selector.pts[selector.iselect, :])
