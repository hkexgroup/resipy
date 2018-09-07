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

    def __init__(self, ax, pts, typ='poly', callback=None):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.pts = pts
        self.vertices = []
        self.path = Path
        self.currentIndex = 1
        self.iselect = np.zeros(len(pts), dtype=bool)
        self.typ = typ
        self.connect()
        self.callback = callback
        self.canvas.mpl_connect('key_press_event', self.keyPress)
        self.cax = self.ax.plot([],[], 'r.-')
        self.caxSelected = self.ax.plot([],[], 'k+')
        self.rect = None
        
    def connect(self):
        if self.typ == 'poly':
            cid = self.canvas.mpl_connect('button_press_event', self.onPressPoly)
            self.connections = [(cid,)]
        elif self.typ == 'rect':
            self.rect = None
            cid1 = self.canvas.mpl_connect('button_press_event', self.onPressRect),
            cid2 = self.canvas.mpl_connect('motion_notify_event', self.onMotionRect)
            self.connections = [cid1, (cid2,)]
        else:
            raise ValueError("Type unknown, choose either 'rect' or 'poly'")
            
    def onPressPoly(self, event):
#        print('Poly: mouse press event detected', event.button, event.xdata, event.ydata)
        if event.button == 1:
            self.vertices.append((event.xdata, event.ydata))
            self.drawPath()
        if event.button == 3:
            self.vertices.append(self.vertices[0])
            self.drawPath()
            self.getPointsInside()
            self.disconnect()
    
    def onPressRect(self, event):
#        print('Rect : mouse press event detected', event.button, event.xdata, event.ydata)
        if self.rect is None:
#            print('no rect yet, will add it')
#            self.rect = self.ax.add_artist(Rectangle([event.xdata, event.ydata],
#                                                0, 0, color='red', alpha=0.3))
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
#        print('on motion ...', event.button, event.xdata, event.ydata)
        if (self.rect is not None) & (event.xdata is not None):
#            self.rect.set_width(event.xdata - self.vertices[0][0])
#            self.rect.set_height(event.ydata - self.vertices[0][1])
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
        path = Path(self.vertices)
        self.iselect = path.contains_points(self.pts)
        self.caxSelected[0].set_xdata(self.pts[self.iselect, 0])
        self.caxSelected[0].set_ydata(self.pts[self.iselect, 1])
        if self.callback is not None:
            print('execute callback')
            self.callback(self.iselect)
            print('end of callback')
        else:
            print(np.sum(self.iselect), 'elements selected')
        self.canvas.draw_idle()
    
    def keyPress(self, event):
        if event.key == 'escape':
#            print('disconnected')
            self.disconnect()
        if event.key == 'e':
#            print('connected')
            self.reset()
            self.connect()
        
        
    def disconnect(self):
        for cid in self.connections:
            self.canvas.mpl_disconnect(cid[0])
        
    def reset(self):
        if self.rect is not None:
            self.rect.set_xdata([])
            self.rect.set_ydata([])
        self.vertices = []
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
    
    def sayHello():
        name = input('your name: ')
        print('Hello', name)
    selector = SelectPoints(ax, pts, typ='rect', callback=sayHello)
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
