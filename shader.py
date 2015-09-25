from __future__ import division
'''
Created on Sep 9, 2015

@author: andrew.seitz
'''


import matplotlib.pyplot as pyplot
from shapely import geometry 
from descartes.patch import PolygonPatch
import numpy as np
import csv

#===============================================================================
# BASE CLASS
#===============================================================================

class General_Polygon(object):
    """This is a general base class for all other object, holds two import parameters
    a shape object and its transparency (although the __init__ is often over written
    during inheritance. In addition it lies the foundation for plotting various aspects
    of the polygon"""
    
    def __init__(self, transparency = 1):
        self.transparency = transparency
        self.shape = ""
        
    def plot_coords(self, ax, ob, color = '#999999'):
        """plots outer coordinates of a polygon"""
        x, y = ob.xy
        ax.plot(x, y, 'o', color=color, zorder=1)
    
    def plot_shape(self, ax, ob, color = '#00ff00', zorder = 0):
        """plots a patch to fill in the exterior line that makes up a polygon"""
        
        #self.plot_coords(ax, ob.exterior)
        try:
            patch = PolygonPatch(ob, facecolor = color, edgecolor = 'none', zorder = zorder, alpha = 1 - self.transparency)
            ax.add_patch(patch)
        except AssertionError:
            print "AssertionError: Object is either not a polygon or your object did not intersect the panel"
            
    def plot_line(self, ax, ob, color='#000000'):
        """"plots the line around the exterior of a polygon"""
        
        x, y = ob.xy
        ax.plot(x, y, color, linewidth=1, solid_capstyle='round', zorder=1, alpha = 1 - self.transparency)
    
    def v_color(self,ob):
        """logic function to check and see if you created a valid polygon, if valid
        it will chose the color blue to plot with, its not used in this program"""
        
        COLOR = {
         True:  '#6699cc',
         False: '#ff3333'
        }
        #will return color blue if shapely object is valid
        #will return color red if shapely object is invalid
        return COLOR[ob.is_valid]

#===============================================================================
# BASIC_OBJECTS
#===============================================================================

class Circle(General_Polygon):
    """This function that creates a circular object with a center and a radius.
    The center is inputed as a tuple such as (x,y) and the transparency must be a
    number between 0 and 1, with 0 be opaque and 1 being transparent. 
    
    Example:
    
    circle = Circle((1, 2), 200, transparency = 1)
    
    """
    def __init__(self, center, radius, transparency = 0):
        if type(center) != tuple:
            raise "Please enter your center as a tuple: (x,y)"
        
        if transparency > 1 or transparency < 0:
            raise Exception("transparency values must be between 0 and 1 ") 

        self.center = center
        self.radius = radius
        self.transparency = transparency
        self.shape = geometry.Point(center[0],center[1]).buffer(radius)

class Rectangle(General_Polygon):
    """This function takes a set of coordinates to define the four points of a 4 sided polygon, this does
    not technically have to be a square or rectangle as it can be any 4 sided polygon. The 4 coordinates
    are defined in an array as a series of tuples [(lower left), (upper left), (upper right), (lower right), (lower left)].
    The transparency must be a int or float between 0 or 1.
    
    Example:
    
    rect = Rectangle([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)], transparency = .5)
    
    """
    #takes a set of coordinates 
    #take a fill_factor with 0 being transparent and 1 being opaque 
    def __init__(self,coords, transparency = 0):
        if len(coords) != 5:
            raise Exception("Squares require 5 points to be created such as [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]")
        
        if transparency > 1 or transparency < 0:
            raise Exception("transparency values must be between 0 and 1 ") 
        
        self.coords = coords
        self.shape = geometry.Polygon(coords)
        #0 is completely opaque. 1 = complete transmission 
        self.transparency = transparency

class Rectangle2P(General_Polygon):
    """This function creates a rectangle from just two points, these two points will form a line from which a 
    rectangle will be drawn either above or below the line with a specific width. The direction is dependent on
    the right hand rule, meaning, the order you define your points and the value you set the kwarg direct to
    will change the draw direction. The width of the rectangle is inserted as a int or float and the transparency must be
    a value between 0 and 1.
    
    Example:
    
    rect = Rectangle2P([(0, 1), (1, 1)], 200, direct = True, transparency = .2)
    
    """
    def __init__(self, coords, width, direct = False, transparency = 0):
        if len(coords) != 2:
            raise Exception("Function require 2 points to be created such as [(1, 0), (0, 1)]")
        
        if transparency > 1 or transparency < 0:
            raise Exception("transparency values must be between 0 and 1 ") 
        
        self.direct = direct
        self.width = width
        self.coords = self.make_coords(coords[0][0],coords[0][1],coords[1][0],coords[1][1], self.direct, self.width)
        self.shape = geometry.Polygon(self.coords)
        #0 is completely opaque. 1 = complete transmission 
        self.transparency = transparency            
    
    def make_coords(self,x1,y1,x2,y2, direct, width):
        """function that handles the calculation of the right handle rule 
        and the vector math to draw the rectangle at any angle that the two 
        points make.
        """
        v = np.array([x2-x1,y2-y1])
        if direct == True:
            #90 degree rotation of the vector
            transform = np.array([[0,-1], [1,0]])
        else:
            #270 degree rotation of the vector
            transform = np.array([[0,1], [-1,0]])
        
        #find perpendicular vector and normalize it 
        perp = np.dot(transform,v)
        length = np.sqrt(perp[0]**2 + perp[1]**2)
        self.norm = perp / length
        
        self.R1 = (x1,y1)
        self.R2 = (x2,y2)
        self.R3 = (x1 + self.norm[0]*width, y1 + self.norm[1]*width)
        self.R4 = (x2 + self.norm[0]*width, y2 + self.norm[1]*width)
        
        #order the points go into making a rectangle matters 
        if direct == True:
            return [self.R1, self.R3, self.R4, self.R2, self.R1]
        else:
            return [self.R3, self.R1, self.R2, self.R4, self.R3]
        
#===============================================================================
# GRADIENT OBJECTS      
#===============================================================================

class Gradient_Circle(General_Polygon):
    
    """The function creates a series of rings that create a single circle of your define radius with
    a transparency gradient. The functions takes a center, a radius and a gradient. The 
    center is inputed as a tuple such as (x,y), the radius a float or int and the gradient a numpy array between 0 and 1, 
    the number of point defines of the slope of the gradient, and direction of the array such as from 0 to 1 or from 1 to 0, 
    defines the direction of the gradient. The gradient can be as complex or simple as you like, however you will lose resolution
    at a certain level of a gradient. 
    
    Example:
    
    grad_obj = Gradient_Circle((500,500), 300, np.linspace(0,1,30))
    
    """
    
    def __init__(self, center, radius, gradient):
        if type(center) != tuple:
            raise "Please enter your center as a tuple: (x,y)"  
        if np.max(gradient) > 1 or np.min(gradient) < 0:
            raise Exception("gradient must be entered as a np array of transparencies between 0 and 1 corresponding to the slope of the gradient") 
 
        self.center = center
        self.radius = radius
        self.gradient = gradient 
        self.circle_list = []
        self.object_list = []
        
        #define multi circle parameters
        ring_width = self.radius / len(self.gradient)
        #generate initial circles to be subtracted
        for x in range(len(self.gradient)):
            rad = (x+1) * ring_width
            circle = Circle((center[0], center[1]), rad, transparency = self.gradient[x])
            self.circle_list.append(circle)
        
        #add first circle to object list as it is not subtracted away
        self.object_list.append(self.circle_list[0])
        
        for x in range(len(self.gradient) - 1):
            ring = self.circle_list[x+1].shape.difference(self.circle_list[x].shape)
            dummy = General_Polygon()
            dummy.shape = ring
            dummy.transparency = self.circle_list[x+1].transparency
            self.object_list.append(dummy)
                        
class Gradient_Rect(General_Polygon):
    """This function creates a series of rectangles with various transparencies that define a larger 
    single rectangle with a gradient of transparencies. This function takes an initial set of two coordinates
    a width, a gradient, and kwarg direct that defines the draw direction, this work identically to Rectangle2P, 
    as in it follows the right hand rule. The init coords are inputed as an array of tuples such as [(x1,y1),(x2,y2)],  
    the width as a float or int and the gradient a numpy array between 0 and 1, the number of point defines of the slope of the gradient, 
    and direction of the array such as from 0 to 1 or from 1 to 0, defines the direction of the gradient. 
    The gradient can be as complex or simple as you like, however you will lose resolution at a certain level of a gradient. 
    
    Example:
    
    grad_rect = Gradient_Rect([(300,0),(0, 300)], 400, np.linspace(0,1,30), direct = True)
    
    """
    def __init__(self, init_coords, width, gradient, direct = False):
        if len(init_coords) != 2:
            raise Exception("Function require 2 points to be created such as [(1, 0), (0, 1)]")
        if np.max(gradient) > 1 or np.min(gradient) < 0:
            raise Exception("gradient must be entered as a np array of transparencies between 0 and 1 corresponding to the slope of the gradient") 
        
        self.init_coords = init_coords
        self.width = width
        self.gradient = gradient
        self.direct = direct 
        self.object_list = []
        self.counter = len(gradient) 
        self.individual_width = self.width / len(gradient)
        
        self.generate(self.init_coords)
                         
    def generate(self, coords):
        self.counter -= 1 
        current = Rectangle2P(coords, self.individual_width, direct = self.direct, transparency = self.gradient[self.counter])
        self.object_list.append(current)
        if self.counter == 0:
            return False
        else: 
            self.generate([current.R3,current.R4])

class Gradient_Object(General_Polygon):
    """
    A gradient object is a flexible function with various ways to draw a shape.
    This function will take a set of two coordinates to use a the base of rectangle and a width to find the other points.
    A rectangle will be draw and its intersection with the panel will be used as the generated shape.
    A series of rectangles with variable transparencies will be fit into the generated shape to give it a gradient.
    
    Key Parameters:
    gradient: a numpy array between 0 and 1, the number of point defines of the slope of the gradient,
              and direction of the array such as from 0 to 1 or from 1 to 0, defines the direction of the gradient
    direct_grad: defines the direction that the gradient will be drawn within the shape:
                 3 options:
                 - 'up'   : the gradient will face either up or down in reference to the panel
                 - 'right': the gradient will face either left or right in reference to the panel
                 - 'norm' : normal to the face 
                 
    direct: this variable defines to witch side of the init_coords the rectangle will be drawn (there are two options)
             The direction depend on the right hand rule, True will follow the right hand rule and False will be the opposite
             
    The order you define your two initial points also matters in regard to the direct kwarg (follows right hand rule)
    
    Example:
    grad_obj = Gradient_Object(panel, [(600,2033),(945,1800)], 200, np.linspace(1,0,100), direct = False, direct_grad = 'right')
    
    """
    def __init__(self, panel, init_coords, width, gradient, direct = False, direct_grad = 'up'):
        
        if len(init_coords) != 2:
            raise Exception("Function require 2 points to be created such as [(1, 0), (0, 1)]")
        if np.max(gradient) > 1 or np.min(gradient) < 0:
            raise Exception("gradient must be entered as a np array of transparencies between 0 and 1 corresponding to the slope of the gradient") 
        
        accept_grad = ["up", "right", "norm"]
        if direct_grad not in accept_grad:
            raise Exception("Choose from valid gradient directions: 'up', 'right, 'normF, 'normB' ")
        
        self.panel = panel
        self.init_coords = init_coords
        self.width = width
        self.gradient = gradient
        self.direct = direct
        self.direct_grad = direct_grad 
        self.object_list = []
                
        init_shape = Rectangle2P([self.init_coords[0],self.init_coords[1]], self.width, direct = self.direct, transparency = 0)    
        panel_rect = Rectangle([(0,0), (0,self.panel.panel_height), (self.panel.panel_width, self.panel.panel_height), (self.panel.panel_width, 0), (0,0)])
        #poly is our intersection of our inital shape with the panel
        #this is what will be used to intersect with our gradient object
        poly = init_shape.shape.intersection(panel_rect.shape)

        #poly.bounds return a boundary rectangle of any shape in a tuple (xmin,ymin,xmax,ymax)
        minX = poly.bounds[0]
        minY = poly.bounds[1]
        maxX = poly.bounds[2]
        maxY = poly.bounds[3]
        
        if self.direct_grad == 'up':
            width = maxY - minY
            grad_rect = Gradient_Rect([(minX,minY),(maxX, minY)], width, self.gradient, direct = self.direct)
            #must check the draw direction of the gradient rectangle
            #if we draw in the wrong direction, the intersection of the first rectangle 
            #will either be just a point or line (if the x or y coordinates of the Rectangle2P are equal) 
            if grad_rect.object_list[0].shape.intersection(poly).geom_type == 'Point' or grad_rect.object_list[0].shape.intersection(poly).geom_type == 'LineString':
                grad_rect = Gradient_Rect([(minX,minY),(maxX, minY)], width, self.gradient, direct = not self.direct)
                
        if self.direct_grad == 'right':
            width = maxX - minX
            grad_rect = Gradient_Rect([(minX,maxY),(minX,minY)], width, self.gradient, direct = self.direct)
            #must check the draw direction of the gradient rectangle
            #if we draw in the wrong direction, the intersection of the first rectangle 
            #will either be just a point or line (if the x or y coordinates of the Rectangle2P are equal) 
            if grad_rect.object_list[0].shape.intersection(poly).geom_type == 'Point' or grad_rect.object_list[0].shape.intersection(poly).geom_type == 'LineString':
                grad_rect = Gradient_Rect([(minX,maxY),(minX,minY)], width, self.gradient, direct = not self.direct) 
        
        if self.direct_grad == 'norm':
            grad_rect = Gradient_Rect([self.init_coords[0],self.init_coords[1]], width, self.gradient, direct = self.direct)
        
        #dummy objects are created so that each new created shapley shape
        #will have the draw properties of the base class       
        for item in grad_rect.object_list:
            dummy = General_Polygon(transparency = item.transparency)
            dummy.shape = item.shape.intersection(poly)
            if dummy.shape.is_empty:
                pass
            else:
                self.object_list.append(dummy)

#===============================================================================
# INTERFACE CLASSES: where the magic happens                       
#===============================================================================

class Panel(object):
    """This class serves to represent and create a panel comprised of rectangles to define
    every cell in a panel. The class takes m number of vertical cells and n number of horizontal cells,
    in regards to a tetris panel, but this can be modified for various of panels and orientations.
    You must define the overlap of your cells, the spacing between hypercells (or strings of cells)
    and how you divide a wafer (the size of the wafer is hard coded but this can be changed"
    
    Next a series of rectangles are created to comprise the panel and transmission list, which
    keeps track of how much a single rectangle is shaded by another rectangle. Each rectangle
    is stored in a dictionary with a key comprised of its "x_y" positions in the panel, so the cell 
    in the lower right corner is defined as "0_0".
    
    The panel can be plotted through plot_panel. Object collisions with the panel can be calculated
    using collision, however, this function is limited to case with object colliding with the panel 
    do not overlap with each other.
    
    Lastly, the tranmission list, or record of shading of each cell can be printed and recorded to a 
    csv file in the folder you store this library, this can be change or edited in any way seen fit.
    
    Example:
    
    panel = Panel(83, 6, 1.5, 2, wafer_divisor = 6)
 
    """
    def __init__(self, m, n, overlap, hc_spacing, wafer_divisor = 6):        
        #number of cells in a hc
        self.m = int(m)
        #number of hyper cells
        self.n = int(n)
        #dimensions of the initial wafer
        self.wafer_length = 156 #mm
        self.wafer_height = 156 #mm
        #calculation for exposed cell height
        self.cell_height = (self.wafer_height / wafer_divisor) - overlap 
        #spacing between the hypercells 
        self.hc_spacing = hc_spacing
        #area of one cell
        self.cell_area = self.wafer_length * self.cell_height
        #create a panel
        self.num_cells = self.m * self.n 
        self.object_container = {}
        self.coords_dict = {}
        #panel dimensions
        self.panel_height = self.m * self.cell_height
        self.panel_width = self.n * self.wafer_length + (self.n - 1)*hc_spacing
        
        #collision parameters
        #pocs = percent of cell shaded 
        self.transmission = {}
        for y in range(self.n):
            for x in range(self.m):
                self.transmission["%d_%d" % (x,y)] = 1
        
        #4 corners of rectangle: lower left, lower right, upper left, upper right
        ll = [0,0]
        lr = [0,0]
        ul = [0,0]
        ur = [0,0]
        for y in range(self.n):
            for x in range(self.m):
                ll[0] = y*self.wafer_length + y*self.hc_spacing
                ll[1] = x*self.cell_height
                lr[0] = y*self.wafer_length + self.wafer_length + y*self.hc_spacing
                lr[1] = x*self.cell_height
                ul[0] = y*self.wafer_length + y*self.hc_spacing
                ul[1] = x*self.cell_height + self.cell_height 
                ur[0] = y*self.wafer_length + self.wafer_length + y*self.hc_spacing
                ur[1] = x*self.cell_height + self.cell_height
                
                coords = ([(ll[0],ll[1]),
                           (ul[0],ul[1]),
                           (ur[0],ur[1]),
                           (lr[0],lr[1]),
                           (ll[0],ll[1])])
                
                self.coords_dict["%d_%d" % (x,y)] = coords
                self.object_container["%d_%d" % (x,y)] = Rectangle(coords, transparency = 0)

    def plot_panel(self,ax):
        """Plots the panel"""
        for item in self.object_container:
            self.object_container[item].plot_shape(ax,self.object_container[item].shape)
            self.object_container[item].plot_line(ax,self.object_container[item].shape.boundary)
            
    def collision(self, ob_class, ob):
        """This function will first check and see if a collision happen for an object
        and for an of the cells of panel, and then will calculate the area of cell it
        intersected and its affect on the light reaching the cell"""
        
        #initially check for a collision and if so add it to a new dict
        self.coll = {}
        for item in self.object_container:
            if self.object_container[item].shape.intersects(ob) == True:
                self.coll[item] = self.object_container[item]
         
        #calculate the overlap for the collisions that exist with the panel
        for item in self.coll:
            area_intersected = self.coll[item].shape.intersection(ob).area
            ratio_of_area = area_intersected / self.cell_area
            self.transmission[item] = self.transmission[item] - ratio_of_area * (1- ob_class.transparency)
            
    def transmission_list(self):
        """calculate the transmission list based on the collisions that happened"""
         
        holder = np.empty([self.m,self.n])
        for y in range(self.n):
            for x in range(self.m):
                holder[self.m - x - 1][y] = self.transmission["%d_%d" % (x,y)]
        
        with open('transmission_list.csv','w') as f:
            csvwriter = csv.writer(f, delimiter=',', lineterminator = '\n')
            for row in holder:
                csvwriter.writerow(row)
                        
        print holder
        
class WorkSpace(object):
    """This class serves as a location to store a panel and all the shading objects
    created. In addition, it handle two different types of shading calculations
    and will plot your entire workspace to see visually how the shade is laid out
    against a panel. 
    
    The two main types of collision calculations are geometric and visual: 
        Sum_collisions: handles geometric collisions using shapley and the panel
                        collision function. Main draw back is the function does not 
                        handle when two shaded objects intersect each other, however
                        this function is very fast, so for single objects or for multiple
                        non-overlapping objects it is recommened to use.
        
        Sum_overlap_collisions: visual inspect of collisions with a panel based on
                        matplotlibs ability to plot transparency. Each cell is plotted
                        and then converted a greyscale, summed and averaged to produce
                        a tramission of light for the cell. Handle any case, however,
                        this is a slower routine, might take 1-5 minutes to run.
    
    Example:
    
    panel = Panel(83, 6, 1.5, 2, wafer_divisor = 6)
    grad_obj = Gradient_Circle((500,500), 300, np.linspace(0,1,30))
    grad_obj2 = Gradient_Circle((250,250), 400, np.linspace(0,1,30))
    workspace = WorkSpace(panel)
     
    for item in grad_obj.object_list:
        workspace.add(item)
       
    for item in grad_obj2.object_list:
        workspace.add(item)
    
    workspace.plot_workspace()
    workspace.sum_overlap_collisions()    
    
    """    
    def __init__(self, Panel):
        self.panel = Panel
        self.shapes = []
    
    def add(self, obj):
        self.shapes.append(obj)
    
    def sum_collisions(self):
        #THIS FUNCTION IS DESIGNED FOR CALCUATION of NON-OVERLAPING OBJECTS
        #THESE OBJECT CANNOT BE OVERLAPPING 
        #this FASTER than sum_collisions
    
        #WARNING: use one or the other, please do not use both  
        print "WARNING: objects intersecting is not supported in this function"    
                                           
        for item in self.shapes:
            self.panel.collision(item, item.shape)
        
        print self.panel.transmission_list()
        
    def plot_workspace(self):
        self.fig = pyplot.figure(1, dpi=90)
        self.fig.hold(True)
        ax = self.fig.add_subplot(111)
        ax.set_title('Tetris Panel')
        ax.set_ylim([0, self.panel.panel_height])
        ax.set_xlim([0, self.panel.panel_width])
        ax.set_aspect('equal', adjustable = 'box')
        self.panel.plot_panel(ax)
        
        if len(self.shapes) == 0:
            print "please add some shapes"
        else:
            for item in self.shapes:
                item.plot_shape(ax, item.shape, color = "#000000",zorder=2)
        
        pyplot.show()

    def sum_overlap_collisions(self):
        #THIS FUNCTION IS DESIGNED FOR CALCUATION OVERLAP OF MULTIPLE OBJECTS
        #THESE OBJECT CAN BE OVERLAPPING 
        #this is SLOWER than sum_collisions
    
        #WARNING: use one or the other, please do not use both  
            
        #initial collision detection that just stores the location of cells
        #that have intersections
        coll = {}
        for sh in self.shapes:
            for item in self.panel.object_container:
                if self.panel.object_container[item].shape.intersects(sh.shape) == True:
                    if item not in coll:
                        coll[item] = 1
                    else:
                        pass
                        
        #now cycle through all the cells and do some image processing               
        for key in coll:
            cell = pyplot.figure(1, dpi=10, facecolor= '#FFFFFF')
            ax = cell.add_axes([0,0,1,1])
            if len(self.shapes) == 0:
                print "please add some shapes"
            else:
                for item in self.shapes:
                    item.plot_shape(ax, item.shape, color = "#000000",zorder=2)
            
            #find the coordinates of the cell and set the bounds of the image the coords  
            coords = self.panel.coords_dict[key]
            ylim_low = coords[0][1]
            ylim_high = coords[1][1]
            xlim_low = coords[0][0]
            xlim_high = coords[3][0]
            ax.set_ylim([ylim_low, ylim_high])
            ax.set_xlim([xlim_low, xlim_high])
            ax.set_aspect('auto', adjustable = 'box') 
            pyplot.axis('off')
            #draw the plot, but do not display it
            cell.canvas.draw()
            #store the image as a numpy array in RGB format
            data = np.fromstring(cell.canvas.tostring_rgb(), dtype=np.uint8, sep = '')
            data = data.reshape(cell.canvas.get_width_height()[::-1] + (3,))
            #remove figure from metadata for next round
            pyplot.close(cell)  
            
            #convert image to greyscale and then to a scale from 0 (black) to 1 (white)
            grey_list = []
            for item in data:
                for color in item:
                    r, g, b = color[0], color[1], color[2]
                    grey = (0.2989 * r + 0.5870 * g + 0.1140 * b) / 255
                    grey_list.append(grey)
            
            ave = np.average(grey_list)  
            self.panel.transmission[key] = ave
            
            print "%s : %f" % (key, ave)
                             
        print self.panel.transmission_list()

#===============================================================================
# TEST CODE
#===============================================================================

def main():

    panel = Panel(83, 6, 1.5, 2, wafer_divisor = 6)
    #grad_obj = Gradient_Object(panel, [(600,2033),(945,1800)], 200, np.linspace(1,0,100), direct = False, direct_grad = 'right')
    grad_obj = Gradient_Circle((500,500), 300, np.linspace(0,1,30))
    grad_obj2 = Gradient_Circle((250,250), 400, np.linspace(0,1,30))
 
    workspace = WorkSpace(panel)
     
    for item in grad_obj.object_list:
        workspace.add(item)
       
    for item in grad_obj2.object_list:
        workspace.add(item)
    
    workspace.plot_workspace()
    
    ######### -> use this method for single objects as it is faster 
    #workspace.sum_collisions()
    ######### -> use this method for multiple objects, this is slower but can handle overlapping shapes
    workspace.sum_overlap_collisions()
      
if __name__ == '__main__':
    main()



