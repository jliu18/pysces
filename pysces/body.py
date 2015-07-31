import numpy as np
from .motion import RigidMotion

__all__ = ['Body', 'cylinder', 'flat_plate', 'naca_airfoil']

class Body(object):
    """Base class for representing bodies
    """
    def __init__(self, points):
        """Create a body with nodes at the given points

        Parameters
        ----------
        points : 2d array, shape (n,2)
            Array of points defining the boundary of the body
            For a closed body, the boundary curve should be positively oriented
            (counter-clockwise around outside of body), starting from trailing
            edge
        """
        self._time = 0
        self._points = points

    @property
    def time(self):
        """The time used to specify the body's motion"""
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

    def get_points(self, body_frame=False):
        return self._points

    def get_body(self):
        """Return the Body object in the body-fixed frame"""
        return self
                   
    def area(self):
        """Returns the area of the body"""
        #see http://paulbourke.net/geometry/polygonmesh/ 
        #"Calculating the area and centroid of a polygon"
        #check for other references?
        points = self._points
        area = 0
        for i in range(-1, len(points)-1):
            p1 = points[i]
            p2 = points[i+1]
            area += p1[0]*p2[1] - p2[0]*p1[1] 
        area = np.absolute(area/2)
        return area
        
        
    def centroid(self):
        """Returns the centroid of the body in the body-fixed frame
            (assuming uniform density)
        """
        points = self._points
        cx = 0
        cy = 0
        for i in range(-1, len(points)-1):
            p1 = points[i]
            p2 = points[i+1]
            cx += (p1[0] + p2[0]) * (p1[0]*p2[1] - p2[0]*p1[1])
            cy += (p1[1] + p2[1]) * (p1[0]*p2[1] - p2[0]*p1[1])    
        A = self.area()
        cx = 1/(6*A) * cx
        cy = 1/(6*A) * cy
        c = np.array((cx, cy))
        return c

def cylinder(radius, num_points):
    """Return a circular Body with the given radius and number of points"""
    th = np.linspace(0, 2 * np.pi, num_points)
    points = radius * np.array([np.cos(th), np.sin(th)]).T
    return Body(points)

def flat_plate(num_points):
    """Return a flat plate with the given number of points"""
    x = np.linspace(1, 0, num_points)
    y = np.zeros_like(x)
    return Body(np.array([x, y]).T)

def naca_airfoil(code, num_points, zero_thick_te=False, uniform=False):
    """Return a NACA 4-digit series airfoil"""
    # extract parameters from 4-digit code
    code_str = "%04d" % int(code)
    if len(code_str) != 4:
        raise ValueError("NACA designation is more than 4 digits")
    max_camber = 0.01 * int(code_str[0])
    p = 0.1 * int(code_str[1])  # location of max camber
    thickness = 0.01 * int(code_str[2:])
    if uniform:
        x = np.linspace(0, 1, num_points)
    else:
        # closer spacing near leading edge
        theta = np.linspace(0, 0.5 * np.pi, num_points)
        x = 1 - np.cos(theta)

    # thickness
    coefs = [-0.1015, 0.2843, -0.3516, -0.1260, 0, 0.2969]
    if zero_thick_te:
        coefs[0] = -0.1036
    y_thick = 5 * thickness * (np.polyval(coefs[:5], x) +
                               coefs[5] * np.sqrt(x))

    # camber
    front = np.where(x <= p)
    back = np.where(x > p)
    y_camber = np.zeros_like(x)
    if p:
        y_camber[front] = max_camber * x[front] / p**2 * (2 * p - x[front])
        y_camber[back] = max_camber * ((1. - x[back])/(1. - p)**2 *
                                       (1 + x[back] - 2 * p))
    x = np.hstack([x[-1:0:-1], x])
    y = np.hstack([y_camber[-1:0:-1] + y_thick[-1:0:-1],
                   y_camber - y_thick])
    return Body(np.array([x, y]).T)



    

def constrain_x():
    pass

def constrain_y():
    pass

def constrain_theta():
    pass


