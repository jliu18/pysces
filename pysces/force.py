import numpy as np
#from .timestepper import Timestepper

__all__ = ['PressureMethod', 'ImpulseMethod']

class PressureMethod(object):
    """Finds the forces on the body by integrating the pressure"""


class ImpulseMethod(object):
    """Finds the force on the body by using vortex impulse"""
        
#    def __init__(self):
#        self._force = np.array([0,0,0], dtype=np.float64)
#    
#    @property
#    def force(self):
#        """The force on the body"""
#        return self._force
#        
#    def _compute_impulse(self, wake):
#        impulse = np.array([0,0,0], dtype=np.float64)
#        if wake is None:
#            return impulse
#        else:
#            k = np.array([0, 0, 1], dtype=np.float64)
#            impulse = np.array([0,0,0], dtype=np.float64)
#            for i in range(wake.__len__()):
#                impulse += wake.strengths[i] * np.cross(wake.positions[i], k) #check indexing here
#            return impulse
#    
#    #computes force by using the change in impulse           
#    def compute_force(self, old_wake, new_wake):
#        impulse1 = self._compute_impulse(old_wake)
#        impulse2 = self._compute_impulse(new_wake)
#        self._force = impulse2 - impulse1 #need minus sign?
        
    def compute_force(self, wake):
        dImpulse = np.array([0,0,0], dtype=np.float64) #change in impulse with each time step
        k = np.array([0,0,1], dtype=np.float64)
        #need to calculate the induced velocity, or access it if its stored
        #induced_velocity = 
        for i in range(wake.__len__()):
            dImpulse += wake.strengths[i] * np.cross(induced_velocity[i], k)
        return -dImpulse #the force is the negative change in impulse
