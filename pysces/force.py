import numpy as np

#__all__ = ['PressureMethod', 'ImpulseMethod']
#
#class Force(object):
#    def __init__(self, force=None):
#        self._force = np.array(force) 
#        
#    @property
#    def force(self):
#        """Force on the body"""
#        return self._force
#        
#    @force.setter
#    def force(self, value):
#        self._force = np.array(value, dtype=np.float64)
#        
#    def append(self, force):
#        force = np.array(force, ndmin=2)
#        if self._force is None:
#            self._force = force
#        else:
#            self._force = np.vstack((self._force, force))
#            
#    def compute_force(self): #do i need this?
#        pass #defer to subclasses
#
#class PressureMethod(Force):
#calculates the total potential due to the vortices
def _vortex_potential(vortex):
   theta = np.arctan2(vortex.positions[:,1], vortex.positions[:,0])
   phi = vortex.strengths * theta / (2 * np.pi)
   return np.sum(phi)

def compute_force_pressure(bound, vortex_old, vortex_new, pref=0):
    """Finds the forces on the body by integrating the pressure"""
    #is this right? This is just a constant then for all panels 
    dPhi = _vortex_potential(vortex_new) - _vortex_potential(vortex_old) 
    p = - dPhi + pref
    return p
    
#class ImpulseMethod(Force):
       
def compute_force_wake(wake, induced_velocity):
    """Finds the force on the body by using vortex impulse of the wake"""
    #force is the negative change in impulse
    #change in impulse is the sum of gamma_i * (induced_velocity_i cross k)
    #take the cross of the induced velocity with the unit vector in the z direction   
    #q = (v, -u) 
    q = np.transpose(np.array([induced_velocity[:,1], -induced_velocity[:,0]]))
    #multiply by the vortex strength
    dI = q * wake.strengths[:,np.newaxis]
    dI = np.sum(dI, axis=0)
    return -dI
    
def compute_force_bound(vortex_old, vortex_new): #bound vortices
    """Finds the force on the body by using vortex impulse of the bound vortices"""
    dgam = vortex_new.strengths - vortex_old.strengths
    dpos = vortex_new.positions - vortex_old.positions #positions need to be in inertial frame?
    a = np.transpose(np.array([vortex_new.positions[:,1], -vortex_new.positions[:,0]]))
    a = a * dgam[:,np.newaxis]
    b = np.transpose(np.array([dpos[:,1], -dpos[:,0]]))
    b = b * vortex_new.strengths[:,np.newaxis]
    dI = a + b
    dI = np.sum(dI, axis=0)
    return -dI
    
