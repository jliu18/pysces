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

#fix this/check, also consider Uinfty other than (1,0)
#computes the tangential velocity at each collocation point
#figure out which way tangents should point
def _compute_tan_velocity(bound, Uinfty=(1,0), wake=None):
    # get collocation points and normals
    motion = bound._body.get_motion()
    if motion:
        xcoll_inertial = motion.map_position(bound.collocation_pts)
        tangents_inertial = motion.map_vector(bound.tangents)
    else:
        xcoll_inertial = bound.collocation_pts
        tangents_inertial = bound.tangents
    # velocity induced by wake
    if wake:
        vel = wake.induced_velocity(xcoll_inertial)
    else:
        vel = np.zeros((bound.num_panels, 2))
    # assume body is not deforming: only motion is translation/rotation
    if motion:
        vel -= motion.map_velocity(bound.collocation_pts)
    vel += np.array(Uinfty)
    # compute v . t
    return np.sum(vel * tangents_inertial, 1)

#gets the magnitude kinematic velocity at each collocation point (only for rigid bodies)
def _compute_kin_velocity(bound):
    motion = bound._body.get_motion()
    if motion:
        v = motion.map_velocity(bound.collocation_pts)
    else:
        v = np.array([0,0])
    return np.linalg.norm(v, axis=1)

#calculates the total potential due to the vortices
def _vortex_potential(vortex):
   theta = np.arctan2(vortex.positions[:,1], vortex.positions[:,0])
   phi = vortex.strengths * theta / (2 * np.pi)
   return np.sum(phi)

#doesnt work for flat plate/zero thickness yet
def compute_force_pressure(bound, vortex_old, vortex_new, wake, pref=0): #vortex_old, vortex_new refer to bound vortices
    """Finds the forces on the body by integrating the pressure"""
    #q is the local velocity at collocation point (tangential component only)
    q = _compute_tan_velocity(bound, wake=wake)
    #vref is the kinematic velocity
    vref = _compute_kin_velocity(bound)
    #is dPhi right? This is just a constant then for all panels 
    dPhi = _vortex_potential(vortex_new) - _vortex_potential(vortex_old) 
    p = -(q*q)/2 + (vref*vref)/2 - dPhi + pref
    f = p * bound.lengths 
    f = f[:,np.newaxis] * bound.normals
    force = -np.sum(f, axis=0)
    return force
    
#this is only for a flat plate (thin airfoil)/debugging purposes (Katz & Plotkin pg. 415)
def compute_force_flat_plate(bound, vortex_old, vortex_new, wake):
    q = _compute_tan_velocity(bound, wake=wake)
#    dPhi #still need to implement for unsteady motion
    p = q * bound.vortices.strengths / bound.lengths 
    f = p * bound.lengths 
    f = f[:,np.newaxis] * bound.normals
    force = -np.sum(f, axis=0)
    return force
    
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
    
