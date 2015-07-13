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
#    def compute_force(self): 
#        pass #defer to subclasses
#
#class PressureMethod(Force):

#computes the velocity at each collocation point (should be tangential only)
def _compute_velocity(bound, Uinfty, wake):
    # get collocation points and tangents
    motion = bound._body.get_motion()
    if motion:
        xcoll_inertial = motion.map_position(bound.collocation_pts)
    else:
        xcoll_inertial = bound.collocation_pts
    # velocity induced by wake
    if wake:
        vel = wake.induced_velocity(xcoll_inertial)
    else:
        vel = np.zeros((bound.num_panels, 2))   
    vel += bound.induced_velocity(xcoll_inertial) 
    # assume body is not deforming: only motion is translation/rotation
    if motion:
        vel -= motion.map_velocity(bound.collocation_pts)
    vel += np.array(Uinfty)  
    return vel
#    return np.linalg.norm(vel, axis=1)

#gets the magnitude of kinematic velocity at each collocation point (only for rigid bodies)
def _compute_kin_velocity(bound, Uinfty):
    motion = bound._body.get_motion()
    if motion:
        v = motion.map_velocity(bound.collocation_pts)
    else:
        v = np.array([0,0])
    v += np.array(Uinfty) 
    return np.linalg.norm(v, axis=1)

#calculates the total potential due to the vortices
#at collocation points
def _vortex_potential(bound, wake):
    motion = bound._body.get_motion()
    if motion:
        pos_inertial = motion.map_position(bound.collocation_pts)
    else: 
        pos_inertial = bound.collocation_pts
    phi = bound.vortices.induced_potential(pos_inertial, motion)
    phi += wake.induced_potential(pos_inertial) 
    #need to include incoming flow potential?
    return phi
    
#doesnt work for flat plate/zero thickness 
def compute_force_pressure(bound_old, wake_old, bound_new, wake_new, pref=0, Uinfty=(1,0)): 
    """Finds the forces on the body by integrating the pressure (from unsteady Bernoulli equation)"""
    #q is the local velocity at collocation point (tangential component only)
    motion = bound_new._body.get_motion()    
    if motion:
        norm_inertial = motion.map_vector(bound_new.normals)
    else:
        norm_inertial = bound_new.normals
    vel = _compute_velocity(bound_new, Uinfty, wake_new)
    q = np.linalg.norm(vel, axis=1)
    #vref is the kinematic velocity
    vref = _compute_kin_velocity(bound_new, Uinfty)
    #dPhi = _vortex_potential(bound_new, wake_new) - _vortex_potential(bound_old, wake_old) 
    #dPhi = np.sum(dPhi)
    p = -(q*q)/2 + (vref*vref)/2 
    #- dPhi + pref
    Cp = 1 - (q**2)
    f = p * bound_new.lengths 
    f = f[:,np.newaxis] * norm_inertial
    force = -np.sum(f, axis=0)
    return force, vel 

def _vortex_potential_thin(strengths):
    phi = np.zeros_like(strengths)
    n = len(strengths)
    for i in range(n):
        if i == 0:
            phi[n-1] = strengths[n-1]
        else:
            phi[n-i-1] = phi[n-i] + strengths[n-i]
    return phi       
    
#this is only for a flat plate (zero-thickness airfoil)/debugging purposes (Katz & Plotkin pg. 415)
def compute_force_flat_plate(bound_old, bound_new, wake, Uinfty=(1,0)):
    vel = _compute_velocity(bound_new, Uinfty, wake)
    q = np.linalg.norm(vel, axis=1)
    dPhi = _vortex_potential_thin(bound_new.vortices.strengths) - _vortex_potential_thin(bound_old.vortices.strengths) 
    p = q * bound_new.vortices.strengths / bound_new.lengths + dPhi
    f = p * bound_new.lengths 
    f = f[:,np.newaxis] * bound_new.normals
    force = -np.sum(f, axis=0)
    return force, vel
    
#class ImpulseMethod(Force):
       
def compute_force_wake(wake, induced_velocity):
    """Finds the force on the body by using vortex impulse of the wake"""
    #force is the negative change in impulse of the wake
    #change in impulse is the sum of gamma_i * (induced_velocity_i cross k)
    #take the cross of the induced velocity with the unit vector in the z direction   
    #q = (v, -u) 
    q = np.transpose(np.array([induced_velocity[:,1], -induced_velocity[:,0]]))
    #multiply by the vortex strength
    dI = q * wake.strengths[:,np.newaxis]
    dI = np.sum(dI, axis=0)
    return -dI

#not working yet
def compute_force_bound(bound, vortex_old, vortex_new): #bound vortices
    """Finds the force on the body by using vortex impulse of the bound vortices"""
    dgam = vortex_new.strengths - vortex_old.strengths
    motion = bound._body.get_motion()
    #this doesnt help    
    if motion:
       vnew_inertial = motion.map_position(vortex_new.positions)
       vold_inertial = motion.map_position(vortex_old.positions) #careful! need old motion? 
    dpos = vnew_inertial - vold_inertial 
    a = np.transpose(np.array([vortex_new.positions[:,1], -vortex_new.positions[:,0]]))
    a = a * dgam[:,np.newaxis]
    b = np.transpose(np.array([dpos[:,1], -dpos[:,0]]))
    b = b * vortex_new.strengths[:,np.newaxis]    
    dI = a + b
    tot =  np.sum(vortex_new.strengths)
    dI = np.sum(dI, axis=0) + np.array([0, tot]) #pU*gam? 
    return -dI

#computes virtual momentum of the body, i.e. pressure impulse (Saffman ch. 4)
def compute_virt_force(bound, vortex_old, vortex_new):
    pass