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
    return np.sum(vel * tangents_inertial, axis=1)

#gets the magnitude of kinematic velocity at each collocation point (only for rigid bodies)
def _compute_kin_velocity(bound):
    motion = bound._body.get_motion()
    if motion:
        v = motion.map_velocity(bound.collocation_pts)
    else:
        v = np.array([0,0])
    return np.linalg.norm(v, axis=1)

#calculates the total potential due to the vortices
#at collocation points
def _vortex_potential(bound, wake):
    motion = bound._body.get_motion()
    if motion:
        pos_inertial = motion.map_position(bound.collocation_pts)
    else: 
        pos_inertial = bound.collocation_pts
    phi = bound.vortices.induced_potential(pos_inertial) #must map vortices into inertial frame too then?
    phi += wake.induced_potential(pos_inertial) 
    #need to include incoming flow potential?
    return phi

#def _correct_last_panel_length(points):
#    m = points[1] - points[0]
#    m = m[1] / m[0]
#    b = 6.3e-4 - m
#    x = - b / m
#    length = np.sqrt((points[1][0]-x)**2 + points[1][1]**2)
#    return length
    
#doesnt work for flat plate/zero thickness 
def compute_force_pressure(bound_old, wake_old, bound_new, wake_new, pref=0): 
    """Finds the forces on the body by integrating the pressure (from unsteady Bernoulli equation)"""
    #q is the local velocity at collocation point (tangential component only)
    q = _compute_tan_velocity(bound_new, wake=wake_new)
    #vref is the kinematic velocity
    vref = _compute_kin_velocity(bound_new)
    dPhi = _vortex_potential(bound_new, wake_new) - _vortex_potential(bound_old, wake_old) 
    dPhi = np.sum(dPhi)
    p = -(q*q)/2 + (vref*vref)/2 - dPhi + pref
    p = p * bound_new.lengths 
#    p[0] = p[0] * _correct_last_panel_length(bound._body.get_points(body_frame=True))
#    p[bound.num_panels-1] = p[bound.num_panels-1] *  _correct_last_panel_length(bound._body.get_points(body_frame=True))
    f = p[:,np.newaxis] * bound_new.normals
    force = -np.sum(f, axis=0)
    return force


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
def compute_force_flat_plate(bound_old, bound_new, wake):
    q = _compute_tan_velocity(bound_new, wake=wake)
    dPhi = _vortex_potential_thin(bound_new.vortices.strengths) - _vortex_potential_thin(bound_old.vortices.strengths) 
    p = q * bound_new.vortices.strengths / bound_new.lengths + dPhi
    f = p * bound_new.lengths 
    f = f[:,np.newaxis] * bound_new.normals
    force = -np.sum(f, axis=0)
    return force
    
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