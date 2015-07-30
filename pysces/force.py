import numpy as np

#computes the velocity at each collocation point (should be tangential to body)
def _compute_velocity(bound, Uinfty, wake):
    #velocity induced by the bound vortices (in the body frame)
    vel_body = bound.vortices.induced_velocity(x=bound.collocation_pts)
    #map to the inertial frame
    motion = bound._body.get_motion()
    if motion:
        vel = motion.map_vector(vel_body)
        xcoll_inertial = motion.map_position(bound.collocation_pts)
    else:
        vel = vel_body
        xcoll_inertial = bound.collocation_pts
    #velocity induced by wake (in the inertial frame)
    if wake:        
        vel += wake.induced_velocity(xcoll_inertial)  
    vel += np.array(Uinfty)  
    return vel

#computes the magnitude of kinematic velocity at each collocation point 
def _compute_kin_velocity(bound, Uinfty):
    motion = bound._body.get_motion()
    if motion:
        v = motion.map_velocity(bound.collocation_pts)
    else:
        v = np.array([0,0])
    v += np.array(Uinfty) 
    return np.linalg.norm(v, axis=1)

#computes the angle of the collocation points relative to the horizontal 
#(in the body frame)   
def _compute_alpha(bound):
    xcoll = bound.collocation_pts
    alpha = np.arctan2(xcoll[:,1], xcoll[:,0])
    alpha[ alpha<0 ] += 2*np.pi
    return alpha

#calculates the total potential due to the vortices at collocation points
def _vortex_potential(bound, wake=None, Uinfty=(1,0)):
    #in the body frame
    alpha = _compute_alpha(bound)
    phi = bound.vortices.induced_potential(bound.collocation_pts, alpha)
    #wake not considered yet, may be problematic for not crossing branch cuts
    #see "induced_potential" in .vortex
    motion = bound._body.get_motion()
    if motion:
        pos_inertial = motion.map_position(bound.collocation_pts)
    else: 
        pos_inertial = bound.collocation_pts
    phi+= np.linalg.norm(Uinfty) * pos_inertial[:,0]
    return phi
    
#for airfoil with some thickness
def compute_force_pressure(bound_old, wake_old, bound_new, wake_new, pref=0, Uinfty=(1,0)): 
    """Finds the forces on the body by integrating the pressure (from unsteady Bernoulli equation)"""
    #q is the local velocity at collocation point 
    vel = _compute_velocity(bound_new, Uinfty, wake_new)
    q = np.linalg.norm(vel, axis=1)
    #vref is the kinematic velocity
    vref = _compute_kin_velocity(bound_new, Uinfty)
    if bound_old and wake_old:
        dPhi = _vortex_potential(bound_new, wake_new) - _vortex_potential(bound_old, wake_old) 
        dPhi = np.sum(dPhi)
    else:
        dPhi = 0 
    p = -(q**2)/2 + (vref**2)/2 - dPhi + pref
    Cp = 1 - ((2*q)**2) #factor of two due to singularity in velocity caluclation
    #note: look at Saffman ch. 2 for better calculation of velocity at surface?
    f = p * bound_new.lengths 
    motion = bound_new._body.get_motion()    
    if motion:
        norm_inertial = motion.map_vector(bound_new.normals)
    else:
        norm_inertial = bound_new.normals
    f = f[:,np.newaxis] * norm_inertial
    force = -np.sum(f, axis=0)
    return force, Cp 

def _vortex_potential_thin(strengths):
    phi = np.zeros_like(strengths)
    n = len(strengths)
    for i in range(n):
        if i == 0:
            phi[n-1] = strengths[n-1]
        else:
            phi[n-i-1] = phi[n-i] + strengths[n-i]
    return phi       
    
#for zero-thickness airfoil/flat plate (Katz & Plotkin pg. 415)
def compute_force_flat_plate(bound_old, bound_new, wake, Uinfty=(1,0)):
    vel = _compute_velocity(bound_new, Uinfty, wake)
    q = np.linalg.norm(vel, axis=1)
    if bound_old:
        dPhi = _vortex_potential_thin(bound_new.vortices.strengths) - _vortex_potential_thin(bound_old.vortices.strengths)
    else:
        dPhi = 0
    p = q * bound_new.vortices.strengths / bound_new.lengths + dPhi
    Cp = 1 - (q**2)
    f = p * bound_new.lengths 
    motion = bound_new._body.get_motion()
    if motion:
        norm_inertial = motion.map_position(bound_new.normals)
    else:
        norm_inertial = bound_new.normals
    f = f[:,np.newaxis] * norm_inertial
    force = -np.sum(f, axis=0)
    return force, Cp
    
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
def compute_force_bound(bound_old, bound_new, Uinfty=(1,0)): #bound vortices
    """Finds the force on the body by using vortex impulse of the bound vortices"""
    dgam = bound_new.vortices.strengths - bound_old.vortices.strengths
    motion = bound_new._body.get_motion()  
    if motion:
       vnew_inertial = motion.map_position(bound_new.vortices.positions)
       vold_inertial = motion.map_position(bound_old.vortices.positions) #need old motion? 
    dpos = vnew_inertial - vold_inertial 
    a = np.transpose(np.array([bound_new.vortices.positions[:,1], -bound_new.vortices.positions[:,0]]))
    a = a * dgam[:,np.newaxis]
    b = np.transpose(np.array([dpos[:,1], -dpos[:,0]]))
    b = b * bound_new.vortices.strengths[:,np.newaxis]    
    dI = a + b
    tot =  np.sum(bound_new.strengths)
    Uinfty = np.array(Uinfty)
    dI = np.sum(dI, axis=0) +  np.transpose(np.array([Uinfty[:,1], -Uinfty[:,0]]))* tot[:,np.newaxis] #last term is rho*Uxgam? 
    return -dI

#computes pressure impulse, i.e. virtual momentum (Saffman ch. 4)
#added mass is the change in pressure impulse Ip?
#not tested yet
def compute_virt_momentum(bound, wake=None):
    phi = _vortex_potential(bound, wake)
    motion = bound._body.get_motion()
    if motion:
        norm_inertial = motion.map_position(bound.normals)
    else:
        norm_inertial = bound.normals
    Ip = phi[:,np.newaxis] * norm_inertial
    Ip = np.sum(Ip, axis=0)
    return Ip