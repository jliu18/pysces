"""A module to easily set up and manage a simulation"""
import numpy as np
from .vortex import Vortices
import force
import copy

__all__ = ['ExplicitEuler', 'RungeKutta2', 'RungeKutta4']

class Timestepper(object):
    """Base class for timesteppers for unsteady boundary element simulation"""

    def __init__(self, dt, Uinfty=(1,0), bound=None, wake=None, need_force=False):
        """Initialize a simulation"""
        self._dt = dt
        self._Uinfty = np.array(Uinfty)
        self._bound = bound
        self._has_body = (bound is not None)
        self._need_force = need_force
        if need_force:
            self._force = np.array([0,0], dtype=np.float64)
        self.initialize(wake)

    def initialize(self, wake=None):
        """Initialize a timestepper

        Solve for panel strengths so that surface boundary conditions are
        satisfied, and shed a particle into the wake so that overall circulation
        is zero.

        """
        self._time = 0
        if wake is None:
            self._wake = Vortices()
        else:
            self._wake = Vortices(wake.positions, wake.strengths)
        if self._has_body:
            self._bound.time = 0
            self._bound.update_strengths_unsteady(self._dt, self._Uinfty)
            self._wake.append(*self._bound.get_newly_shed())

    def advance(self, dt=None):
        """Advance the simulation for one timestep"""
        if not dt:
            dt = self._dt
        x = self._wake.positions
        #calculate force on body
#        if self._need_force:
#            vel = self._wake_velocity() #dont compute twice?
#            f = force.compute_force_wake(self.wake, vel)
#            self._force = np.vstack((self._force, f)) 
        v1 = copy.deepcopy(self.bound.vortices) #need deepcopy?
#        #move the wake
        self._wake_advance(x, dt)    # defer to subclass
        v2 = self.bound.vortices
        #f = force.compute_force_bound(v1, v2)
        f = force.compute_force_pressure(self._bound, v1, v2, self._wake)
        #f = force.compute_force_flat_plate(self._bound)
        self._force = np.vstack((self._force, f)) 
        #move the body   
        #self._body_advance(force, dt)   #also defer to subclass? Can we use things like Runge Kutta on this?
        
    @property
    def time(self):
        """Current simulation time"""
        return self._time

    @property
    def bound(self):
        """Body panels used in the simulation"""
        return self._bound

    @property
    def wake(self):
        """Wake vortices used in the simulation"""
        return self._wake

    @property
    def dt(self):
        """Timestep for the simulation"""
        return self._dt
    
    @property
    def force(self):
        """Force on the body"""
        if self._need_force:
            return self._force 
        else:
            print 'No force calculated'
        
    
    def _wake_velocity(self, pos=None, dt=0):
        """Compute the induced velocity at each of the wake vortices

        This is the right-hand side for the timestepper that advances the
        positions of the wake vortices

        Parameters
        ----------
        pos : array, optional
            Array (shape (n,2)) of positions of wake vortices.  Default is the
            current positions of wake vortices.
        dt : float, optional
            Timestep between current simulation time, and time at which the
            velocity is to be computed (default is 0).

        Returns
        -------
        vel : array, shape (n,2)
            Induced velocities at the specified locations of wake vortices.

        Notes
        -----
        If ``pos`` is specified, the wake positions in the Timestepper object
        are updated, and the strengths of bound elements are updated as well to
        satisfy the no-flow-through boundary condition.

        """
        wake = self._wake
        bound = self._bound
        if pos is None:
            pos = wake.positions
            shed = None
        else:
            wake.positions = pos
            if self._has_body:
                # update body position and strengths of surface elements
                bound.time = self._time + dt
                bound.update_strengths_unsteady(dt, self._Uinfty, wake)
                shed = Vortices(*bound.get_newly_shed())
        vel = wake.induced_velocity()
        vel += self._Uinfty
        if self._has_body:
            vel += bound.induced_velocity(pos)
            if shed:
                vel += shed.induced_velocity(pos)
        return vel

    def _update_flow(self, wake_pos, dt):
        """Update the flow with new positions of wake vortices

        Parameters
        ----------
        wake_pos : array
            The new locations of wake vortices
        dt : float
            The amount by which the time should be incremented

        Notes
        -----
        The body motion is updated to the new time, the strengths of the
        bound elements are updated to enforce the no-flow-through boundary
        condition, and a newly shed vortex is added to the wake.

        """
        self._wake.positions = wake_pos
        self._time += dt
        if self._has_body:
            self._bound.time = self._time
            self._bound.update_strengths_unsteady(dt, self._Uinfty, self._wake)
            self._wake.append(*self._bound.get_newly_shed())


class ExplicitEuler(Timestepper):
    """Timestepper using the explicit Euler method"""

    def _wake_advance(self, x, dt):
        vel = self._wake_velocity()
        self._update_flow(x + vel * dt, dt)
        
#    def _body_advance(self, force, dt):
#        pass
        #add stuff here
        #see _update_flow(), also has body motion update stuff

class RungeKutta2(Timestepper):
    """Timestepper using 2nd-order Runge Kutta"""

    def _wake_advance(self, x, dt):
        k1 = self._wake_velocity()
        k2 = self._wake_velocity(x + dt/2 * k1, dt/2)
        self._update_flow(x + dt * k2, dt)

class RungeKutta4(Timestepper):
    """Timestepper using 4th-order Runge Kutta"""

    def _wake_advance(self, x, dt):
        k1 = self._wake_velocity()
        k2 = self._wake_velocity(x + dt/2 * k1, dt/2)
        k3 = self._wake_velocity(x + dt/2 * k2, dt/2)
        k4 = self._wake_velocity(x + dt * k3, dt)
        self._update_flow(x + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4), dt)
