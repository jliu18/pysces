import numpy as np

__all__ = ['Vortices']

class Vortices(object):
    core_radius = 1.e-3

    def __init__(self, positions=None, strengths=None):
        if positions is None:
            self._positions = None
        else:
            self._positions = np.array(positions, ndmin=2, dtype=np.float64)

        self._circulation = 0
        if strengths is None:
            if positions is None:
                self._strengths = None
            else:
                self._strengths = np.zeros(self._positions.shape[0],
                                           dtype=np.float64)
        else:
            self._strengths = np.array(strengths, ndmin=1, dtype=np.float64)
            self._circulation = np.sum(self._strengths)

    @property
    def positions(self):
        return self._positions

    @positions.setter
    def positions(self, value):
        self._positions = np.array(value, dtype=np.float64)

    @property
    def strengths(self):
        return self._strengths

    @strengths.setter
    def strengths(self, value):
        strengths = np.array(value, ndmin=1, dtype=np.float64)
        self._strengths = strengths
        self._circulation = np.sum(strengths)

    @property
    def circulation(self):
        return self._circulation

    def __len__(self):
        if self._positions is None:
            return 0
        return self._positions.shape[0]

    def __iter__(self):
        if self._positions is None:
            return iter([])
        return iter(zip(self._positions, self._strengths))

    def append(self, position, strength):
        position = np.array(position, ndmin=2)
        strength = np.array(strength, ndmin=1)
        if self._positions is None:
            self._positions = position
            self._strengths = strength
            self._circulation = np.sum(strength)
        else:
            self._positions = np.append(self._positions, position, axis=0)
            self._strengths = np.append(self._strengths, strength)
            self._circulation += np.sum(strength)

    def induced_velocity_single(self, x, xvort, gam):
        """Compute velocity induced at points x by a single vortex

        Parameters
        ----------
        x : 2d array
            Locations at which to compute induced velocity.  Expressed as
            column vectors (i.e., shape should be (n,2))
        xvort : 1d array
            Location of vortex (shape should be (2,))
        gam : float
            Strength of vortex

        Notes
        -----
        Induced velocity is

        .. math:: u_\theta = \frac{\Gamma}{2 \pi r}

        where r is the distance between the point and the vortex.  If this
        distance is less than :class:`core_radius` :math:`r_0`, the velocity is
        regularized as solid-body rotation, with

        .. math:: u_\theta = \frac{\Gamma r}{2\pi r_0^2}
        """
        r = np.array(x, ndmin=2) - np.array(xvort)
        rsq = np.maximum(np.sum(r * r, 1), self.core_radius**2)
        # alternative regularization (Krasny, Eldredge)
        # rsq = np.sum(r * r, 1) + self.core_radius**2
        vel = np.transpose(np.array([-r[:,1], r[:,0]]))
        vel = gam / (2 * np.pi) * vel / rsq[:,np.newaxis]
        return np.squeeze(vel)

    def induced_velocity(self, x=None, motion=None):
        """Compute the induced velocity at the given point(s)"""
        if motion is None:
            positions = self._positions
        else:
            positions = motion.map_position(self._positions)
        if x is None:
            x = self._positions
        else:
            x = np.array(x)
        vel = np.zeros_like(x, dtype=np.float64)
        for xvort, gam in zip(positions, self._strengths):
            vel += self.induced_velocity_single(x, xvort, gam)
        return vel
        
        
#    #still in progress
#    def induced_potential_single(self, x, xvort, gam, alpha=None):
#        """Compute the induced potential at point(s) x due to a single vortex
#        
#        Parameters
#        ----------
#        x : 1d array
#            Locations at which to compute induced potential.  Expressed as
#            column vectors (i.e., shape should be (n,2))
#        xvort : 1d array
#            Location of vortex (shape should be (2,))
#        gam : float
#            Strength of vortex
#        alpha: location of branch cut for each vortex; the potential is 
#            multi-valued, so we make theta (the angle between the vortex and 
#            the x) between [-pi + alpha, pi + alpha]
#            For a convex body, alpha can be defined as the angle of vortex 
#            relative to the horizontal. For a wake vortex, the branch cut 
#            points radially away from the body
#            Unclear for concave body...
#            Need to be careful not to cross the branch, so this might break
#            down for say the multiple body case. 
#
#        """
#        r = np.array(x, ndmin=2) - np.array(xvort)
#        theta = np.arctan2(r[:,1], r[:,0]) 
#        if alpha:
#            theta += alpha
#        phi = self._strengths / (2*np.pi) * theta
#        return phi
#        
#
#    def induced_potential(self, x=None, motion=None):
#        """Compute the induced potential at the given point(s)"""
#        if motion is None:
#            positions = self._positions
#        else:
#            positions = motion.map_position(self._positions)
#        if x is None:
#            x = self._positions
#        else:
#            x = np.array(x)
#        if x.shape==(2,):
#            length = 1
#        else:
#            length = x.shape[0]
#        phi = np.zeros(length)
#        for xvort, gam in zip(positions, self._strengths):
#            phi += self.induced_potential_single(x, xvort, gam) #need way of putting alpha in
#        return phi
    
    #still in progress, seems to work, see cylinder.py
#     #alternative where angle alpha is defined for x
    def induced_potential(self, x, alpha=None, motion=None):
        """Compute the induced potential at point(s) x due to the vortices
        
        Note: this is different from the order that the induced velocity 
              functions compute things
        
        Parameters
        ----------
        x : 1d array
            Locations at which to compute induced potential.  Expressed as
            column vectors (i.e., shape should be (n,2))
        alpha: location of branch for each x; the potential is multi-valued, 
            so we make theta (the angle between the vortex and the x) between
            [-pi + alpha, pi + alpha]
            For a convex body, alpha can be defined as the angle of x 
            relative to the horizontal. May not work adding wake...
            Unclear for concave body... 
            Need to be careful not to cross the branch. 

        """
        if motion is None: 
            positions = self._positions
        else:
            positions = motion.map_position(self._positions)
        if x.shape==(2,):
            length = 1
        else:
            length = x.shape[0]
        phi = np.zeros(length)
        for i in range(length):
            r = x[i] - np.array(positions)
            theta = np.arctan2(r[:,1], r[:,0]) 
            if alpha is not None:
                theta[ theta<(alpha[i]-np.pi) ] += 2*np.pi
            phi_single = self._strengths / (2*np.pi) * theta
            phi_single = np.sum(phi_single)
            phi[i] = phi_single
        return phi
