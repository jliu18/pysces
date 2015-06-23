import unittest
import force
from pysces.vortex import Vortices

class TestForce(unittest.TestCase):
        
    def two_vortices_downward(self):
        v1 = (0,0)
        v2 = (1,0)
        g1 = -1
        g2 = 1
        vort = Vortices([v1, v2], [g1, g2])
        vel = vort.induced_velocity
        forces = force.compute_force(vort, vel)
        self.assertEqual(forces, [0,0,0])
        
    def two_vortices_counterrotating(self):
        v1 = (0,0)
        v2 = (1,0)
        g1 = 1
        g2 = 1
        vort = Vortices([v1, v2], [g1, g2])
        vel = vort.induced_velocity
        forces = force.compute_force(vort, vel)
        self.assertEqual(forces, [0,0,0])
        
        