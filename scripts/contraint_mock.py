#mock main method for describing motions with constraints 
from pysces import *
import numpy as np

airfoil = naca_airfoil("0006", num_points=20)
#pitching with free swimming in the x direction
fy = ConstPosition(c=0) #why not working?
ftheta = Sinusoidal(amplitude, freq, phase=0.) #in radians
#don't constrain x
airfoil = airfoil.ycontraint(fy) #constrain y and theta
airfoil = airfoil.thetaconstraint(ftheta)

t = 0
dt = 0.01

bound = BoundVortices(airfoil) #panel object
flow = ExplicitEuler(dt, Uinfty, bound, need_force='wake_impulse') #timestepper object
#etc.

#inside timestepper object
# 