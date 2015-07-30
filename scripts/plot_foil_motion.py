from pysces import *
import numpy as np
import matplotlib.pyplot as plt


dt = 0.01
Uinfty = (0,0)
Vortices.core_radius = dt

airfoil = naca_airfoil("0006", 20, zero_thick_te=True) 
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0)) #can only do 0 displacement for now
bound = BoundVortices(airfoil)

flow = ExplicitEuler(dt, Uinfty, bound, need_force=False)
pos = airfoil.get_points()[0]

num_steps = 4
for i in range(num_steps):
    flow.advance()
    pos = np.vstack((pos,airfoil.get_points()[0]))

#print flow.wake.positions

#plot position
steps = np.arange(0, num_steps*dt + dt, dt)
x = np.polyfit(steps, pos[:,0], 2)
print x
p = np.poly1d(x)
s = np.linspace(0, 2, 100)
plt.figure(1)
plt.plot(steps, pos[:,0], 'rx', label='position')
plt.plot(s, p(s), 'k', label='quadratic fit')
plt.xlabel('time')
plt.ylabel('position')
plt.legend()
#plt.axis('equal')
plt.show()

#plot velocity
vel = np.diff(pos, axis=0) / dt
steps = np.arange(0, num_steps*dt, dt)
xdot = np.polyfit(steps, vel[:,0], 1)
print xdot
p2 = np.poly1d(xdot)
plt.figure(2)
plt.plot(steps, vel[:,0], 'bx', label='velocity')
plt.plot(s, p2(s), 'k', label='linear fit')
plt.xlabel('time')
plt.ylabel('velocity')
plt.legend()
#plt.axis('equal')
plt.show()
