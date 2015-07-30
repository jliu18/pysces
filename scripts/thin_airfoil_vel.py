from __future__ import division
import numpy as np
from pysces import *
import matplotlib.pyplot as plt

alpha_deg = 9 # degrees
alpha = alpha_deg * np.pi / 180
Uinfty = (1,0)

num_points = 12
airfoil = naca_airfoil("0012", num_points)
airfoil = TransformedBody(airfoil, angle=alpha_deg, displacement=(1,2))
panels = BoundVortices(airfoil, Uinfty)
panels.update_strengths()
xvort, gam = panels.vortices.positions, panels.vortices.strengths
q = airfoil.get_points()
ds = -np.linalg.norm(np.diff(q, axis=0), axis=1)
dgam = gam / ds
s = np.sqrt(xvort[:,0]**2 + xvort[:,1]**2)

q = airfoil.get_points()
coll_body = panels.collocation_pts
vel_body = panels.vortices.induced_velocity(x=coll_body)
motion = airfoil.get_motion()
coll = motion.map_position(coll_body)
vel = motion.map_velocity(coll_body, vel_body)
vel += Uinfty
plt.figure(1)
plt.plot(q[:,0], q[:,1], 'k-')
plt.quiver(coll[:,0], coll[:,1], vel[:,0], vel[:,1])
plt.axis('equal')
plt.show()

v = np.linalg.norm(vel, axis=1)
Cp = 1 - v**2
Cp1, Cp2 = np.split(Cp, 2)
x = panels.collocation_pts
x1, x2 = np.split(x, 2)
plt.figure(2)
plt.scatter(x1[:,0], Cp1, c='b', edgecolors='none', label='upper')
plt.scatter(x2[:,0], Cp2, c='r', edgecolors='none', label='lower')
plt.legend()
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()

