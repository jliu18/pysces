from __future__ import division
import numpy as np
from pysces import *
import matplotlib.pyplot as plt

#plots the forces/pressures on body for the steady case


alpha_deg = 9 # degrees
Uinfty = (1,0)

num_points = 32
airfoil = naca_airfoil("0012", num_points, zero_thick_te=True)
airfoil = TransformedBody(airfoil, angle=alpha_deg, displacement=(1,2))
panels = BoundVortices(airfoil, Uinfty)

panels.update_strengths()
xvort, gam = panels.vortices.positions, panels.vortices.strengths
q = airfoil.get_points()

vel = force._compute_velocity(panels, Uinfty, wake=None)
#force, Cp = force.compute_force_pressure(bound_old=None, wake_old=None, bound_new=panels, wake_new=None)
v = np.linalg.norm(vel, axis=1)
Cp = 1 - 4*v**2 #extra factor of 2

coll_body = panels.collocation_pts
motion = airfoil.get_motion()
coll = motion.map_position(coll_body)
plt.figure(1)
plt.plot(q[:,0], q[:,1], 'k-')
plt.quiver(coll[:,0], coll[:,1], vel[:,0], vel[:,1])
plt.axis('equal')

#plot pressure
x = panels.collocation_pts
x1, x2 = np.array_split(x, 2)
Cp1, Cp2 = np.array_split(Cp, 2)
#print('The total force is:')
#print force

plt.figure(2)
plt.scatter(x1[:,0], Cp1, c='b', edgecolors='none', label='upper')
plt.scatter(x2[:,0], Cp2, c='r', edgecolors='none', label='lower')
plt.legend()
plt.xlabel('s/c')
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()

plt.figure(3)
gam = panels.vortices.strengths
g1, g2 = np.array_split(gam, 2)
plt.scatter(x1[:,0], g1, c='b', edgecolors='none', label='upper')
plt.scatter(x2[:,0], g2, c='r', edgecolors='none', label='lower')
# reverse g1
g1 = g1[::-1]
plt.scatter(x1[:,0], g1+g2, c='g', edgecolors='none', label='difference')
plt.legend()
plt.xlabel('s/c')
plt.ylabel('gam')
plt.grid(True)
plt.show()

#check Kutta condition
#print panels.vortices.induced_velocity((1,0))