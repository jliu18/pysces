from pysces import *
import matplotlib.pyplot as plt

angle = 0
airfoil = naca_airfoil("0006", 40)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = TransformedBody(airfoil, angle) # rotate by 5 degrees about 1/4 chord

#body = cylinder(0.1, 13)

bound = BoundVortices(airfoil)
#bound = BoundVortices(body)

num_steps = 200
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, bound, need_force=True)

for i in range(1,num_steps):
    flow.advance()

#plot pressure
x = flow._bound.collocation_pts
x1, x2 = np.array_split(x, 2)
Cp = flow._pressure
Cp1, Cp2 = np.array_split(Cp, 2)

plt.figure(1)
plt.scatter(x1[:,0], Cp1, c='b', edgecolors='none', label='upper')
plt.scatter(x2[:,0], Cp2, c='r', edgecolors='none', label='lower')
plt.legend()
plt.xlabel('s/c')
plt.ylabel('Cp')
plt.gca().invert_yaxis()
plt.grid(True)
plt.show()
