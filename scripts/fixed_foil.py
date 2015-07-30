from pysces import *
import matplotlib.pyplot as plt

airfoil = naca_airfoil("0006", 20)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = TransformedBody(airfoil, angle=5) # rotate by 10 degrees about 1/4 chord
bound = BoundVortices(airfoil)

num_steps = 400
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, bound, need_force='None')

for i in range(1,num_steps):
    flow.advance()

vort = flow.wake.positions
gam = flow.wake.strengths
q = airfoil.get_points()
plt.plot(q[:,0], q[:,1], 'k-')
maxval = 0.01
plt.scatter(vort[:,0], vort[:,1], c=gam,
            cmap='bwr', vmin=-maxval, vmax=maxval, edgecolor='none')
plt.axis('equal')
plt.grid(True)
plt.show()
