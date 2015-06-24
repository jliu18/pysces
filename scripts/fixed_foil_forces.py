from pysces import *
import matplotlib.pyplot as plt

airfoil = naca_airfoil("0006", 20)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
#airfoil = TransformedBody(airfoil, angle=0) # rotate by 0 degrees about 1/4 chord
airfoil = TransformedBody(airfoil, angle=10) # rotate by 5 degrees about 1/4 chord
bound = BoundVortices(airfoil)

num_steps = 200
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = RungeKutta2(dt, Uinfty, bound, need_force=True)

for i in range(1,num_steps):
    flow.advance()

#print flow.force

#plot force
f = flow.force
steps = np.arange(0, num_steps*dt, dt);
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(steps, 2*f[:,0], c='b', edgecolors='none', label='Cd')
ax1.scatter(steps, 2*f[:,1], c='r', edgecolors='none', label='Cl')
plt.legend(loc='upper left');
plt.grid(True)
plt.show()
