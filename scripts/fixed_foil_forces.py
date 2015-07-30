from pysces import *
import matplotlib.pyplot as plt

angle = 5
airfoil = naca_airfoil("0006", 20, zero_thick_te=True)  
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = TransformedBody(airfoil, angle)
bound = BoundVortices(airfoil)

num_steps = 1500
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, bound, need_force='wake_impulse')

for i in range(1,num_steps):
    flow.advance()

#plot force
f = flow.force
steps = np.arange(0, num_steps*dt, dt)
expected_Cl = (np.pi * np.pi / 90) * angle
expected = np.array([expected_Cl, expected_Cl])
steps2 = np.array([0, num_steps*dt])
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(steps, 2*f[:,0], c='b', edgecolors='none', label='Cd')
ax1.scatter(steps, 2*f[:,1], c='r', edgecolors='none', label='Cl')
ax1.plot(steps2, expected, c='g', label='expected Cl')
plt.legend(loc='lower right');
plt.xlabel('time')
plt.grid(True)
plt.show()
