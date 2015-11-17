from pysces import *
import matplotlib.pyplot as plt
import numpy as np

#define body
#airfoil = flat_plate(10)
airfoil = naca_airfoil("0006", 20)     
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
#define motion of body
freq = 0.3 * 2 * np.pi
airfoil = Pitching(airfoil, 10, freq, phase=90)
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)
bound = BoundVortices(airfoil)

#set up simulation and flow rate
num_steps = 327 
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, bound, need_force='wake_impulse')

#advance simulation
for i in range(1,num_steps):
    flow.advance()

#plot
f = flow.force
steps = np.arange(0, num_steps*dt, dt);

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(steps, 2*f[:,0], c='b', edgecolors='none', label='Cd')
ax1.scatter(steps, 2*f[:,1], c='r', edgecolors='none', label='Cl')
plt.xlabel('time')
plt.legend(loc='upper left');
plt.grid(True)
plt.show()

average_thrust = np.average(flow.force[:,0])
print average_thrust

average_lift = np.average(flow.force[:,1])
print average_lift