from pysces import *
import matplotlib.pyplot as plt
import numpy as np
#lift should equal 2*pi*alpha where alpha is angle of attack in radians
#=pi^2/90 * angle of attack (deg)

angle = 10
flat_plate = flat_plate(20)
flat_plate = TransformedBody(flat_plate, displacement=(-0.25, 0))
flat_plate = TransformedBody(flat_plate, angle)
bound = BoundVortices(flat_plate)

num_steps = 400
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
expected_Cl = (np.pi * np.pi / 90) * angle
expected = np.array([expected_Cl, expected_Cl])
steps2 = np.array([0, num_steps*dt])
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(steps, 2*f[:,0], c='b', edgecolors='none', label='Cd')
ax1.scatter(steps, 2*f[:,1], c='r', edgecolors='none', label='Cl')
ax1.plot(steps2, expected, c='g', label='expected Cl')
plt.legend(loc='upper left');
plt.grid(True)
plt.show()

