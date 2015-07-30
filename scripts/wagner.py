#check the results for lift from impulsive start against wagner function
#not working yet
from pysces import *
import matplotlib.pyplot as plt
import numpy as np

angle = 5
flat_plate = flat_plate(20)
flat_plate = TransformedBody(flat_plate, displacement=(-0.25, 0))
flat_plate = TransformedBody(flat_plate, angle)
bound = BoundVortices(flat_plate)

num_steps = 500
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, bound, need_force='wake_impulse')

for i in range(1,num_steps):
    flow.advance()

#plot force
f = flow.force
t = np.arange(0, num_steps*dt, dt);
wagner1 = 1 - 0.165*np.exp(-0.0455*t) - 0.335*np.exp(-0.3*t) #approximation of wagner function
wagner2 = (t+2)/(t+4) #approximation of wagner function
expected_Cl = (np.pi * np.pi / 90) * angle
expected_L = 0.5 * expected_Cl
plt.figure(1)
plt.plot(t, f[:,1]/expected_L, 'ro', label='Lift/Lift_final')
plt.plot(t, wagner1, 'g--', label='Wagner1')
plt.plot(t, wagner2, 'b--', label='Wagner2')
plt.legend(loc='upper left');
plt.grid(True)
plt.show()