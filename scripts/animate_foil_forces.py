from pysces import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

dt = 0.01
Uinfty = (0,0)
Vortices.core_radius = dt

airfoil = naca_airfoil("0006", 20, zero_thick_te=True) 
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0)) #why no wake for case with no heaving/pitching, bug, also can only do 0 displacement right now
freq = 2*np.pi
airfoil = Pitching(airfoil, 20, freq, phase=90) 
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)
bound = BoundVortices(airfoil)

flow = ExplicitEuler(dt, Uinfty, bound, need_force='wake_impulse')
#flow = ExplicitEuler(dt, Uinfty, bound, need_force=False)

fig, ax = plt.subplots()
ax.axis('equal')
ax.axis([-4, 4,-4,4])
ax.grid(True)
q = airfoil.get_points()
line, = ax.plot(q[:,0], q[:,1], '-k')
maxval = dt
pts = ax.scatter(0, 0, c=0,
                 cmap='bwr', vmin=-maxval, vmax=maxval, edgecolors='none')

def gen_points():
    flow.initialize()
    num_steps = 200
    for i in range(num_steps):
        flow.advance()
        yield airfoil.get_points(), flow.wake.positions, flow.wake.strengths

def redraw(data):
    q, xvort, gam = data
    line.set_data(q[:,0], q[:,1])
    pts.set_offsets(np.array([xvort[:,0], xvort[:,1]]).T)
    pts.set_array(gam)

movie = animation.FuncAnimation(fig, redraw, gen_points, interval=50)
plt.show()