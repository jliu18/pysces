from pysces import *
import matplotlib.pyplot as plt

angle = 9
airfoil = naca_airfoil("0012", 20)   # NACA 2412 airfoil with 20 points per side
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
airfoil = TransformedBody(airfoil, angle) 

bound = BoundVortices(airfoil)

num_steps = 300
Uinfty = (1,0)
dt = 0.01
Vortices.core_radius = dt

flow = ExplicitEuler(dt, Uinfty, bound, need_force='airfoil')

for i in range(1,num_steps):
    flow.advance()

motion = bound._body.get_motion()

vort = motion.map_position(bound.vortices.positions)
coll = motion.map_position(bound.collocation_pts)
v = flow._velocity 
tan = motion.map_vector(bound.tangents)
foil = airfoil.get_points()
plt.plot(foil[:,0], foil[:,1], 'k-+')
plt.plot(vort[:,0], vort[:,1], 'ro', label="vortices")
plt.plot(coll[:,0], coll[:,1], 'bx', label="collocation pts")
wake = flow.wake
plt.plot(wake.positions[:,0], wake.positions[:,1], 'go')
plt.quiver(coll[:,0], coll[:,1], v[:,0], v[:,1])
#plt.quiver(coll[:,0], coll[:,1], tan[:,0], tan[:,1], color='b')
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()