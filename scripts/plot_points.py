from pysces import *
import matplotlib.pyplot as plt
import numpy as np

body = naca_airfoil("00012", 20, zero_thick_te=True)
body = TransformedBody(body, angle=9)

#body = naca_airfoil("0012", 6)

#body = cylinder(0.1, 13)

#body = flat_plate(10)
#body = TransformedBody(body, angle=10)

body_panels = BoundVortices(body)

motion = body_panels._body.get_motion()

vort = motion.map_position(body_panels.vortices.positions)
coll = motion.map_position(body_panels.collocation_pts)
norm = motion.map_vector(body_panels.normals)
tan = motion.map_vector(body_panels.tangents)
foil = body.get_points(body_frame=False)
plt.plot(foil[:,0], foil[:,1], 'k-+')
plt.plot(vort[:,0], vort[:,1], 'ro', label="vortices")
plt.plot(coll[:,0], coll[:,1], 'bx', label="collocation pts")
plt.quiver(coll[:,0], coll[:,1], norm[:,0], norm[:,1])
plt.quiver(coll[:,0], coll[:,1], tan[:,0], tan[:,1])
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
