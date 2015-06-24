from pysces import *
import matplotlib.pyplot as plt
import numpy as np

body = naca_airfoil("0006", 6)
body = TransformedBody(body, angle=10)

#body = naca_airfoil("0012", 6)

#body = cylinder(0.1, 13)

#body = flat_plate(2)
#body = TransformedBody(body, angle=10)

body_panels = BoundVortices(body)

vort = body_panels.vortices.positions
coll = body_panels.collocation_pts
norm = body_panels.normals
tan = body_panels.tangents
foil = body.get_points(body_frame=True)
plt.plot(foil[:,0], foil[:,1], 'k-+')
plt.plot(vort[:,0], vort[:,1], 'ro', label="vortices")
plt.plot(coll[:,0], coll[:,1], 'bx', label="collocation pts")
plt.quiver(coll[:,0], coll[:,1], norm[:,0], norm[:,1])
plt.quiver(coll[:,0], coll[:,1], tan[:,0], tan[:,1])
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
