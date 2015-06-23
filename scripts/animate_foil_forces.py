from pysces import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#FIX
airfoil = naca_airfoil("0006", 20) # NACA 0012 airfoil with 20 points per side
# airfoil = flat_plate(20)
airfoil = TransformedBody(airfoil, displacement=(-0.25, 0))
freq = 0.3 * 2*np.pi
airfoil = Pitching(airfoil, 20, freq, phase=90)
airfoil = Heaving(airfoil, (0,0.2), freq, phase=0)
#airfoil = XMotion(airfoil, (0,0), (0,0))
bound = BoundVortices(airfoil)


