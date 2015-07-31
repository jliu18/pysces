#mock main method for describing motions with constraints

airfoil = naca_airfoil("0006", num_points) # NACA 0012 airfoil with 20 points per side
airfoil = airfoil.constrain_x(x_displacement)
airfoil = airfoil.constrain_y(y_displacement, freq)
airfoil = airfoil.constrain_theta(angle, freq)



bound = BoundVortices(airfoil) #panel object
flow = ExplicitEuler(dt, Uinfty, bound, need_force=None) #timestepper object
#etc. 