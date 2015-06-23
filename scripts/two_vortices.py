from pysces import *
#import matplotlib.pyplot as plt
import numpy as np
        
v1 = (0,0)
v2 = (1,0)
s1 = -1
s2 = 1
vort = Vortices([v1,v2],[s1,s2])

dt = 0.1
Uinfty = (0,0)
euler = ExplicitEuler(dt, Uinfty, wake=vort, need_force=True)

num_steps = 200
for i in range(1,num_steps):
    euler.advance()

print euler.force
