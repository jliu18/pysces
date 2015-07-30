from pysces import *
import matplotlib.pyplot as plt

Uinfty = (1,0)
R = 0.5 #radius of the cylinder
num_pts = 201
body = cylinder(R, num_pts)
#body = TransformedBody(body, angle=45) #rotated cylinder not working
bound = BoundVortices(body, Uinfty)

bound.update_strengths(Uinfty)

q = body.get_points()

vel = force._compute_velocity(bound, Uinfty, wake=None)
    
coll = bound.collocation_pts
motion = bound._body.get_motion()
if motion:
    coll = motion.map_position(bound.collocation_pts)

#checking the velocity
#plt.figure(1)
#plt.plot(q[:,0], q[:,1], 'k-')
#plt.quiver(coll[:,0], coll[:,1], vel[:,0], vel[:,1])
#plt.axis('equal')
#plt.show()

theta = np.arctan2(coll[:,1], coll[:,0])
theta[ theta<0 ] += 2*np.pi
#Vr = vel[:,0]*np.cos(theta) + vel[:,1]*np.sin(theta)
#Vtheta = -vel[:,0]*np.sin(theta) + vel[:,1]*np.cos(theta) 
#Vtheta = Vtheta / np.linalg.norm(Uinfty) #normalize the velocity
#s = np.linspace(0,2*np.pi,100)
#vel_theor = -2 * np.sin(s)
#plt.figure(2)
#plt.plot(s, vel_theor, label="theoretical")
#plt.plot(theta, Vtheta, 'x', label="actual") #I'm getting half of the velocity I should be getting
#plt.ylabel('V_theta/U')
#plt.xlabel('theta (rad)')
#plt.legend()
#plt.show()

#velocity along the y axis
#p = np.linspace(0,1,500)
#p = np.vstack((np.zeros_like(p),p))
#p = np.transpose(np.array([p[0,:], p[1,:]]))
#v = bound.vortices.induced_velocity(p, motion)
#v += Uinfty
#plt.figure(1)
#plt.plot(p[:,1], v[:,0], 'x')
#y = np.linspace(0,2,100)
#R = R*np.ones_like(y)
#plt.plot(R, y, label='r=R')
#plt.xlabel('r')
#plt.ylabel('Vy')
#plt.legend()
#plt.show()

#checking the potential on the cylinder surface
s = np.linspace(0,2*np.pi,100)
phi_theor = np.linalg.norm(Uinfty) * (2*R) * np.cos(s)
phi = force._vortex_potential(bound)
plt.figure(1)
plt.plot(s, phi_theor)
plt.plot(theta, phi, 'x')
plt.ylabel('potential')
plt.xlabel('theta (rad)')
plt.show()

print force.compute_virt_momentum(bound)

##potential along the +y axis
#p_yaxis = np.linspace(0,1,500)
#p_yaxis = np.vstack((np.zeros_like(p_yaxis),p_yaxis))
#p_yaxis = np.transpose(np.array([p_yaxis[0,:], p_yaxis[1,:]]))
#plt.figure(2)
#phi_y = bound.vortices.induced_potential(p_yaxis, motion=None)
#plt.plot(p_yaxis[:,1], phi_y)
#plt.show()
#
##potential along the +x axis
#x = 5
#xaxis = np.linspace(0,x,500)
#p_xaxis = np.vstack((np.zeros_like(xaxis), xaxis))
#p_xaxis = np.transpose(np.array([p_xaxis[1,:], p_xaxis[0,:]]))
#s = np.linspace(0.001,x,100)
#phi_x_theor = np.linalg.norm(Uinfty) * (s + R**2/s)
#plt.figure(3)
#phi_x = bound.vortices.induced_potential(p_xaxis, motion=None)
#phi_x += np.linalg.norm(Uinfty) * xaxis
#plt.plot(p_xaxis[:,0], phi_x, 'bx')
#plt.plot(s, phi_x_theor, 'g')
#plt.axis([0,x,-1,5])
#plt.show()
