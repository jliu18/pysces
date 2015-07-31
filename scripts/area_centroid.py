from pysces import *
import matplotlib.pyplot as plt


#check the area and centroid of a square
sq = np.array([[-1, -1], [-1, 1], [1, 1], [1, -1]])
square = Body(sq)
print square.area() #should be 4
print square.centroid() #should be (0,0)


#check the area and centroid of a triangle
tri = np.array([[0, 0], [1, 0], [0.5, 0.75]])
triangle = Body(tri)
print triangle.area() #should be 0.375
print triangle.centroid() #should be (0.5, 0.25)

#check area and centroid of concave polygon
conc = np.array([[0, 0], [2, 0], [2, 1], [1, 0.25], [0, 1]])
concave = Body(conc)
print concave.area() #should be 1.25
print concave.centroid() #should be (1, 0.35)

#flip concave polygon about x axis
concave_x = np.transpose(np.array([-conc[:,0], conc[:,1]]))
concave_x = Body(concave_x)
print concave_x.area() #should be 1.25
print concave_x.centroid() #should be (1, -0.35)

#flip concave polygon about y axis
concave_y = np.transpose(np.array([conc[:,0], -conc[:,1]]))
concave_y = Body(concave_y)
print concave_y.area() #should be 1.25
print concave_y.centroid() #should be (-1, 0.35)