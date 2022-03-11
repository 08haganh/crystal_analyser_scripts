# Class for creating minimum volume ellipsod around a set of points. IF the molecule if planar, points will be added above and below the centroid at a 
# distance of 1A
# This is a Python port from a matlab function which uses the Khachiyan algorithm
# https://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python
# Below is the link to the original matlab implementation
# https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

pi = np.pi
sin = np.sin
cos = np.cos

class MinimumVolumeEllipsoid(object):

    def __init__(self,points,tol=0.001):
        """
        Find the minimum volume ellipse.
        Return A, c where the equation for the ellipse given in "center form" is
        (x-c).T * A * (x-c) = 1
        """
        self.points = points
        points = np.asmatrix(points)
        N, d = points.shape
        Q = np.column_stack((points, np.ones(N))).T
        err = tol+1.0
        u = np.ones(N)/N
        while err > tol:
            # assert u.sum() == 1 # invariant
            X = Q * np.diag(u) * Q.T
            M = np.diag(Q.T * la.inv(X) * Q)
            jdx = np.argmax(M)
            step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
            new_u = (1-step_size)*u
            new_u[jdx] += step_size
            err = la.norm(new_u-u)
            u = new_u
        c = u*points
        A = la.inv(points.T*np.diag(u)*points - c.T*c)/d    
        self.A = np.asarray(A)
        self.centre = np.squeeze(np.asarray(c))
        
    def get_components(self):
        # Returns the centre of the ellipsoid, the axes of the ellipsoid, and the corresponding rotation matrix
        _, D, V = la.svd(self.A)    
        return self.centre, 1./np.sqrt(D), V

    def plot_axes(self,ax=None,rotate=True):
        c, radii, rotation = self.get_components()
        if ax is None:
            fig = plt.figure(figsize=[8,5],dpi=100)
            ax = fig.add_subplot(111,projection="3d")
        else:
            pass
        x, y, z = c
        ax.scatter(x,y,z,s=200)
        rx,ry,rz = radii
        x_vect_p = np.array([rx,0,0])
        x_vect_m = np.array([-rx,0,0])
        y_vect_p = np.array([0,ry,0])
        y_vect_m = np.array([0,-ry,0])
        z_vect_p = np.array([0,0,rz])
        z_vect_m = np.array([0,0,-rz])
        if rotate:
            x_vect_p = np.dot(x_vect_p,rotation)
            x_vect_m = np.dot(x_vect_m,rotation)
            y_vect_p = np.dot(y_vect_p,rotation)
            y_vect_m = np.dot(y_vect_m,rotation)
            z_vect_p = np.dot(z_vect_p,rotation)
            z_vect_m = np.dot(z_vect_m,rotation)
        x_vect_p += c
        x_vect_m += c
        y_vect_p += c
        y_vect_m += c
        z_vect_p += c
        z_vect_m += c
        ax.plot(xs=[x_vect_p[0],x_vect_m[0]],ys=[x_vect_p[1],x_vect_m[1]],zs=[x_vect_p[2],x_vect_m[2]],color='blue')
        ax.plot(xs=[y_vect_p[0],y_vect_m[0]],ys=[y_vect_p[1],y_vect_m[1]],zs=[y_vect_p[2],y_vect_m[2]],color='blue')
        ax.plot(xs=[z_vect_p[0],z_vect_m[0]],ys=[z_vect_p[1],z_vect_m[1]],zs=[z_vect_p[2],z_vect_m[2]],color='blue')

        return ax

    def plot(self,ax=None):
        if ax is None:
            fig = plt.figure(figsize=[8,5],dpi=100)
            ax = fig.add_subplot(111,projection="3d")
        else:
            pass
        c, radii, rotation = self.get_components()
        rx, ry, rz = radii
        # Create grid of points
        u, v = np.mgrid[0:2*pi:20j, -pi/2:pi/2:10j]
        # Ellipsoid 
        x = rx*cos(u)*cos(v)
        y = ry*sin(u)*cos(v)
        z = rz*sin(v)
        E = np.dstack([x,y,z])
        # Rotate points
        E = np.dot(E,rotation) + c
        x, y, z = np.rollaxis(E, axis = -1)
        ax.plot_surface(x, y, z, cstride = 1, rstride = 1, alpha = 0.05)

        return ax

        


