import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

sin = np.sin
cos = np.cos
tan = np.tan
arcsin = np.arcsin

class Plane():
    # https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    
    def __init__(self,points):
        """
        p, n = planeFit(points)

        Given an array, points, of shape (d,...)
        representing points in d-dimensional space,
        fit an d-dimensional plane to the points.
        Return a point, p, on the plane (the point-cloud centroid),
        and the normal, n.
        """
        if len(points) == 2:
            centre = np.average(points,axis=0)
            print(centre,points)
            points = np.concatenate([points,centre],axis=0)
        if points.shape[0] >= points.shape[1]:
            points = np.vstack([points[:,0],points[:,1],points[:,2]])
        points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
        assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
        self.ctr = points.mean(axis=1)
        x = points - self.ctr[:,np.newaxis]
        M = np.dot(x, x.T) # Could also use np.cov(x) here.
        vect = la.svd(M)[0][:,-1]
        self.a, self.b, self.c = vect
        # ax + by + cz + d = 0
        self.d = (points[0,0]*self.a + points[1,0]*self.b + points[2,0]*self.c)*-1

    def plane_angle(self, plane):
        a1,b1,c1 = self.a,self.b, self.c
        a2,b2,c2 = plane.a,plane.b, plane.c
            
        d = ( a1 * a2 + b1 * b2 + c1 * c2 )
        e1 = np.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
        e2 = np.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
        d = d / (e1 * e2)
        if d >= 0:
            A = np.degrees(np.arccos(min(1,d))) # ensure d is never greater than 1 because can only do arccos[-1,1]
        else:
            A = np.degrees(np.arccos(max(-1,d))) # ensure d is never less than -1 because can only do arccos[-1,1]
        if A > 90:
            A = 180 - A
        return A

    def unit_normal(self):
        mag = np.sqrt(self.a**2 + self.b**2 + self.c**2)
        unit_norm = np.array([self.a,self.b,self.c]) / mag
        return unit_norm

    def vector_angle(self,vector,degrees=False):
        dot_product = np.round(np.dot(vector,self.unit_normal()),6)
        v_product = np.round(np.dot(vector,vector),6)
        theta = np.arcsin(dot_product / np.sqrt(v_product))

        return np.degrees(theta) if degrees else theta

    def point_distance(self,coordinates): 
        x1, y1, z1 = coordinates[0], coordinates[1], coordinates[2]
        d = np.abs((self.a * x1 + self.b * y1 + self.c * z1 + self.d)) 
        e = (np.sqrt(self.a * self.a + self.b * self.b + self.c * self.c))
        return d/e

    def test_planarity(self,atoms = None):
        if atoms == None:
            devs = [self.point_distance(atom) for atom in self.atoms]
            if len(np.where(np.array(devs)>2)[0]) >= 1:
                return False
            else:
                return True
        else:
            devs = [self.point_distance(atom) for atom in atoms]
            if len(np.where(np.array(devs)>2)[0]) >= 1:
                return False
            else:
                return True

    def get_planar_basis(self):
        normal = np.array(self.a,self.b,self.c)

    def project(self,point):
        # Projects an n by 3 array of points onto itself
        point = point.reshape(1,3)
        print(point)
        normal = self.unit_normal().reshape(1,3)
        disp = (point - self.ctr).reshape(1,3)
        dist = np.sqrt(np.dot(disp,disp.T))
        theta = np.arcsin(np.dot(normal,disp.T) / (np.sqrt(np.dot(disp,disp.T)) * np.sqrt(np.dot(normal,normal.T))))
        proj = dist*sin(theta)*self.unit_normal()
        print(proj.reshape(-1))

        return proj.reshape(-1)

    def plot(self,ax=None,unit_norm=True,x_min=0,x_max=10,y_min=0,y_max=10,density=100):
        if ax is None:
            fig = plt.figure(figsize=[8,5],dpi=100)
            ax = fig.add_subplot(111,projection="3d")
        else:
            pass
        normal = self.unit_normal()
        normal = np.array([self.a,self.b,self.c])
        x = np.linspace(x_min, x_max, density)
        y = np.linspace(y_min, y_max, density)
        xx, yy = np.meshgrid(x, y)
        z = (-normal[0] * xx - normal[1] * yy - self.d) * 1. /normal[2]
        ax.plot_surface(xx, yy, z, alpha=0.2)
        if unit_norm:
            new_point = self.ctr + normal
            print(new_point)
            new_x, new_y, new_z = new_point
            ax.scatter(new_x,new_y,new_z)
            ax.scatter(normal[0]+self.d,normal[1]+self.d,normal[2]+self.d)
            ax.plot(xs=[new_x,self.ctr[0]],ys=[new_y,self.ctr[1]],zs=[new_z,self.ctr[2]],label='Unit Normal',color='black')

        return ax