import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Ian if you read this, it produces the incorrect plane.
# For the life of me I cannot figure out why.
# The Mathematica notebook (.nb) produces the figures used in my computational report, I hope that is ok.
# Please let me know if you see why this produces the wrong plane. 

#Define coefficients
a = 0.0390336
b = 0.363382
c = 0.0780671


#Define points
#A = np.array([1,2,3])
A = (1,2,3)
B = (-3,2,5)
C = (np.pi, np.exp(1), -(2)**(1/2))


def f(x,y):
    z = (1/c) - (a/c)*x - (b/c)*y
    return z




xx = np.arange(-4,4,.01)
yy = np.arange(0,4,.01)
X, Y = np.meshgrid(xx, yy)

Z = f(X, Y)


if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.title('Plot of plane containing 3 desired points')
    ax.plot_surface(X, Y, Z)
    ax.scatter(A, B, C)
    plt.show()
    #ax.legend()
    plt.savefig('planeGraph.png')
    #os.system('cp planeGraph.png /mnt/c/Users/nicho/Pictures/WSL_Pictures/AM_213a')

