import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    x = np.linspace(0,10,50)
    y = x**2
    X , Y = np.meshgrid(x,y)
    z = np.sin(X) + np.cos(Y)

    fig = plt.figure(figsize=plt.figaspect(2.))
    ax = fig.add_subplot(4,1,(3,4), projection='3d')
    
    surface = ax.plot_surface(X, Y, z, color='blue', rstride=1, cstride=1, cmap='coolwarm')
    plt.colorbar(surface,shrink=1.0,aspect=20)
    
    ax1 = fig.add_subplot(4,1,1)
    ax1.plot(y,z)

    ax2 = fig.add_subplot(4,1,2)
    ax2.plot(x,y)

    fig.tight_layout()  

    plt.show()