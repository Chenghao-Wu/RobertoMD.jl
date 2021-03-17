from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
data=pd.read_csv("data",sep="\s+",header=None)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = data[0]
y = data[1]
z = data[2]
c = data[3]
ax.view_init(45,60)
img = ax.scatter(x, y, z, c=c, cmap=plt.hot())
fig.colorbar(img)
fig.savefig("plot.png")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x,y = np.meshgrid(data[0], data[1])
z = data[2]
c = data[3]

ax.view_init(45,60)

# here we create the surface plot, but pass V through a colormap
# to create a different color for each patch
ax.plot_surface(x, y, z, facecolors=cm.Oranges(c))
fig.savefig("plot2.png")
