import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import dill as pickle


def plotfun(plotdata):
    ax = plt.gca()
    ax.plot_surface(plotdata['xs'], plotdata['ys'], plotdata['zs'])


xs = np.linspace(0, 1, 100)
ys = np.linspace(0, 4*np.pi, 100)
zs = np.zeros((100, 100))
for j, phase in enumerate(ys):
    zs[j, :] = xs*np.sin(phase)
xs, ys = np.meshgrid(xs, ys)

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
plotdata = {'xs': xs, 'ys': ys, 'zs': zs, 'plotfun': plotfun}
plotdata['plotfun'](plotdata)
with open("test.pickle", 'wb') as file:
    pickle.dump(plotdata, file)