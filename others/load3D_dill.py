import dill as pickle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


plt.figure().add_subplot(111, projection='3d'),
plotdata = pickle.load(open('test.pickle', 'rb'))
plotdata['plotfun'](plotdata)
plt.show()