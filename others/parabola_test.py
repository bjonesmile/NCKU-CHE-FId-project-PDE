import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    f = 100
    v1 = 240/2
    v2 = 225
    x = np.linspace(0,240,240)
    y = -1/(4*f)*(x-v1)**2+v2
    plt.plot(x,y)
    plt.show()
