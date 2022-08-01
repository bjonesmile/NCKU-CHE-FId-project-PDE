import numpy as np
import pandas as pd
import math
import plotly.graph_objects as go

def Ackley_func(x, y):
    pi = math.pi
    val = -20*np.exp(-0.2*np.sqrt(0.5*(x**2+y**2)))
    val = val - np.exp(0.5*(np.cos(2*pi*x)+np.cos(2*pi*y)))
    val = val + math.e+20
    return val
x_ary = np.linspace(-5,5,50,endpoint=True)
y_ary = np.linspace(-5,5,50,endpoint=True)
X, Y = np.meshgrid(x_ary, y_ary)
Z = Ackley_func(X, Y)
fig = go.Figure(data=[go.Surface(z=Z, x=x_ary, y=y_ary)])
fig.update_layout(title='Mt Bruno Elevation', autosize=False,
                  width=500, height=500,
                  margin=dict(l=65, r=50, b=65, t=90))

"""
z_data = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/api_docs/mt_bruno_elevation.csv')
z = z_data.values
sh_0, sh_1 = z.shape
print(z.shape)
x, y = np.linspace(0, 1, sh_0), np.linspace(0, 1, sh_1)
fig = go.Figure(data=[go.Surface(z=z, x=x, y=y)])
fig.update_layout(title='Mt Bruno Elevation', autosize=False,
                  width=500, height=500,
                  margin=dict(l=65, r=50, b=65, t=90))
"""

fig.show()