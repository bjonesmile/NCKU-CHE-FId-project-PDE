import numpy as np
from numpy.core.fromnumeric import shape

rnd = np.random.uniform(0,1)

def check_input(sol,v):
    p = sol.argsort()
    sol[:] = sol[p]
    v[:] = v[p]
    print(p)
    print(sol)
    print(v)


if __name__ == '__main__':
    a = np.array([3,37,23,350,242])
    v = np.array([1,2,3,4,5])
    check_input(a,v)
    print(a)
    print(v)

    a = np.array([4.1,1.7,3.45,5.1,2.6])
    a = np.around(a)
    print(a)
    a = np.sort(a)
    print(a)
    a = a.astype('int32')
    print(a)
    print(a.dtype)