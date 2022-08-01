import numpy as np
from multiprocessing import Array

SIM_TIME = 350

nz= 5
org_vertex = (1,35,100,150,250,'t')
min_vertex = Array('i',np.zeros(nz*2,dtype=np.int32))

for i in range(len(min_vertex)):
    print(i,min_vertex[i])

for i in range(len(org_vertex)):
    if org_vertex[i] != 't': 
        min_vertex[i] = org_vertex[i]
    else:
        min_vertex[i] = SIM_TIME+1
for i in range(len(org_vertex),len(min_vertex)):
    min_vertex[i] = -1

for i in range(len(min_vertex)):
    print(i,min_vertex[i])