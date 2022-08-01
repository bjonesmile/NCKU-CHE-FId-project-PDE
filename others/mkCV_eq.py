import numpy as np
import sys

# Single
# time_H: 800
# time set: m
# CVeq-Name: ['qout']

# alcoholicCSTR
# time_H: 240
# time set: m
# CVeq-Name: ['Fin', 'Fout']

# CHS
# time_H: 350
# time set: T
# CVeq-Name: ['BC_L', 'BC_R']

# CHS
# time_H: 240
# time set: t
# CVeq-Name: ['m']

model_name = "ptc"
nz = int(sys.argv[1])
time_H = 240
time_set_name = 't'
eq_name_list = ['m']
ary = np.linspace(0,time_H,nz+1,dtype=np.int32)
ary = np.delete(ary,0)

print(ary)

file_name = model_name+"_nz"+str(nz)+'.txt'
with open(file_name,'w') as out_file:
    for eq_name in eq_name_list:
        for t in range(len(ary)):
            out_file.write(f"\t{eq_name}{t+1}({time_set_name})\n")
    out_file.write("\n")
    for eq_name in eq_name_list:
        out_file.write(f"\t{eq_name}1({time_set_name})$(ord({time_set_name}) lt {ary[0]}).. {eq_name}({time_set_name}) =e= {eq_name}({time_set_name}+1);\n" )
        for t in range(len(ary)-1) :
            out_file.write(f"\t{eq_name}{t+2}({time_set_name})$((ord({time_set_name}) gt {ary[t]}) and (ord({time_set_name}) lt {ary[t+1]}))..\
 {eq_name}({time_set_name}) =e= {eq_name}({time_set_name}+1);\n" )