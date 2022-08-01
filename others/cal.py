import numpy as np
import math


def cal_Ht(q,Tf):
    return math.pow(q,0.8)*(2.17e6-5.01e4*Tf+453*Tf**2-1.64*Tf**3+2.1e-3*Tf**4)

def cal_Hl(I,Tf):
    return 0.00249*(Tf-25)-0.06133

def cal_density(T):
    return -0.90797*T+0.00078116*T**2-2.367e-6*T**3+1083.25

def cal_cp(T):
    return 0.002414*T+5.9591e-6*T**2-2.9879e-8*T**3+4.4172e-11*T**4+1.498

def cal_k(T):
    return -8.19477e-5*T-1.92257e-7*T**2+2.5034e-11*T**3-7.2974e-15*T**4+0.1337746

def cal_Enthalpy(T):
    return -18.17+1.496*T-0.000147*T**2

def cal_k_viscosity(T):
    return math.exp(544.1449/(T+114.43)-2.59578)

if __name__ == "__main__":
    I = 200
    G = 0.575
    Tm = 400
    Tf = 250
    Da = 0.06
    A = 0.00196
    q = 0.6
    m = 3

    start_T = 175
    end_T = 225

    sum = np.zeros(6)
    for t in range(start_T,end_T+1):
        var = [cal_density(t), cal_cp(t), cal_k(t), cal_Enthalpy(t), cal_k_viscosity(t)]
        sum[0] += cal_density(t)
        sum[1] += cal_cp(t)*1000
        sum[2] += cal_k(t)
        sum[3] += cal_Enthalpy(t)
        sum[4] += cal_k_viscosity(t)/1000000
        Re = Da*m/(cal_k_viscosity(t)/1000000)
        Pr = (cal_k_viscosity(t)/1000000)*cal_density(t)*cal_cp(t)*1000/cal_k(t)
        hp = 0.023*math.pow(Re,0.8)*math.pow(Pr,0.4)*cal_k(t)/Da
        sum[5] += hp
        print(t,': ',var)

    

    print(f"avg({start_T} to {end_T})")
    print(sum/(end_T-start_T))

        


