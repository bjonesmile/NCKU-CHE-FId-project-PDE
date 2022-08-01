import random 
import math
import copy
import numpy as np

import ConvectionDiffusionwithHeatSink as problem
import mytimer as timer

proble_Nz = 2
Nz = proble_Nz*2
MIN_POINT_DISTANCE = 5
MIN_RANGE_TO_NEXT_ACT = 25
SIM_TIME = 350
AVAILABLE_CHOICE_NUMBER = 0
GAMMA = 1.0
try_limit = 10
FId_Precision = None
change_temp = False
# uncertain Vx_d, alpha_d, T_out_d
Uncertain_parameter_name = 'T_out_d'
if Uncertain_parameter_name == 'T_out_d':
    set_model_file = 'CDHS_modelinstance_T_nz'+str(proble_Nz)+'.gms'

V_d = 1

solver = problem.Temp_distribute1D(0.001,0.05)
solver.set_Nz(proble_Nz)
solver.setStartDirect(V_d)
solver.set_GAMSfile(set_model_file)
solver.set_temp_base_value(0.4)
solver.set_IC_X(0.3,3,1,ctype='eq',showfig=False)
solver.set_upper_bound()
solver.set_lower_bound()
solver.set_GAMSModelInstance(uncertain_para=[Uncertain_parameter_name])

class BatAlgorithm():
    def __init__(self, D, NP, N_Gen, A, r, Qmin, Qmax, Lb, Ub, function):
        self.D = D  #dimension
        self.NP = NP  #population size 
        self.N_Gen = N_Gen  #generations
        self.A = A*np.ones(self.NP,dtype=np.float32) #loudness
        self.r0 = r*np.ones(self.NP,dtype=np.float32)
        self.r = r*np.ones(self.NP,dtype=np.float32)  #pulse rate
        self.Qmin = Qmin  #frequency min
        self.Qmax = Qmax  #frequency max

        self.f_min = 0.0  #minimum fitness
        
        self.Lb = np.copy(Lb)  #lower bound
        self.Ub = np.copy(Ub)  #upper bound
        self.Q = np.zeros(self.NP)  #frequency

        #self.v = [[0 for i in range(self.D)] for j in range(self.NP)]  #velocity
        self.v = np.zeros((self.NP,self.D))
        self.Sol = np.zeros((self.NP,self.D))  #population of solutions
        self.Fitness = np.zeros(self.NP)  #fitness
        self.best = np.zeros(self.D) #best solution
        self.Fun = function

        self.alpha = 0.99
        self.gamma = 0.9

    def best_bat(self):
        j = np.argmin(self.Fitness)
        self.best[:] = self.Sol[j][:]
        self.f_min = self.Fitness[j]

    def init_bat(self):
        for i in range(self.NP):
            self.Q[i] = 0
            for j in range(self.D):
                rnd = np.random.uniform(0, 1)
                self.v[i][j] = 0.0
                self.Sol[i][j] = self.Lb[j] + (self.Ub[j] - self.Lb[j]) * rnd
            self.Fitness[i] = self.Fun(self.D, self.Sol[i])
        self.best_bat()

    def simplebounds(self, val):
        for i in range(self.D):
            if val[i] < self.Lb[i]:
                val[i] = self.Lb[i]
            if val[i] > self.Ub[i]:
                val[i] = self.Ub[i]

    def move_bat(self):
        S = [[0.0 for i in range(self.D)] for j in range(self.NP)]

        self.init_bat()

        for t in range(self.N_Gen):
            temp_Sol = np.zeros(self.D)
            A_avg = np.average(self.A)
            for i in range(self.NP):
                beta = np.random.uniform(0, 1)
                self.Q[i] = self.Qmin + (self.Qmax - self.Qmin) * beta
                self.v[i] = self.v[i] + (self.Sol[i]-self.best) * self.Q[i]
                temp_Sol = self.Sol[i] + self.v[i]
                    
                self.simplebounds(temp_Sol)

                rnd1 = np.random.random_sample()

                if rnd1 > self.r[i]:
                    """
                    for j in range(self.D):
                        temp_Sol[j] = self.best[j] + 0.001 * random.gauss(0, 1)
                    """
                    for j in range(self.D):
                        eplison = np.random.uniform(-1,1)
                        temp_Sol[j] = self.best[j]+eplison*A_avg
                    self.simplebounds(temp_Sol)

                        
                Fnew = self.Fun(self.D, temp_Sol)

                rnd2 = np.random.random_sample()

                if (Fnew <= self.Fitness[i]) and (rnd2 < self.A[i]):
                    self.A[i] = self.alpha*self.A[i]
                    self.r[i] = self.r0[i]*(1-math.exp(-self.gamma*t))
                self.Sol[i][:] = temp_Sol[:]
                self.Fitness[i] = Fnew
            
            for i in range(self.NP):
                if self.Fitness[i] <= self.f_min:
                    for j in range(self.D):
                        self.best[j] = self.Sol[i][j]
                    self.f_min = self.Fitness[i]

        print(self.f_min)
        print(self.best)

def func(D,input):
    if len(input)!=D:
        raise IndexError()
    FId = solver.solve(input)
    print("vertex:",input)
    return FId

if __name__ == '__main__':
    Ub = np.array([350,350,350])
    Lb = np.array([0,0,0])
    
    timer1 = timer.MyTimer()
    Algorithm = BatAlgorithm(3, 40, 200, 0.95, 0.9, 0.0, 1.0, Lb, Ub, func)
    Algorithm.move_bat()
    t = timer1.getTime(kind='real')