from unittest import result
import pandas as pd
import numpy as np
import os
import sys
import pickle
import itertools
import random
import copy
import matplotlib.pyplot as plt
from datetime import datetime

import PTC as model
import mytimer as timer

try:
    loadFileName = sys.argv[1]
    with open(loadFileName, 'rb') as f:
        loadInfo= pickle.load(f)
except IndexError:
    print("must input rough result pickle file")

Nz = loadInfo['Nz']
Uncertain_parameter_name = []
for var_name in loadInfo['StartDirect']:
    Uncertain_parameter_name.append(var_name)
StartDirect = loadInfo['StartDirect']
vertex = loadInfo['vertex']
for key in vertex.keys():
    vertex[key] = list(vertex[key])

print(f"load Info from  {loadFileName} suc")
print('Nz',Nz)
print('Uncertain_parameter_name',Uncertain_parameter_name)
print('StartDirect',StartDirect)
print('Vertex',vertex)

segNum = len(Uncertain_parameter_name)
GAMSfile = "PTC_modelinstance_I_nz"+str(Nz)+"_"+str(segNum)+"seg.gms"
gdxInGAMSFILE = "PTC_modelinstance_I_nz"+str(Nz)+"_"+str(segNum)+"seg_gdxin.gms"

solver = model.PTC(dt=30,dx=100,nfe=5)
solver.set_Nz(Nz)
solver.set_Seg(segNum)
solver.set_uncertain_para(Uncertain_parameter_name)
solver.setStartDirect(StartDirect)
solver.set_GAMSfile(GAMSfile)
solver.set_upper_bound("T_o")
solver.set_lower_bound("T_o")
solver.set_gdxInGAMSFILE(gdxInGAMSFILE)
solver.set_GAMSModelInstance()
solver.is_show_solve_result = False
print(solver.uncertainParaList)

if solver.total_sim_time is None:
    raise AttributeError("modle time len of solver {} is not exist".format(solver.model_name))

def gen_deviated_list(time_len):
    deviated_list = []
    p = [0.01, 0.015, 0.02, 0.025]
    for ip in p:
        t= int(time_len*ip)
        if t <= 5 :
            break
        else:
            deviated_list.append(t)
    deviated_list.extend([i for i in range(5,0,-1)])
    deviated_list = sorted(deviated_list,reverse=True)

    return deviated_list

deviated_t = gen_deviated_list(solver.total_sim_time)
deviated_t = [1, 3, 5, 10]
#deviated_t = sorted(deviated_t,reverse=True)
print("List deviated_t:",deviated_t)

def shake(vertex, up_list,k):
    new_vertex = {}
    for up in up_list:
        pool = []
        dt = deviated_t[k]
        for vt in vertex[up] :
            pool.append([vt-dt, vt, vt+dt])
        up_vertex = []
        for item in pool:
            t = random.choice(item)
            up_vertex.append(t)
        new_vertex[up] = up_vertex
    return new_vertex

def LocalSearch(vertex, up_list, FirstMove=True):
    # bulid
    pool = {}
    vertex_pool = {}
    for up in up_list:
        pool[up] = []
        for vt in vertex[up] :
            pool[up].append([vt-1, vt, vt+1])
        vertex_pool[up] = itertools.product(*pool[up])
        # if transform iterator to list will cause problem to product
        #print(list(vertex_pool[up]))
    
    print("\nlocalsearch:",vertex)
    # best
    #min_f = 1000
    #min_v = None
    # first
    min_f = solver.solve(vertex)
    min_v = vertex
    for item in itertools.product(*(vertex_pool.values())):
        scen = {}
        for i in range(len(up_list)):
            scen[up_list[i]] = item[i]
        #print(scen)
        f = solver.solve(scen)
        if f <  min_f:
            min_f = f
            min_v = copy.deepcopy(scen)
            print("local opt:",min_f,min_v)
            if FirstMove :
                return min_v, min_f
    return min_v, min_f

def NeighbourhoodChange(obj_min, sol_min, obj, sol, k):
    if obj < obj_min:
        return obj, sol, 0
    else:
        return obj_min, sol_min, k+1

def variable_neighborhood_search(vertex, kmax, uplist):
    k = 0
    sol = vertex
    best_f = solver.solve(vertex)
    last_f = best_f
    while True:
        print("\rk: {} vertex: {}.".format(k, sol), end="")
        sys.stdout.flush()
        sol_t = shake(sol,uplist, k)
        sol_tt, f_tt = LocalSearch(sol_t,uplist,FirstMove=True)
        best_f, sol, k = NeighbourhoodChange(best_f, sol, f_tt, sol_tt, k)
        if k == kmax and last_f == best_f:
            break
        if k == kmax:
            k = kmax-1
        last_f = best_f

    return best_f, sol

def simplfy(vertex, up_list):
    org_f = solver.solve(vertex)
    err_std = 0.00001*org_f
    vertex_R = copy.deepcopy(vertex)

    for key in vertex_R.keys():
        #print(vertex_R[key])
        try:
            vertex_R[key].remove('t')
        except:
            pass
    #print("after remove 't':",vertex_R)

    k = sum([len(vertex[up]) for up in up_list])
    while k != 0 :
        refined_vertex = copy.deepcopy(vertex_R)
        for _ in range(k):
            if sum([len(refined_vertex[up]) for up in up_list]) == 0:
                break
            max_ver_t = max([max(refined_vertex[up], default=0) for up in up_list])
            key = next(up for up in up_list if max_ver_t in refined_vertex[up])
            #print(f"remove ver_t {key}: {max_ver_t}")
            #print(refined_vertex[key])
            refined_vertex[key].remove(max_ver_t) 

        f = solver.solve(refined_vertex)
        if abs(org_f-f) <= err_std:
            #vertex_R =  refined_vertex
            #k = sum([len(vertex_R[up]) for up in up_list])
            return refined_vertex
        else:
            k -= 1

    return vertex_R


if __name__ == '__main__':
    kmax = len(deviated_t)
    
    up_list = Uncertain_parameter_name
    loaded_vertex = vertex

    sim_result = simplfy(loaded_vertex, up_list)
    print("original vertex: ",loaded_vertex)
    print("simplfy vertex: ",sim_result)
    print("")
    input()
    VNS_timer = timer.MyTimer()
    VNS_result = variable_neighborhood_search(sim_result,kmax,up_list)
    usingTime = VNS_timer.getTime(kind='real')
    print(f"\nNz: {Nz}\tStatrDirect: {StartDirect}\noriginal vertex: {sim_result}\nVNS result: {VNS_result}")
    
    saveFileName = os.getcwd()+'\\VNS-Nz'+str(Nz)+'-'+str(len(Uncertain_parameter_name))+'seg-result'+str(datetime.now().strftime('%Y-%m-%d'))+'.pkl'
    with open(saveFileName, 'wb') as f:
        saveVertexInfo = {}
        saveVertexInfo['Nz'] = Nz
        saveVertexInfo['StartDirect'] = StartDirect
        saveVertexInfo['vertex'] = VNS_result[1]
        saveVertexInfo['FId'] = VNS_result[0]
        saveVertexInfo['using-time'] = usingTime
        pickle.dump(saveVertexInfo, f)
        print(f"save file {saveFileName} finished.")

    exit()