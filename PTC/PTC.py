from msilib.schema import Directory
from multiprocessing import managers
from operator import le
from unicodedata import name
from cfgmcc import ptrTocfgHandle
from gams import *
import sys
import os
import pickle
import numpy as np
import math
import copy
from datetime import datetime

import mytimer as timer

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from numpy.core.numeric import NaN

class PTC():
    def __init__(self,dt=12,dx=100,nfe=10):
        self.total_sim_time = 240

        self.ws = GamsWorkspace(working_directory = os.getcwd())
        self.model_name = "PTC"

        self.is_show_solve_result = False

        self.dt = dt
        self.dx = dx
        self.nfe =nfe
        self.ocfe = 5
        self.drawing = 1

        self.t_grid = np.linspace(0,self.total_sim_time,self.total_sim_time,endpoint=True)
        self.roots = None
        if self.ocfe == 5:
            roots = [0.0, 0.1127, 0.5, 0.8873, 1.0]
        x_grid = [x*dx for x in roots]
        for i in range(1,self.nfe):
            x_grid.extend([x*dx+dx*i for x in roots[1:]])
        self.x_grid = np.array(x_grid)

        # upper-bound
        self.upper_bound = {"T_a": None,"T_o": None,"T_e": None}
        self.lower_bound = {"T_a": None,"T_o": None,"T_e": None}

        self.Nz = 1
        self.Seg = 1
        self.uncertainParaList = ["I"]
        self.startDirect = {"I": 1}
        self.GAMSModelInstance = None
        self.cp = None
        self.job = None
        self.cp_file = None

        self.gdxInGAMSFILE = None
        self.gdxInDir = "ptc_nz"+str(self.Nz)+"_seg"+str(self.Seg)

    def set_Nz(self,nz):
        self.Nz = nz
        self.gdxInDir = "ptc_nz"+str(self.Nz)+"_seg"+str(self.Seg)

    def set_Seg(self,seg):
        self.Seg = seg
        self.gdxInDir = "ptc_nz"+str(self.Nz)+"_seg"+str(self.Seg)

    def setStartDirect(self,d):
        self.startDirect = d
        for p in self.uncertainParaList:
            if p not in d:
                raise AttributeError("Uncertain Parameter ", p," not in startdirect dict.")
            p_direct = self.startDirect[p]
            if not (p_direct == 1 or p_direct == -1):
                raise ValueError("Uncertain Parameter ", p, " startdirect is ", p_direct, ", it should be +1 or -1.")
        print("current Uncertain Parameter start direction",self.startDirect)

    def set_uncertain_para(self,uncertain_para):
        if isinstance(uncertain_para,str):
            self.uncertainParaList = [uncertain_para]
        elif isinstance(uncertain_para,list):
            self.uncertainParaList = uncertain_para
        else:
            raise TypeError("uncertain paramaters setting must be str or list type.")

    def set_GAMSfile(self,file=None):
        if file is None:
            file_name = self.model_name+"_modelinstance"
            for p in self.uncertainParaList:
                file_name += '_'+p
            file_name += '_nz'+str(self.Nz)+'_'+str(self.Seg)+'seg_gdxin.gms'
        else:
            file_name = file
        print("GAMSgdxInfile:", file_name)
        if os.path.isfile('./'+file_name):
            self.jobGAMSfile = file_name
        else:
            raise FileNotFoundError("GAMSgdxInfile file :"+file_name+"not found")

    def set_gdxInGAMSFILE(self,file=None):
        if file is None:
            file_name = self.model_name+"_modelinstance"
            for p in self.uncertainParaList:
                file_name += '_'+p
            file_name += '_nz'+str(self.Nz)+'.gms'
        else:
            file_name = file
        print("GAMSfile:", file_name)
        if os.path.isfile('./'+file_name):
            self.gdxInGAMSFILE = file_name
        else:
            raise FileNotFoundError("GAMSModelInstance file :"+file_name+"not found")

    def set_gdxInDir(self,dir=None):
        if dir is None :
            dir_name = "ptc_nz"+str(self.Nz)+"_seg"+str(self.Seg)
        else:
            if os.path.isdir(dir_name):
                self.gdxInDir = dir_name
            else:
                raise FileNotFoundError("GAMS gdxIn dir :"+dir_name+"not found")
        print("set gdxInDir:", self.gdxInDir)

    def set_GAMSModelInstance(self,solver="conopt"):
        self.cp = self.ws.add_checkpoint()
        self.job = self.ws.add_job_from_file(self.jobGAMSfile)
        self.job.run(checkpoint=self.cp)
        mi = self.cp.add_modelinstance()

        #modifiers = [GamsModifier(load_value) for load_value in load_values]

        opt = self.ws.add_options()
        opt.all_model_types = solver
        opt.solprint = 1
        print("fdopt:\t\t",opt.fdopt)
        print("fddelta:\t",opt.fddelta)

        modifiers = []

        to = mi.sync_db.add_variable("T_o", 3 ,VarType.Free)
        up_limit = mi.sync_db.add_parameter("up_limit", 3, "upper limit plane constraint")
        modifiers.append(GamsModifier(to, UpdateAction.Upper, up_limit))
        
        lo_limit = mi.sync_db.add_parameter("lo_limit", 3, "lower limit plane constraint")
        modifiers.append(GamsModifier(to, UpdateAction.Lower, lo_limit))

        for p in self.uncertainParaList:
            vertex = mi.sync_db.add_parameter(p+'_vertex', 1, p+"_d direction vertex")
            modifiers.append(GamsModifier(vertex))

        mi.instantiate(self.model_name+" using nlp maximizing obj",modifiers,opt)
        for t in range(1,self.total_sim_time+1):
            i = 0
            for x in range(1,self.nfe+1):
                for j in range(1,self.ocfe+1):
                    up_limit.add_record((str(t),str(x),str(j))).value = self.upper_bound["T_o"][t-1][i]
                    lo_limit.add_record((str(t),str(x),str(j))).value = self.lower_bound["T_o"][t-1][i]
                    if j != 5:
                        i += 1
        
        self.GAMSModelInstance = mi

    def set_upper_bound(self,sv):
        f = 80
        v1 = self.total_sim_time/2
        v2 = 225

        X , Y = np.meshgrid(self.x_grid,self.t_grid)
        Z = -1/(4*f)*(Y-v1)**2+v2
        self.upper_bound[sv] = Z+100

    def set_lower_bound(self,sv):
        f = 80
        v1 = self.total_sim_time/2
        v2 = 225

        X , Y = np.meshgrid(self.x_grid,self.t_grid)
        Z = -1/(4*f)*(Y-v1)**2+v2
        self.lower_bound[sv] = Z-20

    def update_vertex(self,p_name,shift_list):
        parameter = self.GAMSModelInstance.sync_db.get_parameter(p_name+'_vertex')
        parameter.clear()

        v_len_V = len(shift_list)
        d_cur_V = self.startDirect[p_name]
        shift_t = 0
        for rec in self.job.out_db.get_parameter(p_name+'_vertex'):
            #print(rec.keys)
            key = "".join([char for char in rec.keys[0] if char.isdigit()])
            if shift_t < v_len_V:
                if shift_list[shift_t] == 't':
                    shift_t = 10000
                elif int(key) >=shift_list[shift_t] :
                    d_cur_V = -d_cur_V
                    shift_t += 1
            parameter.add_record(key).value = d_cur_V

    def solve(self,vertex):
        mi = self.GAMSModelInstance

        for p in self.uncertainParaList:
            self.update_vertex(p,vertex[p])
        
        mi.solve()
        while mi.model_status != 2 or mi.solver_status != 1:
            #print("re-solve vertex=" + str(vertex) + ":")
            rf = self.solve_by_gdxIn(vertex)
            if rf != 0.0:
                return rf

        if self.is_show_solve_result:
            print("Scenario vertex=" + str(vertex) + ":")
            print("  Modelstatus: " + str(mi.model_status))
            print("  Solvestatus: " + str(mi.solver_status))
            print("  Obj: " + str(mi.sync_db.get_variable("obj")[()].level))
        result = mi.sync_db.get_variable("obj")[()].level

        return result

    def solve_by_gdxIn(self,vertex):
        w_dir = os.path.join(".",self.gdxInDir)
        ws = GamsWorkspace(w_dir)

        db = ws.add_database(database_name="Db_vertex")
        for p in self.uncertainParaList:
            up = db.add_parameter(p+'_vertex',1)
            v_len_V = len(vertex[p])
            d_cur_V = self.startDirect[p]
            shift_t = 0
            for i in range(1,self.total_sim_time+1):
                if shift_t < v_len_V:
                    if vertex[p][shift_t] == 't':
                        shift_t = 10000
                    elif i >= vertex[p][shift_t] :
                        d_cur_V = -d_cur_V
                        shift_t += 1
                up.add_record(str(i)).value = d_cur_V
        
        job = ws.add_job_from_file(self.gdxInGAMSFILE)
        db.export(db.name)

        opt = ws.add_options()
        opt.defines["gdxincname"] = db.name
        opt.all_model_types = "conopt4"
        job.run(opt)
        del opt

        if self.is_show_solve_result:
            print("Scenario vertex=" + str(vertex) + ":")
            print("  Modelstatus: " + str(job.out_db["ms"].find_record().value))
            print("  Solvestatus: " + str(job.out_db["ss"].find_record().value))
            print("  Obj: " + str(job.out_db.get_variable("obj")[()].level))
        result = job.out_db.get_variable("obj")[()].level

        return result

    def solve_outDB(self,vertex):
        mi = self.GAMSModelInstance

        for p in self.uncertainParaList:
            self.update_vertex(p,vertex[p])
        
        mi.solve()
        print("Scenario vertex=" + str(vertex) + ":")
        print("  Modelstatus: " + str(mi.model_status))
        print("  Solvestatus: " + str(mi.solver_status))
        print("  Obj: " + str(mi.sync_db.get_variable("obj")[()].level))
        FId = mi.sync_db.get_variable("obj")[()].level
        
        export_gdx_name = os.getcwd()+'\\Nz'+str(self.Nz)+'-'+str(len(self.uncertainParaList))+'Seg-result.gdx'
        mi.sync_db.export(export_gdx_name)
        print("save gdx:",export_gdx_name)
        return export_gdx_name
        
    def plot_Db(self, gdxfile, sv_list=["T_a","T_o","T_e"]):
        ws = GamsWorkspace(working_directory = os.getcwd())
        db = ws.add_database_from_gdx(gdxfile)
        
        dirPath = os.getcwd()+'\\Nz'+str(self.Nz)+'-result'
        if os.path.isdir(dirPath) == False:
            os.mkdir(dirPath)
        
        dx = self.dx
        roots = None
        if self.ocfe == 5:
            roots = [0.0, 0.1127, 0.5, 0.8873, 1.0]
        x_grid = [x*dx for x in roots]
        for i in range(1,self.nfe):
            x_grid.extend([x*dx+dx*i for x in roots[1:]])

        x_grid = np.array(x_grid)
        t_grid = np.linspace(0,self.total_sim_time,self.total_sim_time)
        X, Y = np.meshgrid(x_grid,t_grid)

        m_data = []
        for rec in db["m"]:
            m_data.append(rec.level)
        x_len_num = len(x_grid)
        t_len_num = len(t_grid)
        I_data = np.empty((t_len_num,x_len_num))
        i, j = 0, 0
        for rec in db["I"]:
            if(rec.key(1) != '1' and rec.key(2) == '1'):
                continue
            I_data[i][j] = rec.level
            j += 1
            if j >= I_data.shape[1] :
                i += 1
                j = 0

        t_len = self.total_sim_time
        m_ub = 4.0
        m_lb = 1.5
        fig_m = plt.figure()
        ax_m = fig_m.add_subplot()
        ax_m.plot(t_grid,m_data,color='steelblue',linestyle='-',label=r'$\dot m$')
        ax_m.hlines(m_ub,xmin=0,xmax=t_len,color='black',linestyle='--')
        ax_m.hlines(m_lb,xmin=0,xmax=t_len,color='black',linestyle='--')
        ax_m.set_xlim(0, 240)
        ax_m.yaxis.set_ticks(np.linspace(1.5,4.0,num=6,endpoint=True))
        ax_m.xaxis.set_ticks(np.linspace(0,240,num=7,endpoint=True))
        ax_m.set_xlabel("time(min)")
        ax_m.set_ylabel("flow rate (kg/s)")
        ax_m.grid(color='y', linestyle='--', linewidth=1, alpha=0.3)
        ax_m.legend()

        plt.tight_layout()
        figName ='Nz'+str(self.Nz)+'-'+str(self.Seg)+'Seg-CV.png'
        figPath = os.path.join(dirPath, figName)
        fig_m.savefig(figPath,dpi=200)
        print("save fig:",figName)

        fig_I = plt.figure()
        ax_I = fig_I.add_subplot()
        c = ax_I.pcolor(X,Y,I_data,cmap='seismic',shading='auto')
        cbar_I = fig_I.colorbar(c,ax=ax_I,shrink=1.0,aspect=20)
        cbar_I.ax.get_yaxis().labelpad = 15
        cbar_I.ax.set_ylabel("uncertain $I(W/m^{2})$", fontsize=12, rotation=270)
        ax_I.set_xlabel("collector length(m)",fontsize=12)
        ax_I.set_ylabel("time(min)",fontsize=12)

        figName ='Nz'+str(self.Nz)+'-'+str(self.Seg)+'Seg-UP.png'
        figPath = os.path.join(dirPath, figName)
        fig_I.savefig(figPath,dpi=200)
        print("save fig:",figName)

        data = {}
        for sv in sv_list:
            data[sv] = np.empty((t_len_num,x_len_num))

        def getDataSV(data, db, sv):
            upper_activated = []
            lower_activated = []

            i, j = 0, 0
            err_limit = 0.00001
            for rec in db[sv] :
                if(rec.key(1) != '1' and rec.key(2) == '1'):
                    continue
                data[i][j] = rec.level
                if abs(rec.level-rec.lower) <= err_limit and rec.key(0) != '1':
                    lower_activated.append((i,j))
                elif abs(rec.level-rec.upper) <= err_limit and rec.key(0) != '1':
                    upper_activated.append((i,j))
                j += 1
                if j >= data.shape[1] :
                    i += 1
                    j = 0
                
            return upper_activated, lower_activated
        
        """db_To_lb = self.ws.add_database()
        db_To_ub = self.ws.add_database()
        To_lb = db_To_lb.add_parameter("T_o_lb", 3, "T_o lower bound")
        To_ub = db_To_ub.add_parameter("T_o_ub", 3, "T_o upper bound")
        for rec in db["T_o"] :
            To_lb.add_record(rec.get_keys()).value = rec.lower
            To_ub.add_record(rec.get_keys()).value = rec.upper
        db_To_lb.export(os.getcwd()+'\\T_o_lb.gdx')
        db_To_ub.export(os.getcwd()+'\\T_o_ub.gdx')"""

        upper_activated = {}
        lower_activated = {}
        for sv in sv_list:
            upper_activated[sv], lower_activated[sv] = getDataSV(data[sv], db, sv)

        #fig.suptitle("Nz"+str(self.Nz)+" FId:"+str(round(FId,4)))
        def plotData(data, svName, lower_bound, upper_bound,lowerAct, upperAct):
            print(svName)
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(projection='3d')
            surface = ax.plot_surface(X, Y, data, color='blue', rstride=1, cstride=1, cmap='coolwarm')
            plt.colorbar(surface,shrink=1.0,aspect=20)
            if svName == "T_o":
                ax.plot_surface(X, Y, upper_bound, color='black', rstride=1, cstride=1, alpha=0.3)
                ax.plot_surface(X, Y, lower_bound, color='black', rstride=1, cstride=1, alpha=0.3)

            #ax.set_title(svName, fontsize=12)
            if svName == "T_o":
                z_label_text = r'$T_{o}(°C)$'
            elif svName == "T_a":
                z_label_text = r'$T_{a}(°C)$'
            elif svName == "T_e":
                z_label_text = r'$T_{e}(°C)$'
            
            ax.zaxis.set_rotate_label(False)  # disable automatic rotation
            ax.set_zlabel(z_label_text,fontsize=12,rotation=90)
            ax.set_xlabel("collector length(m)",fontsize=12)
            ax.set_ylabel("time(min)",fontsize=12)

            upper_activated_num = len(upperAct)
            lower_activated_num = len(lowerAct)
            for act in upperAct:
                ax.scatter(x_grid[act[1]],t_grid[act[0]],data[act[0]][act[1]],color="r",marker='2',s=500)
                print(f"upper_activated point x:{act[1]} {x_grid[act[1]]},t:{act[0]} {t_grid[act[0]]} value{data[act[0]][act[1]]}")
                
            for act in lowerAct:
                ax.scatter(x_grid[act[1]],t_grid[act[0]],data[act[0]][act[1]],color="r",marker='1',s=500)
                print(f"lower_activated point x:{act[1]} {x_grid[act[1]]},t:{act[0]} {t_grid[act[0]]} value{data[act[0]][act[1]]}")
            print("upper activate:",upper_activated_num)
            print("lower activate:",lower_activated_num)
            #set view point
            ax.azim = -145
            ax.dist = 10
            ax.elev = 30

            fig.tight_layout()  
            return fig
        
        figs = {}
        for sv in sv_list:
            figs[sv] = plotData(data[sv],sv,self.lower_bound[sv],self.upper_bound[sv],lower_activated[sv],upper_activated[sv])
        
        for svName, fig in figs.items():
            figName ='Nz'+str(self.Nz)+'-'+str(len(self.uncertainParaList))+'Seg-'+str(svName)+'.png'
            figPath = os.path.join(dirPath, figName)
            fig.savefig(figPath,dpi=200)
            print("save fig:",figName)

        plt.tight_layout()
        plt.show()
        plt.close('all')

if __name__ == '__main__' :
    test_Nz = 10
    test_Seg = 5
    GAMSfile = "PTC_modelinstance_I_nz"+str(test_Nz)+"_"+str(test_Seg)+"seg.gms"
    gdxInGAMSFILE = "PTC_modelinstance_I_nz"+str(test_Nz)+"_"+str(test_Seg)+"seg_gdxin.gms"
    loadInfo = None

    if test_Seg == 5:
        #StartDirect = {'I_x1': +1,'I_x2': +1,'I_x3': +1,'I_x4': +1,'I_x5': +1}
        StartDirect = {'I_x1': -1,'I_x2': -1,'I_x3': -1,'I_x4': -1,'I_x5': -1}
        Uncertain_parameter_name = ['I_x1','I_x2','I_x3','I_x4','I_x5']
    elif test_Seg == 3:
        #StartDirect = {'I_x12': +1,'I_x3': +1,'I_x45': +1}
        StartDirect = {'I_x12': -1,'I_x3': -1,'I_x45': -1}
        Uncertain_parameter_name = ['I_x12','I_x3','I_x45']
    elif test_Seg == 1:
        #StartDirect = {'I': +1}
        StartDirect = {'I': -1}
        Uncertain_parameter_name = ['I']

    test = PTC(dt=30,dx=100,nfe=5)
    test.set_Nz(test_Nz)
    test.set_Seg(test_Seg)
    test.set_uncertain_para(Uncertain_parameter_name)
    test.setStartDirect(StartDirect)
    test.set_GAMSfile(GAMSfile)
    
    test.set_upper_bound("T_o")
    test.set_lower_bound("T_o")

    test.set_GAMSModelInstance()
    test.set_gdxInGAMSFILE(gdxInGAMSFILE)
    test.is_show_solve_result = True
    #test.showInitSetting()

    #vertex = {'I': [110]}
    #f = test.solve(vertex)
    #vertex = {'I_x12': [109], 'I_x3': [122], 'I_x45': [126]}
    #f = test.solve(vertex)
    vertex =  {'I_x1': [], 'I_x2': [], 'I_x3': [], 'I_x4': [], 'I_x5': []}
    #f = test.solve(vertex)
    dbName = test.solve_outDB(vertex)
    test.plot_Db(gdxfile=dbName)
    exit()
    
    timer1 = timer.MyTimer()
    min_vertex = None
    min_f = 100
    if test_Seg == 5:
        t1 = 10
        t2 = 11
        t3 = 12
        t4 = 13
        t5 = 14
        dt = 2
        for x1 in range(t1-dt,t1+dt+1):
            for x2 in range(t2-dt,t2+dt+1):
                for x3 in range(t3-dt,t3+dt+1):
                    for x4 in range(t3-dt,t3+dt+1):
                        for x5 in range(t3-dt,t3+dt+1):
                            #vertex = {'I': [x1, x2, x3, 't']}
                            vertex = {'I_x1': [x1, 105, 't'], 'I_x2': [x2, 106], 'I_x3': [x3, 107], 'I_x4': [x4, 108], 'I_x5': [x5, 109]}
                            f = test.solve(vertex)
                            if f < min_f:
                                min_f = f
                                min_vertex = copy.deepcopy(vertex)
                                print(min_f, min_vertex)
    elif test_Seg == 1:
        t1 = 15
        t2 = 105
        dt = 3
        for x1 in range(t1-dt,t1+dt+1):
            for x2 in range(t2-dt,t2+dt+1):
                vertex = {'I': [x1, x2,]}
                f = test.solve(vertex)
                if f < min_f:
                    min_f = f
                    min_vertex = copy.deepcopy(vertex)
                    print(min_f, min_vertex)

    print("min result:",min_f, min_vertex)
    t = timer1.getTime(kind='real')
    exit()
    

    #f = test.solve(vertex)
    #f = test.solve_by_gdxIn(vertex)

    exit()
    
    t = timer1.getTime(kind='real')
    min_FId = 1000
    min_vertex = None

    for t in range(1,241,20):
        vertex = {'I_x12':[t],'I_x3':[],'I_x45':[]}
        f = test.solve(vertex)
        #f, db = test.solve_outDB(vertex)
        if f < min_FId:
            min_FId = f
            min_vertex = copy.deepcopy(vertex)
    
    print(min_FId,min_vertex)
    t = timer1.getTime(kind='real')
    dbName = test.solve_outDB(min_vertex)
    #test.plot_Db(gdxfile=dbName)
    
    if '-s' in sys.argv :
        saveVertexInfo = {}
        saveVertexInfo['Nz'] = test_Nz
        saveVertexInfo['StartDirect'] = StartDirect
        saveVertexInfo['vertex'] = vertex
        saveFileName = os.getcwd()+'\\TEST-Nz'+str(test_Nz)+'-'+str(len(Uncertain_parameter_name))+'seg-result'+str(datetime.now().strftime('%Y-%m-%d'))+'.pkl'
        with open(saveFileName, 'wb') as f:
            pickle.dump(saveVertexInfo, f)
            print(f"save file {saveFileName} finished.")
    
    exit()