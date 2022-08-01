from gams import *
import sys
import os
import pickle
import numpy as np
import math
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

class CHS():
    def __init__(self,dt=0.001,dx=0.1):
        self.total_sim_time = 350
        self.x_grid_num = int(1/dx)+1

        self.ws = GamsWorkspace(working_directory = os.getcwd())
        self.model_name = "CHS"
        self.is_show_solve_result = False

        self.dt = dt
        self.dx = dx
        self.BC_L = 1.0 #B.C. x = x1
        self.BC_R = 1.0 #B.C. x = x20
        self.IC_X = np.zeros((self.x_grid_num)) #I.C. t = t1
        self.temp_lb = 0.4
        self.drawing = 1
        self.plot_t = 1

        # upper-point
        self.U_point_x = 0.5
        self.U_point_y = 0
        self.U_point_z = 0.7
        # upper-left bound co-eff
        self.vec_T = [0,1,0] #T-line
        self.vec_XL = [1,0,0.4475] #X-line
        self.vec_XR = [1,0,-0.4475]  #X-line
        self.UL_vec = np.cross(self.vec_T, self.vec_XL)
        self.UR_vec = np.cross(self.vec_T, self.vec_XR)

        normal_point = np.array([self.U_point_x,self.U_point_y,self.U_point_z])
        self.UL_d = -normal_point.dot(self.UL_vec)
        self.UR_d = -normal_point.dot(self.UR_vec)

        # upper-bound
        self.upper_bound = None
        self.lower_bound = None

        self.Nz = 1
        self.Seg = 1
        self.uncertainParaList = ["Tout"]
        self.startDirect = {"Tout": 1}
        self.GAMSModelInstance = None
        self.cp = None
        self.job = None
        self.cp_file = None

        self.queue_len = 10
        self.min_FId = 1000

    def set_Nz(self,nz):
        self.Nz = nz

    def set_Seg(self,seg):
        self.Seg = seg

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
            file_name += '_nz'+str(self.Nz)+'.gms'
        else:
            file_name = file
        print("GAMSfile:", file_name)
        if os.path.isfile('./'+file_name):
            self.jobGAMSfile = file_name
        else:
            raise FileNotFoundError("GAMSModelInstance file :"+file_name+"not found")

    def get_activated_point(self,rec):
        return int(rec[0][1:])-1, int(rec[1][1:])-1

    def check_bound(self,data):
        x = np.linspace(0,1,self.x_grid_num)
        t = np.linspace(0,self.total_sim_time,self.total_sim_time)
        
        total_num = 0
        over_num = 0
        under_num = 0

        print(data.shape)

        for i in range(0,len(t)):
            for j in range(1,20,1) :
                total_num += 1
                if data[i][j] >= self.upper_bound[i][j]:
                    #print(f"x:{x[j]}, y:{t[i]}, z:{data[i][j]}")
                    over_num += 1
                elif data[i][j] <= self.lower_bound[i][j] :
                    #print(f"x:{x[j]}, y:{t[i]}, z:{data[i][j]}")
                    under_num += 1
                else :
                    pass

        print(f"total point: {total_num}, over point: {over_num}, under point: {under_num}")
        print("over rate:",(over_num/total_num))
        print("under rate:",(under_num/total_num))
        if over_num+under_num > 0:
            return False
        else:
            return True

    def set_IC_X(self,magnitude,theta,omega,ctype='const',showfig=False):
        if ctype == 'const' :
            ary = np.ones(self.x_grid_num)
            ary = ary*magnitude
        elif ctype == 'eq' :
            ary = np.linspace(0,1,num=self.x_grid_num,endpoint=True)
            ary = -ary**2+ary
            ary += self.temp_lb
        elif ctype == 'linear' :
            ary = np.linspace(0.5,magnitude,21)
        elif ctype == 'phasor':
            ary = np.linspace(0,omega*math.pi,num=self.x_grid_num,endpoint=True)
            ic = magnitude*np.cos(ary+theta*math.pi/2)
            ic = ic+self.temp_lb
            ary = ic
        elif ctype == 'rev-phasor':
            ary = np.linspace(0,omega*math.pi,num=self.x_grid_num,endpoint=True)
            ic = magnitude*np.cos(ary+theta*math.pi/2)
            ic = ic+self.temp_lb
            ary = ic
        elif ctype == "FId":
            ary = np.linspace(0,omega*math.pi,num=self.x_grid_num,endpoint=True)
            ary = ary+theta*math.pi/2
        
        self.IC_X = ary
        x = np.linspace(0,1,num=self.x_grid_num)
        fig, ax = plt.subplots(1,1,sharex=True,sharey=True)
        ax.plot(x,self.IC_X)
        ax.set_title('initial condition')
        ax.set_xlabel('x grid point') # X axis data label
        ax.set_ylabel('init temp data') # Y axis data label
        fig.tight_layout()
        if showfig :
            plt.show()
        plt.close('all')

    def set_upper_bound(self):
        x_grid = np.linspace(0,1,self.x_grid_num)
        t_grid = np.linspace(0,self.total_sim_time,self.total_sim_time)

        mid_point = self.x_grid_num//2
        x_left_grid = x_grid[0:mid_point+1]
        X1 , Y1 = np.meshgrid(x_left_grid,t_grid)
        Z_L = (-self.UL_vec[0]*X1-self.UL_vec[1]*Y1-self.UL_d)*1.0/self.UL_vec[2]
        x_right_grid = x_grid[mid_point:]
        X2 , Y2 = np.meshgrid(x_right_grid,t_grid)
        Z_R = (-self.UR_vec[0]*X2-self.UR_vec[1]*Y2-self.UR_d)*1.0/self.UR_vec[2]
        
        self.upper_bound =np.concatenate((Z_L,Z_R[:,1:]),axis=1)
        self.upper_bound += 0.08
    
    def set_lower_bound(self):
        x_grid = np.linspace(0,1,self.x_grid_num)
        t_grid = np.linspace(0,self.total_sim_time,self.total_sim_time)

        mid_point = self.x_grid_num//2
        x_left_grid = x_grid[0:mid_point+1]
        X1 , Y1 = np.meshgrid(x_left_grid,t_grid)
        Z_L = (-self.UL_vec[0]*X1-self.UL_vec[1]*Y1-self.UL_d)*1.0/self.UL_vec[2]
        x_right_grid = x_grid[mid_point:]
        X2 , Y2 = np.meshgrid(x_right_grid,t_grid)
        Z_R = (-self.UR_vec[0]*X2-self.UR_vec[1]*Y2-self.UR_d)*1.0/self.UR_vec[2]
        
        self.lower_bound =np.concatenate((Z_L,Z_R[:,1:]),axis=1)
        self.lower_bound -= 0.15

    def showInitSetting(self):
        x = np.linspace(0,1,num=self.x_grid_num)

        fig, ax = plt.subplots(1,1,sharex=True,sharey=True)
        ax.plot(x,self.IC_X)
        ax.scatter(x,self.IC_X,c='red')
        ax.set_xlabel(r'$X$')
        ax.set_xlim(0.0,1.0)
        ax.set_ylabel(r'$T$')
        ax.set_ylim(0.3,0.8)
        if self.upper_bound is not None:
            ax.plot(x, self.upper_bound[0][:],color='k',linestyle='--',linewidth=2)
        if self.lower_bound is not None:
            ax.plot(x, self.lower_bound[0][:],color='k',linestyle='--',linewidth=2)
        fig.tight_layout()
        plt.show()

        print("init u",self.IC_X)
        print("lower bound of u",self.lower_bound[0][:])
        print("upper bound of u",self.upper_bound[0][:])

        plt.close('all')

    def set_GAMSModelInstance(self):
        self.cp = self.ws.add_checkpoint()
        self.job = self.ws.add_job_from_file(self.jobGAMSfile)
        self.job.run(checkpoint=self.cp)
        mi = self.cp.add_modelinstance()

        opt = self.ws.add_options()
        opt.all_model_types = "conopt"
        opt.solprint = 1
        print("fdopt:\t\t",opt.fdopt)
        print("fddelta:\t",opt.fddelta)

        modifiers = []

        u = mi.sync_db.add_variable("u", 2 ,VarType.Free)
        up_limit = mi.sync_db.add_parameter("up_limit", 2, "upper limit plane constraint")
        modifiers.append(GamsModifier(u, UpdateAction.Upper, up_limit))
        
        lo_limit = mi.sync_db.add_parameter("lo_limit", 2, "lower limit plane constraint")
        modifiers.append(GamsModifier(u, UpdateAction.Lower, lo_limit))
        
        for p in self.uncertainParaList:
            vertex = mi.sync_db.add_parameter(p+'_vertex', 1, p+"_d direction vertex")
            modifiers.append(GamsModifier(vertex))

        mi.instantiate(self.model_name+" using nlp maximizing obj",modifiers,opt)
        for i in range(1,self.total_sim_time+1):
            for j in range(2,self.x_grid_num):
                up_limit.add_record(('x'+str(j),'t'+str(i))).value = self.upper_bound[i-1][j-1]
                lo_limit.add_record(('x'+str(j),'t'+str(i))).value = self.lower_bound[i-1][j-1]

        self.GAMSModelInstance = mi

    def update_vertex(self,p_name,shift_list):
        parameter = self.GAMSModelInstance.sync_db.get_parameter(p_name+'_vertex')
        parameter.clear()

        v_len_V = len(shift_list)
        d_cur_V = self.startDirect[p_name]
        shift_t = 0
        for rec in self.job.out_db.get_parameter(p_name+'_vertex'):
            key = "".join([chat for chat in rec.keys[0] if chat.isdigit()])
            if shift_t < v_len_V:
                if shift_list[shift_t] == 't':
                    shift_t = 10000
                elif int(key) >=shift_list[shift_t] :
                    d_cur_V = -d_cur_V
                    shift_t += 1
            parameter.add_record('t'+key).value = d_cur_V

    def solve(self,vertex):
        mi = self.GAMSModelInstance

        for p in self.uncertainParaList:
            self.update_vertex(p,vertex[p])
        
        mi.solve()
        i = 0
        while mi.model_status != 2 or mi.solver_status != 1:
            print("resolve time: ",i)
            mi.solve()
            i += 1
        if self.is_show_solve_result:
            print("Scenario vertex=" + str(vertex) + ":")
            print("  Modelstatus: " + str(mi.model_status))
            print("  Solvestatus: " + str(mi.solver_status))
            print("  Obj: " + str(mi.sync_db.get_variable("obj")[()].level))
        result = mi.sync_db.get_variable("obj")[()].level
        return result

    def solve_revised(self,ver):
        FId = self.solve(ver)
        re_FId = self.solve(ver)
        while re_FId != FId :
            FId = re_FId
            re_FId = self.solve(ver)
        FId_r = self.solve_result2(ver)

        return FId_r

    def solve_outDB(self,vertex,savefig=False):
        mi = self.GAMSModelInstance

        for p in self.uncertainParaList:
            self.update_vertex(p,vertex[p])
        
        mi.solve()
        print("Scenario vertex=" + str(vertex) + ":")
        print("  Modelstatus: " + str(mi.model_status))
        print("  Solvestatus: " + str(mi.solver_status))
        print("  Obj: " + str(mi.sync_db.get_variable("obj")[()].level))
        FId = mi.sync_db.get_variable("obj")[()].level
        
        if savefig:
            dirPath = os.getcwd()+'\\Nz'+str(self.Nz)+'-result'
            if os.path.isdir(dirPath) == False:
                os.mkdir(dirPath)
            export_gdx_name = 'Nz'+str(self.Nz)+'-'+str(len(self.uncertainParaList))+'Seg-result.gdx'
            gdxPath = os.path.join(dirPath, export_gdx_name)
            mi.sync_db.export(gdxPath)
            print("save gdx:",gdxPath)
            
        t_len = 350
        x_grid = np.linspace(0,1,self.x_grid_num)
        t_grid = np.linspace(0,self.total_sim_time,self.total_sim_time)/1000
        X, Y = np.meshgrid(x_grid,t_grid)

        data = np.empty((t_len,self.x_grid_num))
        i, j = 0, 0
        err_limit = 1e-7
        upper_activated = []
        lower_activated = []
        for rec in mi.sync_db["u"] :
            data[i][j] = rec.level
            i += 1
            if i >= data.shape[0] :
                i = 0
                j += 1
            if abs(rec.level-rec.lower) <= err_limit:
                lower_activated.append((rec.key(1),rec.key(0)))
            elif abs(rec.level-rec.upper) <= err_limit:
                upper_activated.append((rec.key(1),rec.key(0)))

        left_bc = []
        right_bc = []
        for rec in mi.sync_db["bc_L"]:
            left_bc.append(rec.level)
        for rec in mi.sync_db["bc_R"]:
            right_bc.append(rec.level)
        T_out = np.empty((t_len,self.x_grid_num))
        i, j = 0, 0
        for rec in mi.sync_db["Tout"]:
            T_out[i][j] = rec.level
            i += 1
            if i >= T_out.shape[0] :
                i = 0
                j += 1
        if len(self.uncertainParaList) == 1:
            for j in range(1,self.x_grid_num):
                T_out[:,j] = T_out[:,0]
        
        fig_CV, [ax_CV_bcL, ax_CV_bcR] = plt.subplots(nrows=1, ncols=2,figsize=(8, 3))
        ax_CV_List = [ax_CV_bcL,ax_CV_bcR]
        bc_ub = 0.6
        bc_lb = 0.3

        ax_CV_bcR.plot(t_grid,right_bc,color='steelblue',linestyle='-',label=r'bc$_{R}$')
        ax_CV_bcR.hlines(bc_ub,xmin=0,xmax=t_len,color='black',linestyle='--')
        ax_CV_bcR.hlines(bc_lb,xmin=0,xmax=t_len,color='black',linestyle='--')
        ax_CV_bcR.set_ylim(0.25, 0.65)
        ax_CV_bcR.set_xlim(0, 0.35)
        ax_CV_bcR.yaxis.set_ticks(np.linspace(0.25,0.65,num=5,endpoint=True))
        ax_CV_bcR.ticklabel_format(style='sci', axis='x', scilimits=(-3,-3), useMathText=True)
        ax_CV_bcR.set_xlabel(r'$\tau$')
        ax_CV_bcR.set_ylabel("bc$_{R}$")

        ax_CV_bcL.plot(t_grid,left_bc,color='steelblue',linestyle='-',label=r'bc$_{L}$')
        ax_CV_bcL.hlines(bc_ub,xmin=0,xmax=t_len,color='black',linestyle='--')
        ax_CV_bcL.hlines(bc_lb,xmin=0,xmax=t_len,color='black',linestyle='--')
        ax_CV_bcL.set_ylim(0.25, 0.65)
        ax_CV_bcL.set_xlim(0, 0.35)
        ax_CV_bcL.yaxis.set_ticks(np.linspace(0.25,0.65,num=5,endpoint=True))
        ax_CV_bcL.ticklabel_format(style='sci', axis='x', scilimits=(-3,-3), useMathText=True)
        ax_CV_bcL.set_xlabel(r'$\tau$')
        ax_CV_bcL.set_ylabel("bc$_{L}$")

        for ax_cv in ax_CV_List:
            ax_cv.grid(color='y', linestyle='--', linewidth=1, alpha=0.3)
            ax_cv.legend()

        plt.tight_layout()
        figName ='Nz'+str(self.Nz)+'-'+str(self.Seg)+'Seg-CV.png'
        figPath = os.path.join(dirPath, figName)
        fig_CV.savefig(figPath,dpi=200)
        print("save fig:",figName)

        fig_T_out = plt.figure()
        ax_T_out = fig_T_out.add_subplot()
        axtout = ax_T_out.pcolor(X,Y,T_out,cmap='seismic',shading='auto')
        #surface = ax_T_out.plot_surface(X, Y, T_out, color='blue', rstride=1, cstride=1, cmap='coolwarm')
        cbar_Tout= fig_T_out.colorbar(axtout,shrink=1.0,aspect=20)
        cbar_Tout.ax.get_yaxis().labelpad = 15
        cbar_Tout.ax.set_ylabel("uncertain $T_{out}$", fontsize=12, rotation=270)
        #ax_T_out.set_title("uncertain $T_{out}$", fontsize=12)
        ax_T_out.set_xlabel("$X$",fontsize=12)
        ax_T_out.set_xlim(0, 1)
        ax_T_out.ticklabel_format(style='sci', axis='y', scilimits=(-3,-3), useMathText=True)
        ax_T_out.set_ylabel(r'$\tau$',fontsize=12)
        figName ='Nz'+str(self.Nz)+'-'+str(self.Seg)+'Seg-UP.png'
        figPath = os.path.join(dirPath, figName)
        fig_T_out.savefig(figPath,dpi=200)
        print("save fig:",figName)

        fig_T = plt.figure()
        ax = fig_T.add_subplot(projection='3d')

        surface = ax.plot_surface(X, Y, data, color='blue', rstride=1, cstride=1, cmap='coolwarm')
        cbar_T = fig_T.colorbar(surface,shrink=1.0,aspect=20)
        cbar_T.ax.get_yaxis().labelpad = 15
        cbar_T.ax.set_ylabel("slab temp $T$", fontsize=12, rotation=270)
        #ax.set_title("slab temperature distribute", fontsize=16)
        ax.set_xlabel(r'$X$',fontsize=12)
        ax.set_ylabel(r'$\tau$',fontsize=12)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(-3,-3), useMathText=True)
        ax.set_zlabel(r'$T$',fontsize=12)

        upper_activated_num = len(upper_activated)
        lower_activated_num = len(lower_activated)
        if self.upper_bound is not None:
            ax.plot_surface(X, Y, self.upper_bound, color='black', rstride=1, cstride=1, alpha=0.3)
            for act in upper_activated:
                i, j = self.get_activated_point(act)
                ax.scatter(x_grid[j],t_grid[i],data[i][j]+0.001,color="r",marker='2',s=500)
                print(f"upper_activated point x:{j+1} {x_grid[j]},t:{i+1} {t_grid[i]} value{data[i][j]}")
            
        if self.lower_bound is not None:
            ax.plot_surface(X, Y, self.lower_bound, color='black', rstride=1, cstride=1, alpha=0.3)
            for act in lower_activated:
                i, j = self.get_activated_point(act)
                ax.scatter(x_grid[j],t_grid[i],data[i][j]+0.001,color="r",marker='1',s=500)
                print(f"lower_activated point x:{j+1} {x_grid[j]},t:{i+1} {t_grid[i]} value{data[i][j]}")
        print("upper activate:",upper_activated_num)
        print("lower activate:",lower_activated_num)

        #set data lim region
        ax.set_xlim3d(0,1)
        ax.set_ylim3d(0,0.35)
        ax.set_zlim3d(0.2,0.8)
        #set view point
        ax.azim = -74
        ax.dist = 10
        ax.elev = 27

        plt.tight_layout()
        figName ='Nz'+str(self.Nz)+'-'+str(self.Seg)+'Seg-SV.png'
        figPath = os.path.join(dirPath, figName)
        fig_T.savefig(figPath,dpi=200)
        print("save fig:",figName)

        plt.show()
        plt.close('all')

if __name__ == '__main__' :
    test_Nz = 1
    dx = 0.1
    GAMSfile = None
    loadInfo = None
    StartDirect = None
    Uncertain_parameter_name = None
    vertex = None

    if '-i' in sys.argv :
        loadFileIndex = sys.argv.index('-i')+1
        try:
            loadFileName = sys.argv[loadFileIndex]
            with open(loadFileName, 'rb') as f:
                loadInfo= pickle.load(f)
        except:
            raise IndexError(f"loadFileIndex sys.argv:[{loadFileIndex}] not exist")
        test_Nz = loadInfo['Nz']
        Uncertain_parameter_name = []
        for var_name in loadInfo['StartDirect']:
            Uncertain_parameter_name.append(var_name)
        StartDirect = loadInfo['StartDirect']
        vertex = loadInfo['vertex']
        print(f"load Info from  {loadFileName} suc")
        print('Nz',test_Nz)
        print('Uncertain_parameter_name',Uncertain_parameter_name)
        print('StartDirect',StartDirect)
        print('Vertex',vertex)

        var_num = len(Uncertain_parameter_name)
        GAMSfile = "CHS_modelinstance_Tout_nz"+str(test_Nz)+"_"+str(var_num)+"seg.gms"
    else:
        test_Nz = 15
        seg = 3

        if seg == 1:
            Uncertain_parameter_name = ['Tout']
        elif seg == 3:
            Uncertain_parameter_name = ['Tout_x234','Tout_x567','Tout_x8910']
        elif seg == 5:
            Uncertain_parameter_name = ['Tout_x2','Tout_x34','Tout_x567','Tout_x89','Tout_x10']
        StartDirect = {}
        for up in Uncertain_parameter_name:
            StartDirect[up] = -1
        print(StartDirect)
        print(Uncertain_parameter_name)
        GAMSfile = "CHS_modelinstance_Tout_nz"+str(test_Nz)+"_"+str(seg)+"seg.gms"

    test = CHS()
    test.set_Nz(test_Nz)
    test.set_Seg(seg)
    test.set_uncertain_para(Uncertain_parameter_name)
    test.setStartDirect(StartDirect)
    test.set_GAMSfile(GAMSfile)
    test.set_IC_X(0.3,3,1,ctype='eq',showfig=False)
    test.set_upper_bound()
    test.set_lower_bound()
    test.set_GAMSModelInstance()
    test.is_show_solve_result = True

    #vertex = {'Tout': ['t']}
    vertex = {'Tout_x234': [342], 'Tout_x567': [], 'Tout_x8910': [342]}
    
    f = test.solve(vertex)
    #test.solve_outDB(vertex,savefig=True)

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