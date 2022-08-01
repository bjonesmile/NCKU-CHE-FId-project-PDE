import sys
import os
import psutil
import math
import copy
import random
import numpy as np
import itertools
from datetime import datetime

import PTC as problem
import mytimer as timer

proble_Nz = 6
segNum = 5
Nz = proble_Nz
nz = 1

isPDE = True
MIN_POINT_DISTANCE = 5
MIN_RANGE_TO_NEXT_ACT = None
SIM_TIME = 240
GAMMA = 1.0
try_limit = 10
FId_Precision = None

StartDirect = {}
if(segNum == 1):
    Uncertain_parameter_name = ['I']
if(segNum == 3):
    Uncertain_parameter_name = ['I_x12','I_x3','I_x45']
if(segNum == 5):
    Uncertain_parameter_name = ['I_x1','I_x2','I_x3','I_x4','I_x5']
nzLimitList = {}
for up in Uncertain_parameter_name:
    #StartDirect[up] = random.choice([+1,-1])
    StartDirect[up] = +1
    #nzLimitList[up] = Nz*nz+1
    nzLimitList[up] = 4

GAMSfile = "PTC_modelinstance_I_nz"+str(proble_Nz)+"_"+str(segNum)+"seg.gms"
gdxInGAMSFILE = "PTC_modelinstance_I_nz"+str(proble_Nz)+"_"+str(segNum)+"seg_gdxin.gms"

solver = problem.PTC(dt=30,dx=100,nfe=5)
solver.set_Nz(proble_Nz)
solver.set_Seg(segNum)
solver.set_uncertain_para(Uncertain_parameter_name)
solver.setStartDirect(StartDirect)
solver.set_GAMSfile(GAMSfile)
solver.set_upper_bound("T_o")
solver.set_lower_bound("T_o")
solver.set_gdxInGAMSFILE(gdxInGAMSFILE)
solver.set_GAMSModelInstance()
solver.is_show_solve_result = True
print(solver.uncertainParaList)

#x10: ['Tout_x2', 'Tout_x3', 'Tout_x4', 'Tout_x5', 'Tout_x6', 'Tout_x7', 'Tout_x8', 'Tout_x9', 'Tout_x10']
class State(object):
    def __init__(self,ucps = ['I_x1','I_x2','I_x3','I_x4','I_x5']):
        self.current_value = 0.0
        # For the first root node, the index is 0 and the game should start from 1
        self.current_round_index = 0
        self.cumulative_choices = []
        for ucp in ucps:
            self.cumulative_choices.append({'var':ucp,'choices':[]})
        self.available_switch_point = None
        self.available_switch_point_num = None

    def get_current_value(self):
        return self.current_value

    def set_current_value(self, value):
        self.current_value = value

    def get_current_round_index(self):
        return self.current_round_index

    def set_current_round_index(self, turn):
        self.current_round_index = turn

    def get_cumulative_choices(self):
        return self.cumulative_choices

    def set_cumulative_choices(self, choices):
        self.cumulative_choices = choices

    def is_terminal(self):
        # check whether current state match termination conditions
        is_finish = False

        for choices in self.cumulative_choices :
            if 't' in choices['choices'] :
                return True

        switch_times = 0
        nz_list = {}
        for item in self.cumulative_choices :
            nz_list[item['var']] = len(item['choices'])
            switch_times += len(item['choices'])
        if isPDE :
            num_limit = 0
            for key in nz_list.keys():
                if nz_list[key] == nzLimitList[key]:
                    num_limit += 1
                elif nz_list[key] > nzLimitList[key]:
                    raise ValueError(f"uncertain parameter {key} swich times: {nz_list[key]} over limit Nz {nzLimitList[key]}")
            if num_limit == len(nzLimitList):
                is_finish = True
        else:
            if sum(nz_list.values()) >= Nz :
                is_finish = True
            
        return is_finish

    def compute_reward(self):
        input_vertex = {}
        for choices in self.cumulative_choices :
            if choices['var'] in solver.uncertainParaList :
                input_vertex[choices['var']] = choices['choices']
            
        FId = solver.solve(input_vertex)
        if FId_Precision != None:
            return -round(math.pow(FId,2),FId_Precision)
        else:
            return -math.pow(FId,2)

    def set_available_switch_point(self):
        is_finish = self.is_terminal()
        max_st_point = -1
        for item in self.cumulative_choices:
            st_point = None

            cur_var = item['var']
            cur_choices = item['choices']
            if cur_var not in Uncertain_parameter_name:
                raise ValueError("uncertain parameter ",cur_var," not in uncertain list.")

            if 't' in cur_choices:
                is_finish = True
                break
            if len(cur_choices) == 0:
                st_point = MIN_POINT_DISTANCE
            else:
                st_point = cur_choices[-1]
            
            if st_point > max_st_point:
                max_st_point = st_point

        av_point = 1
        if MIN_RANGE_TO_NEXT_ACT != None:
            av_point = max_st_point+MIN_RANGE_TO_NEXT_ACT
        else:
            if max_st_point != 1:
                av_point = max_st_point+MIN_POINT_DISTANCE

        available_switch_point = []
        if is_finish == False:
            for up in Uncertain_parameter_name:
                up_item = next((item for item in self.cumulative_choices if item['var'] == up), None)
                if isPDE and len(up_item['choices']) >= nzLimitList[up]:
                    continue
                cur_available_switch_point = list(np.arange(av_point,SIM_TIME+1,MIN_POINT_DISTANCE,dtype=np.int32))
                # add max_st_point action to the uncertain parameter which hasnt concerned
                if len(up_item['choices']) == 0 :
                    cur_available_switch_point.insert(0,max_st_point)
                    if MIN_POINT_DISTANCE != 1 and max_st_point <= MIN_POINT_DISTANCE:
                        cur_available_switch_point.insert(0,1)
                else:
                    if up_item['choices'][-1] != max_st_point:
                        cur_available_switch_point.insert(0,max_st_point)
                available_switch_point.extend([{'var':up,'point':t} for t in cur_available_switch_point])

        for up in Uncertain_parameter_name:
            up_item = next((item for item in self.cumulative_choices if item['var'] == up), None)
            if isPDE and len(up_item['choices']) >= nzLimitList[up]:
                continue
            available_switch_point.append({'var':up,'point':'t'})

        self.available_switch_point = available_switch_point
        self.available_switch_point_num = len(available_switch_point)
        #print(self.available_switch_point)
        #print("set available_switch_point finished")

    def get_available_switch_point(self):
        if self.available_switch_point_num == None :
            self.set_available_switch_point()
        return self.available_switch_point_num

    def get_next_state_with_random_choice(self):
        var_list = solver.uncertainParaList
        max_t_list = []
        for item in self.cumulative_choices:
            if 't' in item['choices']:
                max_t_list.clear()
                max_t_list.append(SIM_TIME+1)
            else:
                max_t_list.append(max(item['choices'],default=-1))
        av_t = max(max_t_list)
        
        if av_t <= 0:
            av_t = 1
        elif av_t > SIM_TIME:
            raise ValueError(f"get random next state but current max time av_t {av_t} is OVER total SimTime {SIM_TIME}")

        act_list = []
        for item in self.cumulative_choices:
            if isPDE and len(item['choices']) == nzLimitList[item['var']]:
                continue
            cur_act_list = list(np.arange(av_t+MIN_POINT_DISTANCE,SIM_TIME+1,MIN_POINT_DISTANCE,dtype=np.int32))
            if len(item['choices']) == 0 :
                cur_act_list.insert(0,av_t)
            elif item['choices'][-1] != av_t:
                cur_act_list.insert(0,av_t)
            act_list.extend([{'var':item['var'],'point':t} for t in cur_act_list])
            act_list.append({'var':item['var'],'point':'t'})

        random_choice = copy.deepcopy(random.choice(act_list))
        del act_list

        next_state = State(solver.uncertainParaList)
        next_state.set_current_value(random_choice)
        next_state.set_current_round_index(self.current_round_index + 1)
        var = random_choice['var']
        current_cumulative_choices = copy.deepcopy(self.cumulative_choices)
        spec_item = next((item for item in current_cumulative_choices if item['var'] == var), None)
        spec_item['choices'] = spec_item['choices'] + [random_choice['point']]
        next_state.set_cumulative_choices(current_cumulative_choices)

        #print("random next state: ",next_state)
        return next_state

    def get_next_state_with_spec_choice(self,choice):
        #print("spec choice", choice)
        next_state = State(solver.uncertainParaList)
        next_state.set_current_value(choice)
        next_state.set_current_round_index(self.current_round_index + 1)
        var = choice['var']
        current_cumulative_choices = copy.deepcopy(self.cumulative_choices)
        spec_item = next((item for item in current_cumulative_choices if item['var'] == var), None)
        spec_item['choices'] = spec_item['choices'] + [choice['point']]
        next_state.set_cumulative_choices(current_cumulative_choices)

        #print("spec next state: ",next_state)
        return next_state

    def __repr__(self):
        return "State: {}, value: {}, round: {}, choices: {}".format(
            hash(self), self.current_value, self.current_round_index,
            self.cumulative_choices)

    def showState(self):
        string_choices = ""
        for item in self.cumulative_choices:
            if len(item['choices']) != 0:
                string_choices += str(item) + ", "
        print("State info: value: {}, round: {}, choices: {}".format(
            self.current_value, self.current_round_index, string_choices))

class Node(object):
    """
    Node of MCTS's tree structure, incluinf of parent and children, and function used to calculate quality(UCB) and current state of simulation.
    """

    def __init__(self):
        self.parent = None
        self.children = []

        self.visit_times = 0
        self.quality_value = 0.0
        
        self.best_reward = -10
        self.best_vertex = None
        self.best_action_round = 1

        self.state = None

    def set_state(self, state):
        self.state = state

    def get_state(self):
        return self.state

    def get_parent(self):
        return self.parent

    def set_parent(self, parent):
        self.parent = parent

    def get_children(self):
        return self.children

    def get_visit_times(self):
        return self.visit_times

    def set_visit_times(self, times):
        self.visit_times = times

    def visit_times_add_one(self):
        self.visit_times += 1

    def get_quality_value(self):
        return self.quality_value

    def set_quality_value(self, value):
        self.quality_value = value

    def quality_value_add_n(self, n):
        self.quality_value += n

    def is_all_expand(self):
        return len(self.children) == self.state.get_available_switch_point()

    def add_child(self, sub_node):
        sub_node.set_parent(self)
        self.children.append(sub_node)

    def __repr__(self):
        if self.visit_times == 0:
            ratio_value = 0
        else:
            ratio_value = self.quality_value/self.visit_times
        return "Node: {}, Q/N: {}/{}, ratio: {},state: {}".format(
            hash(self), self.quality_value, self.visit_times, ratio_value, self.state)


def tree_policy(node):
    # Check if the current node is the leaf node
    while node.get_state().is_terminal() == False:
        if node.is_all_expand():
            node = best_child(node, True)
        else:
            # initial consider all 1 level action and expand sub node
            sub_node = init_expand(node,len(node.children))
            return sub_node

        # Return the leaf node
        return node


def default_policy(node):
    """
    simulation step: using unform random method to expand node until the state terminal.
    """

    # Get the state of the game
    current_state = node.get_state()
    # First consider directly terminal state
    if node.get_visit_times() == 0:
        # give a random uncertain parameter and action 't' for first visit
        random_var = random.choice(solver.uncertainParaList)
        current_state.get_next_state_with_spec_choice({'var':random_var,'point':'t'})
    else:    
        # Run until the simulation over
        while current_state.is_terminal() == False:
            # Pick one random action to proceed and get next state
            current_state = current_state.get_next_state_with_random_choice()
            #print("default_policy current state:",current_state)

    final_state_reward = current_state.compute_reward()
    vertex = copy.deepcopy(current_state.get_cumulative_choices())

    if round(final_state_reward,6) >= node.best_reward :
        node.best_reward = final_state_reward
        node.best_vertex = vertex
        node.best_action_round = current_state.get_current_round_index()
    return final_state_reward, vertex

def improve_policy(node,t):
    """
    simulation step: using unform random method to expand node until the state terminal.
    but need to preform better reward than current best reward.
    """
    # if this node already approach to best state
    if t >= try_limit:
        return node.best_reward

    # Get the state of the game
    current_state = node.get_state()
    # First consider directly terminal state
    if node.get_visit_times() == 0:
        random_var = random.choice(solver.uncertainParaList)
        current_state.get_next_state_with_spec_choice({'var':random_var,'point':'t'})
    else:    
        # Run until the game over
        while current_state.is_terminal() == False:
            # Pick one random action to play and get next state
            current_state = current_state.get_next_state_with_random_choice()
            #print("improve_policy current state:",current_state)

    final_state_reward = current_state.compute_reward()

    if round(final_state_reward,6) >= node.best_reward :
        node.best_reward = final_state_reward
        node.best_action_round = current_state.get_current_round_index()
    else :
        return improve_policy(node,t+1)
    return final_state_reward

def init_expand(node,action_index):
    """
    first search all possible action from root node and expand it 
    """
    possible_action = node.get_state().available_switch_point[action_index]
    new_state = node.get_state().get_next_state_with_spec_choice(possible_action)

    sub_node = Node()
    sub_node.set_state(new_state)
    node.add_child(sub_node)

    return sub_node

def expand(node):
    """
    input a node and expand from the node using uniform possibility random choose action add sub-node.
    """

    tried_sub_node_states = [
        sub_node.get_state().get_cumulative_choices() for sub_node in node.get_children()
    ]

    new_state = node.get_state().get_next_state_with_random_choice()

    # Check until get the new state which has the different action from others
    while new_state.get_cumulative_choices() in tried_sub_node_states:
        new_state = node.get_state().get_next_state_with_random_choice()

    sub_node = Node()
    sub_node.set_state(new_state)
    node.add_child(sub_node)

    return sub_node


def best_child(node, is_exploration):
    """
    use UCB alorithm, judge exploration and exploitation and select the subnode with highest value, when predicting just consider Q/N。
    """

    # TODO: Use the min float value
    best_score = -sys.maxsize
    best_sub_node = None

    # Travel all sub nodes to find the best one
    for sub_node in node.get_children():

        # Ignore exploration for inference
        if is_exploration:
            C = 1 / math.sqrt(2.0)
        else:
            C = 0.0

        # UCB = quality / times + C * sqrt(2 * ln(total_times) / times)
        left = sub_node.get_quality_value() / sub_node.get_visit_times()
        right = 2.0 * math.log(node.get_visit_times()) / sub_node.get_visit_times()
        score = left + C * math.sqrt(right)

        if score > best_score:
            best_sub_node = sub_node
            best_score = score
    
    return best_sub_node

def best_child_optimal(node,gamma):
    # TODO: Use the min float value
    best_score = -sys.maxsize
    best_avg = -sys.maxsize
    best_sub_node = None

    # Travel all sub nodes to find the best one
    for sub_node in node.get_children():

        left = (1-gamma)*sub_node.get_quality_value() / sub_node.get_visit_times()
        right = gamma*(sub_node.best_reward)
        
        score = left + right

        #print(left,right)
        #print(sub_node.get_state().get_cumulative_choices(),score)

        if score > best_score:
            best_sub_node = sub_node
            best_score = score
        if left > best_avg:
            best_avg = left
            print("avg best Node",sub_node.get_state().get_cumulative_choices())
    
    #print("choose best child with score:",best_score)
    return best_sub_node

def backup(node, reward):
    """
    MCTS Backpropagation step: input the node need to expend，feedback to expend node and upsteam path nodes and rew new data.
    """

    # Update util the root node
    while node != None:
        # Update the visit times
        node.visit_times_add_one()

        # Update the quality value
        node.quality_value_add_n(reward)

        # Change the node to the parent node
        node = node.parent


def monte_carlo_tree_search(node,computation_budget = 1000):
    """
    MCTS contain four steps: Selection, Expansion, Simulation, Backpropagation。
    first step and scond step: use tree policy find sub-node worth to expolre.
    third step: use default policy to choose random path until terminated and get reward.
    fourth step: Backpropagation, renew reward to all node the random path passed.
    for predicting, just according to Q value to choose exploitation node.
    """
    min_FId = -2
    min_vertex = None

    print("computation budget:",computation_budget)
    # Run as much as possible under the computation budget
    for i in range(computation_budget):
        
        # 1. Find the best node to expand
        expand_node = tree_policy(node)
        #print("\r selected sub node:",expand_node)
        # 2. Random run to add node and get reward
        reward, vertex = default_policy(expand_node)
        # reward = improve_policy(expand_node,0)
        if reward >= min_FId :
            min_vertex = copy.deepcopy(vertex)
            min_FId = reward
        
        if i%10 == 0:
            dispState = []
            for item in min_vertex:
                if len(item['choices']) != 0:
                    dispState.append(item)
            print("\rcompute times: {}/{}\tpercentage: {}\nnow min FId: {}, min vertex: {}"\
            .format(i,computation_budget, str(round(i/computation_budget*100,2))+'%', round(math.sqrt(-min_FId),5), dispState), end="")
            sys.stdout.flush()
        
        # 3. Update all passing nodes with reward
        backup(expand_node, reward)
        
    # N. Get the best next node
    best_next_node = best_child_optimal(node, GAMMA)

    return best_next_node

def show_current_sub_node_state(node,period):
    tried_sub_node_states = [
        sub_node for sub_node in node.get_children()
    ]

    with open("search_subnode_history_round"+str(node.get_state().get_current_round_index()+1)+".txt",'a') as wfile :
        print(f"period: {period}",file=wfile)
        for sub_node in tried_sub_node_states :
            print(sub_node,file= wfile)
        print('\n',file=wfile)

def check_computation_budget(current_node):
    cur_vertex_len = 0
    for uncertain in current_node.get_state().get_cumulative_choices():
        cur_vertex_len += len(uncertain['choices'])
    if isPDE :
        return cur_vertex_len == sum(nzLimitList.values())-1
    else:
        return cur_vertex_len == Nz-1
         
def main():
    const_t = int(10)
    # Create the initialized state and initialized node
    init_state = State(solver.uncertainParaList)
    #init_state.set_cumulative_choices([{'var': 'I_x12', 'choices': [10]}, {'var': 'I_x3', 'choices': []}, {'var': 'I_x45', 'choices': []}])
    #init_state.set_cumulative_choices([{'var': 'I_x1', 'choices': [10]}, {'var': 'I_x2', 'choices': [15]}, {'var': 'I_x3', 'choices': [15]}, {'var': 'I_x4', 'choices': [20]}, {'var': 'I_x5', 'choices': [20]},])
    """init_state.set_cumulative_choices([{'var': 'Tout_x2', 'choices': [150]}, 
                                       {'var': 'Tout_x34', 'choices': [150]}, 
                                       {'var': 'Tout_x567', 'choices': []},
                                       {'var': 'Tout_x89', 'choices': [150]},
                                       {'var': 'Tout_x10', 'choices': [150]}])"""
    init_node = Node()
    init_node.set_state(init_state)
    current_node = init_node

    # Set the rounds to simulate
    i = 1
    while current_node.get_state().is_terminal() == False :
        current_node.get_state().set_available_switch_point()
        possible_action_num = current_node.get_state().get_available_switch_point()
        print("possible actions to next state",possible_action_num)
        print("simulate Nz level: {}".format(i))
        
        computation_budget = int(possible_action_num*const_t)
        if check_computation_budget(current_node):
            computation_budget = possible_action_num+1
        print("this level computation_budget:",computation_budget)
        current_node = monte_carlo_tree_search(current_node,computation_budget=computation_budget)
        print("Choose node: {}".format(current_node))
        i += 1
        #input()

    return current_node.get_state().get_cumulative_choices()

def mk_vertex_generator(Nz):
    ary = np.arange(1,SIM_TIME,step=MIN_POINT_DISTANCE)
    print(len(ary))
    a_list = ary.tolist()
    a_list.append('t')
    vertex_list = itertools.combinations(a_list,Nz)

    return vertex_list

def ex_search(ex_nz, limit=None):
    min_FId = 1000
    min_vertex = None
    limit_nz = ex_nz
    if limit is not None :
        limit_nz = limit

    for nz in range(1,ex_nz+1):
        if nz == limit_nz :
            print("reach limit:",limit_nz)
            break
        cur_generator = mk_vertex_generator(nz)
        for item in cur_generator :
            FId = solver.solve(item)
            if FId <= min_FId:
                re_FId = solver.solve(item)
                while re_FId != FId :
                    FId = re_FId
                    re_FId = solver.solve(item)
                min_FId = FId
                min_vertex = list(item).copy()
            print("current min vertex",min_vertex,min_FId)

    print(min_vertex,min_FId)
    return min_vertex

def multi_process_ex_search(ex_nz, limit=None, nr_workers=4):
    limit_nz = ex_nz
    if limit is not None :
        limit_nz = limit

    for nz in range(1,ex_nz+1):
        if nz == limit_nz :
            print("reach limit:",limit_nz)
            break
        cur_generator = mk_vertex_generator(nz)
        solver.scen_solve_main(ver_gen=cur_generator,uncertain_para=[Uncertain_parameter_name])

        print(solver.min_FId)
    return solver.min_FId


if __name__ == "__main__":
    op_cmd = "MCTS"
    print("mode:",op_cmd)
    
    if op_cmd == "MCTS" :
        timer1 = timer.MyTimer()
        MCTS_result = main()
        print("Simulation: Nz",proble_Nz,"Seg",segNum,"suc")
        print(StartDirect)
        print(MCTS_result)
        t = timer1.getTime(kind='real')

        result_vertex = {}
        for item in MCTS_result:
            result_vertex[item['var']] = item['choices']
        fid = solver.solve(result_vertex)
        print(fid, result_vertex)

        import pickle
        saveFileName = os.getcwd()+'\\MCTS-Nz'+str(proble_Nz)+'-'+str(len(Uncertain_parameter_name))+'seg-result'+str(datetime.now().strftime('%Y-%m-%d'))+'.pkl'
        with open(saveFileName, 'wb') as f:
            saveVertexInfo = {}
            saveVertexInfo['Nz'] = proble_Nz
            saveVertexInfo['StartDirect'] = StartDirect
            saveVertexInfo['vertex'] = result_vertex
            saveVertexInfo['FId'] = fid
            saveVertexInfo['using-time'] = t
            pickle.dump(saveVertexInfo, f)
            print(f"save file {saveFileName} finished.")

    elif op_cmd == "test":
        test_ver = {'I_x12': [65, 150, 280], 'I_x3': [70, 145, 285], 'I_x45': [65, 150, 280]} 
        test = solver.solve(test_ver)
        #solver.is_show_solve_result = True
        timer1 = timer.MyTimer()
        tau_k = 10
        x = 320
        min_f = 100
        min_vertex = None
        for x1 in range(x-tau_k, x+tau_k+1):
            for x2 in range(x-tau_k, x+tau_k+1):
                for x3 in range(x-tau_k, x+tau_k+1):
                    test_ver = {'Tout_x234': [x1, x2, 280], 'Tout_x567': [70, x3, 285], 'Tout_x8910': [x1, x2, 280]} 
                    f = solver.solve(test_ver)
                    if f < min_f:
                        min_f = f
                        min_vertex = copy.deepcopy(test_ver)
                        print(min_f, min_vertex)
        t = timer1.getTime(kind='real')
        print(min_f, min_vertex)
        print("suc")

    exit()