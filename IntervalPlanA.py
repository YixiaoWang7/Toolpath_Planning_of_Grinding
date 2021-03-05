import numpy as np
import copy
import elkai
from functools import partial
import random
from GeometricLocalPath import path_arc

def tsp_path(initial_point, path_set, path_l_set, path_sf_set, clock_set):
    n_v = len(path_l_set)
    M = np.zeros(shape = (n_v+1, n_v+1))

    for i in range(0,n_v+1):
        for j in range(0,n_v+1):
            if j != i:
                if i == 0:
                    M[i, j] = ((initial_point[0]-path_sf_set[j-1][0][0])**2+(initial_point[1]-path_sf_set[j-1][0][1])**2)**0.5
                else:
                    if j == 0:
                        M[i, j] = ((initial_point[0]-path_sf_set[i-1][1][0])**2+(initial_point[1]-path_sf_set[i-1][1][1])**2)**0.5
                    else:
                        M[i, j] = ((path_sf_set[i-1][1][0]-path_sf_set[j-1][0][0])**2+(path_sf_set[i-1][1][1]-path_sf_set[j-1][0][1])**2)**0.5

    
    if len(path_set)>2:
        Mint = np.around(M*100).astype(int)
        path_order = elkai.solve_int_matrix(Mint) 
    else:
        path_order = np.arange(0,len(path_set))
    
    TL_tsp = path_l_set[path_order[0]]
    for i in range(len(path_order)-1):
        TL_tsp += M[path_order[i],path_order[i+1]]

    T_path_tsp = []
    T_clock_tsp = []
    if path_order[0] == 0:
        T_path_tsp = [[np.linalg.norm(np.asarray([initial_point[0],initial_point[1]]))], [np.arctan2(initial_point[1],initial_point[0])]]
        T_clock_tsp.append(0)
        for i in range(len(path_order)-1):
            T_path_tsp[0].extend(path_set[path_order[i+1]-1][0])
            T_path_tsp[1].extend(path_set[path_order[i+1]-1][1])
            T_clock_tsp.extend(clock_set[path_order[i+1]-1])
    else:
        print('path_order[0] is not zero.')
        return (TL_tsp, path_order, T_path_tsp)
    return (TL_tsp, path_order, T_path_tsp, T_clock_tsp)


class GA_path():
    def __init__(self, fun, initial_num, variable_size, crossover_rate = 0.8, mutation_rate = 0.01):
        self.fun = fun
        self.initial_num = initial_num
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.variable_size = variable_size
        self.population = []
        self.value_population = []
        p_tem = np.arange(variable_size[1])
        for i in range(self.initial_num):
            tp = []
            for j in range(variable_size[0]):
                tp.append(random.choice(p_tem))
            self.population.append(tp)
            self.value_population.append(self.fun(tp))
        
        
    def crossover(self):
        index = np.arange(len(self.population))
        random.shuffle(index)
        for i in range(int(len(self.population)/2)):
            if random.random() < self.crossover_rate:
                pos = random.randint(0,self.variable_size[0]-1)
                m = index[2*i]
                n = index[2*i+1]
                self.population.append(self.population[m][:pos] + self.population[n][pos:])
                self.population.append(self.population[n][:pos] + self.population[m][pos:])
                self.value_population.append(self.fun(self.population[-2]))
                self.value_population.append(self.fun(self.population[-1]))
    
    def mutation(self):
        for i in range(len(self.population)):
            sp = copy.deepcopy(self.population[i])
            flag = 0
            for j in range(self.variable_size[0]):
                if random.random() < self.mutation_rate:
                    sp[j] = random.choice(np.arange(self.variable_size[1]))
                    flag = 1
            if flag == 1:
                self.population.append(sp)
                self.value_population.append(self.fun(sp))
        
    def get_population(self):
        return self.population
    
    def select(self):
        value_array = np.asarray(self.value_population)
        value_array = np.exp(value_array - min(value_array))
        p = value_array/np.sum(value_array)
        try:
            index = np.random.choice(len(self.population), self.initial_num, replace=False, p=value_array/np.sum(value_array))
        except:
            print(p,value_array)
        max_value = np.max(self.value_population)
        max_pop = copy.deepcopy(self.population[self.value_population.index(max_value)])
        self.population = [self.population[i] for i in index]
        self.value_population = [self.value_population[i] for i in index]
        if np.max(self.value_population) < max_value:
            self.population.append(max_pop)
            self.value_population.append(max_value)
            min_index = self.value_population.index(np.min(self.value_population))
            self.population.pop(min_index)
            self.value_population.pop(min_index)
            
    def iteration(self, loop, end_loop = 20):
        max_value_his = []
        max_value = np.max(self.value_population)
        max_value_his.append(max_value)
        count = 0
        for _ in range(loop):
            self.crossover()
            self.mutation()
            self.select()
            if max_value == np.max(self.value_population):
                count += 1
            else:
                count = 0
            max_value = np.max(self.value_population)
            max_value_his.append(max_value)
            if count == end_loop:
                break
        max_pop = self.population[self.value_population.index(max_value)]
        return max_value, max_value_his, max_pop
    
    
def ga_tsp_path(a, initial_point, line_set):
    path_set = []
    path_l_set = []
    path_sf_set = []
    clock_set = []
    for i in range(len(line_set)):
        [new_path, l, sf, clock ] = path_arc(line_set[i], flag = a[i]+1)
        path_set.append(new_path)
        path_l_set.append(l)
        path_sf_set.append(sf)
        clock_set.append(clock)

    (TL_tsp, _, _,_) = tsp_path(initial_point, path_set, path_l_set, path_sf_set, clock_set)
    return -TL_tsp/1000