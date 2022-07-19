import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
import subprocess
import os

# ========================================== #
# =============== CELL CLASS =============== #
# ========================================== #
class Cell:
    def __init__(self, log_length0, gfp0, lambda0, q0, time0=0., cell_id = 0, parent_id=-1):
        self.parent_id = parent_id
        self.cell_id = cell_id
        self.length = [np.exp(log_length0)]  # s(t)
        self.log_length = [log_length0]      # x(t) = x0 + int lambda dt
        self.gfp = [gfp0]
        self.lt = [lambda0]
        self.qt = [q0]
        self.time = [time0]
        self.segment = []

    def to_df(self, n=1, start=0):
        df = pd.DataFrame({   "cell_id": ([self.cell_id]*len(self.time))[start::n],
                                "time_min": self.time[start::n],
                                "parent_id": ([self.parent_id]*len(self.time))[start::n],
                                "log_length": self.log_length[start::n], 
                                "gfp": self.gfp[start::n],
                                "lt": self.lt[start::n],
                                "qt": self.qt[start::n]})
        if len(self.segment)>0:
            df['segment']=self.segment[start::n]
        return df

# =============== pd dataframe to cells =============== #

def df2cells(dataset, time="time_min", 
            log_length="log_length", gfp="gfp", 
            cell_id="cell_id", parent_id="parent_id", 
            lt=None, qt=None):
    """ 
    dataset (pandas data frame as read from csv file) to list of Cell instances 
    either using columns with measurment noise or without
    """
    cell_list = []
    last_cell = ""
    for _, row in dataset.iterrows(): 
        if row[cell_id] != last_cell:
            if parent_id == None:   p = -1
            else:                   p = row[parent_id]

            if lt == None:          lambda0 = None
            else:                   lambda0 = row[lt]

            if qt == None:          q0 = None
            else:                   q0 = row[qt]

            new_cell = Cell(row[log_length], row[gfp], 
                        lambda0, q0, 
                        time0=row[time],
                        cell_id=row[cell_id], 
                        parent_id=p)
            cell_list.append(new_cell)
        else:
            cell_list[-1].log_length.append(row[log_length])
            cell_list[-1].gfp.append(row[gfp])
            cell_list[-1].time.append(row[time])

            if lt != None: 
                cell_list[-1].lt.append(row[lt])

            if qt != None:
                cell_list[-1].qt.append(row[qt])

        last_cell = row[cell_id]
    return cell_list



def ggp_df2cells(dataset, time="time", 
            log_length="mean_x", gfp="mean_g", 
            lt="mean_l", qt="mean_q",
            cov_xx="cov_xx",
            cov_gg="cov_gg",
            cov_ll="cov_ll",
            cov_qq="cov_qq",
            cell_id="cell_id", 
            parent_id="parent_id", 
            lane=None):
    """ 
    dataset (pandas data frame as read from csv file) to list of Cell instances, m
    written for ggp output
    """
    cell_list = []
    last_cell = ""
    for _, row in dataset.iterrows(): 
        if row[cell_id] != last_cell:
            if lane !=None:
                c = str(row[lane])+'.'+str(row[cell_id])
                p = str(row[lane])+'.'+str(row[parent_id])
            else:
                c = str(row[cell_id])
                p = str(row[parent_id])

            lambda0 = row[lt]
            q0 = row[qt]

            new_cell = Cell(row[log_length], row[gfp], 
                        lambda0, q0, 
                        time0=row[time],
                        cell_id=c, 
                        parent_id=p)
            cell_list.append(new_cell)
            cell_list[-1].cov_xx = []
            cell_list[-1].cov_gg = []
            cell_list[-1].cov_ll = []
            cell_list[-1].cov_qq = []

        else:
            cell_list[-1].log_length.append(row[log_length])
            cell_list[-1].gfp.append(row[gfp])
            cell_list[-1].time.append(row[time])

            cell_list[-1].lt.append(row[lt])
            cell_list[-1].qt.append(row[qt])

        cell_list[-1].cov_xx.append(row[cov_xx])
        cell_list[-1].cov_gg.append(row[cov_gg])
        cell_list[-1].cov_ll.append(row[cov_ll])
        cell_list[-1].cov_qq.append(row[cov_qq])

        last_cell = row[cell_id]
    return cell_list


# ============ Cell tree ============ #

def get_cell_by_cell_id(cell_list, cell_id):
    for cell in cell_list:
        if cell_id == cell.cell_id:
            return cell 
    return None


def get_cell_by_parent_id(cell_list, parent_id, ignore=None):
    dcells = []
    for cell in cell_list:
        if parent_id == cell.parent_id and cell.parent_id != ignore:
            dcells.append(cell) 
    return dcells + [None, None]


def cell_paths(curr_cell, cell_list):
    def cell_paths_recr(curr_cell, cell_list):
        if curr_cell == None:
            return
        current_path.append(curr_cell)
        dcells = get_cell_by_parent_id(cell_list, curr_cell.cell_id)

        if dcells[0] == None:
            path.append(copy.deepcopy(current_path))
        else:
            cell_paths_recr(dcells[0], cell_list)
            cell_paths_recr(dcells[1], cell_list)
        current_path.pop()

    path = []
    current_path = []
    cell_paths_recr(curr_cell, cell_list)
    return path


# =============== create super cells =============== #

def cocanate_cells(cell_list):
    new = copy.deepcopy(cell_list[0])
    new.cell_id = str(new.cell_id )
    for cell in cell_list[1:]:
        new.cell_id += "_"+ str(cell.cell_id)
        new.time = np.append(new.time, cell.time)

        new.log_length = np.append(new.log_length, cell.log_length)
        new.gfp = np.append(new.gfp, cell.gfp)
        new.lt = np.append(new.lt, cell.lt)
        new.qt = np.append(new.qt, cell.qt)
    return new


def cocanate_ggp_cells(cell_list):
    new = copy.deepcopy(cell_list[0])
    new.cell_id = str(new.cell_id )
    for cell in cell_list[1:]:
        new.cell_id += "_"+ str(cell.cell_id)
        new.time = np.append(new.time, cell.time)

        new.log_length = np.append(new.log_length, cell.log_length)
        new.gfp = np.append(new.gfp, cell.gfp)
        new.lt = np.append(new.lt, cell.lt)
        new.qt = np.append(new.qt, cell.qt)

        new.cov_xx = np.append(new.cov_xx, cell.cov_xx)
        new.cov_gg = np.append(new.cov_gg, cell.cov_gg)
        new.cov_ll = np.append(new.cov_ll, cell.cov_ll)
        new.cov_qq = np.append(new.cov_qq, cell.cov_qq)
    return new



# ========================================== #
# =============== SIMULATION =============== #
# ========================================== #

def single_ou_step(dt, mean, gamma, var, x):
    noise = np.sqrt(var) * np.random.normal(loc=0, scale=1) * np.sqrt(dt)
    return x - gamma*(x-mean)*dt + noise

def growth(cell, dt, lt):  
    # calculate next step
    next_step = cell.log_length[-1] + lt*dt

    # save everything
    cell.lt.append(lt)
    cell.log_length.append(next_step)
    cell.length.append(np.exp(next_step))
    return cell


def gfp_production(cell, dt, qt, beta):
    # calculate next step
    next_step = cell.gfp[-1] + cell.length[-1]*qt*dt - cell.gfp[-1]*beta*dt   

    cell.qt.append(qt)
    cell.gfp.append(next_step)
    return cell 

def cell_divsion(cell, var_dx, var_dg, no_cells, cell_division_model, dt):
    if cell_division_model=="binomial":
        x = cell.log_length[-1]
        g = cell.gfp[-1]
        x_d1 = np.random.normal(loc=x-np.log(2), scale=np.sqrt(var_dx)) 
        if var_dg==0:
            g_d1 = g * np.exp(x_d1 - x)
        else:
            if g <= 0:
                g_d1 = 0
            else:   
                g_d1 = np.random.binomial(g/var_dg, np.exp(x_d1 - x))*var_dg

        x_d2 = np.log(np.exp(x)-np.exp(x_d1))
        
        if g <= 0:
            g_d2 = 0
        else:
            g_d2 = g - g_d1

        cell1 = Cell(x_d1,g_d1, cell.lt[-1], cell.qt[-1],
                        time0 = cell.time[-1]+dt,
                        cell_id = no_cells + 1, parent_id=cell.cell_id)

        cell2 = Cell(x_d2,g_d2, cell.lt[-1], cell.qt[-1],
                        time0 = cell.time[-1]+dt,
                        cell_id = no_cells + 2, parent_id=cell.cell_id)

    elif cell_division_model=="gauss":
        cell1 = Cell(np.random.normal(loc=cell.log_length[-1] - np.log(2), scale=np.sqrt(var_dx)),
                        np.random.normal(loc=cell.gfp[-1]/2, scale=np.sqrt(var_dg)), cell.lt[-1], cell.qt[-1],
                        time0 = cell.time[-1]+dt,
                        cell_id = no_cells + 1, parent_id=cell.cell_id)

        cell2 = Cell(np.random.normal(loc=cell.log_length[-1] - np.log(2), scale=np.sqrt(var_dx)),
                        np.random.normal(loc=cell.gfp[-1]/2, scale=np.sqrt(var_dg)), cell.lt[-1], cell.qt[-1], 
                        time0 = cell.time[-1]+dt,
                        cell_id = no_cells + 2, parent_id=cell.cell_id)
    # print(cell1.cell_id, cell.time[-1], cell1.log_length[0], cell.log_length[-1])
    # print(cell2.cell_id, cell.time[-1], cell2.log_length[0], cell.log_length[-1])

    # _cell_div.append( [cell1.log_length[0] - cell.log_length[-1] + np.log(2) ,  cell1.gfp[0] - cell.gfp[-1]/2 ])
    # _cell_div.append( [cell2.log_length[0] - cell.log_length[-1] + np.log(2) ,  cell2.gfp[0] - cell.gfp[-1]/2 ])

    no_cells += 2 
    return cell1, cell2, no_cells

def is_cell_division(cell, mode, division_log_length, division_time, division_addition):
    if mode == "timer":
        if cell.time[-1] - cell.time[0] > division_time:
            return True
        else: 
            return False
    if mode == "sizer":
        if cell.log_length[-1] > division_log_length:
            return True
        else: 
            return False
    if mode == "adder":
        if np.exp(cell.log_length[-1]) - np.exp(cell.log_length[0]) > np.exp(division_addition):
            return True

# ============= FULL SIMULATION ============= # 

def simulate_cells(dt, n_cells, parameter_set, div_mode, 
                    division_log_length=None, 
                    division_time=None, 
                    division_addition=None, 
                    log_length0=None,
                    gfp0=None, 
                    tree=True,
                    tmax=np.inf, 
                    verbose=False, 
                    discard_cells=0, 
                    cell_division_model="binomial"):
    if gfp0 == None:
        gfp0 = 3*parameter_set['mean_q']/parameter_set['mean_lambda']
    if log_length0 == None:
        log_length0 = division_log_length-np.log(2)

    cell_queue = [Cell(log_length0, gfp0, parameter_set['mean_lambda'], parameter_set['mean_q'])]

    cells_simulated = []
    no_cells = 0   # total number of cells (in queue and calculated)
    t_total = 0    # the largest t in the current set of cells

    while len(cells_simulated) < n_cells and t_total < tmax:
        if len(cells_simulated)>0:
            idx = get_current_leafs_idx(cell_queue, cells_simulated)
            cell_index = np.random.choice(idx)
        else:
            cell_index = 0

        # print_simulated(cells_simulated)
        # print_queue(cell_queue,cell_index)

        # simulate a random cell in queue
        cell = copy.deepcopy(cell_queue[cell_index])

        # --------------------------------------------------------------- #
        # Simulation on single cell level
        # --------------------------------------------------------------- #
        while True:
            cell.time.append(cell.time[-1]+dt)

            t_total = np.max([t_total, cell.time[-1]])
            if t_total > tmax:
                break

            q_ou = single_ou_step(dt,   parameter_set['mean_q'], 
                                        parameter_set['gamma_q'], 
                                        parameter_set['var_q'], 
                                        cell.qt[-1]) 

            cell = gfp_production(cell, dt, q_ou, parameter_set['beta'])

            lambda_ou = single_ou_step(dt,  parameter_set['mean_lambda'], 
                                            parameter_set['gamma_lambda'], 
                                            parameter_set['var_lambda'], 
                                            cell.lt[-1]) 
            cell = growth(cell, dt, lambda_ou)

            if is_cell_division(cell, div_mode, division_log_length, division_time, division_addition):
                # save the simulated cell
                if discard_cells==0:
                    cells_simulated.append(cell)
                else:
                    discard_cells -= 1
                
                # calc. new init conditions for 2 daugter cells
                cell1, cell2, no_cells = cell_divsion(cell, parameter_set['var_dx'], 
                                                            parameter_set['var_dg'], 
                                                            no_cells, cell_division_model, dt)

                # remove the simulated cell from queue and add the new ones 
                cell_queue.pop(cell_index)
                cell_queue.append(cell1)
                if tree and discard_cells==0:
                    cell_queue.append(cell2)
                break
            else:
                pass

        if verbose:
            progress_bar_n = np.around(len(cells_simulated)/n_cells*20).astype(int)
            progress_bar = '='*progress_bar_n + ' '*(20-progress_bar_n)
            print("\r|", progress_bar,  "| Progress {:3.0f}%".format(len(cells_simulated)/n_cells*100), " No of cells: ", len(cells_simulated), end='')  
            print('')

    return cells_simulated

def simulate_cells_segments(dt, n_cells, parameter_sets, div_mode, t_segment,
                    division_log_length=None, 
                    division_time=None, 
                    division_addition=None, 
                    log_length0=None,
                    gfp0=None, 
                    tree=True,
                    tmax=np.inf, 
                    verbose=True, 
                    cell_division_model="binomial"):
    parameter_set = parameter_sets[0]

    if gfp0 == None:
        gfp0 = 3*parameter_set['mean_q']/parameter_set['mean_lambda']
    if log_length0 == None:
        log_length0 = division_log_length-np.log(2)

    cell_queue = [Cell(log_length0, gfp0, parameter_set['mean_lambda'], parameter_set['mean_q'])]

    cells_simulated = []
    no_cells = 0   # total number of cells (in queue and calculated)
    t_total = 0    # the largest t in the current set of cells

    segment=0

    while len(cells_simulated) < n_cells and t_total < tmax:

        if len(cells_simulated)>0:
            idx = get_current_leafs_idx(cell_queue, cells_simulated)
            cell_index = np.random.choice(idx)
        else:
            cell_index = 0

        # print_simulated(cells_simulated)
        # print_queue(cell_queue,cell_index)

        # simulate a random cell in queue
        cell = copy.deepcopy(cell_queue[cell_index])

        # --------------------------------------------------------------- #
        # Simulation on single cell level
        # --------------------------------------------------------------- #
        while True:
            cell.time.append(cell.time[-1]+dt)

            if cell.time[-1]>t_segment:
                segment=1
            else:
                segment=0

            parameter_set = parameter_sets[segment]

            cell.segment.append(segment)

            t_total = np.max([t_total, cell.time[-1]])
            if t_total > tmax:
                break

            q_ou = single_ou_step(dt,   parameter_set['mean_q'], 
                                        parameter_set['gamma_q'], 
                                        parameter_set['var_q'], 
                                        cell.qt[-1]) 

            cell = gfp_production(cell, dt, q_ou, parameter_set['beta'])

            lambda_ou = single_ou_step(dt,  parameter_set['mean_lambda'], 
                                            parameter_set['gamma_lambda'], 
                                            parameter_set['var_lambda'], 
                                            cell.lt[-1]) 
            cell = growth(cell, dt, lambda_ou)

            if is_cell_division(cell, div_mode, division_log_length, division_time, division_addition):
                # save the simulated cell
                cells_simulated.append(cell)
                
                # calc. new init conditions for 2 daugter cells
                cell1, cell2, no_cells = cell_divsion(cell, parameter_set['var_dx'], 
                                                            parameter_set['var_dg'], 
                                                            no_cells, cell_division_model, dt)
                cell1.segment.append(segment)
                cell2.segment.append(segment)

                # remove the simulated cell from queue and add the new ones 
                cell_queue.pop(cell_index)
                cell_queue.append(cell1)
                if tree:
                    cell_queue.append(cell2)
                break
            else:
                pass
        progress_bar_n = np.around(len(cells_simulated)/n_cells*20).astype(int)
        progress_bar = '='*progress_bar_n + ' '*(20-progress_bar_n)
        if verbose:
            print("\r|", progress_bar,  "| Progress {:3.0f}%".format(len(cells_simulated)/n_cells*100), " No of cells: ", len(cells_simulated), end='')  
    print('')
    return cells_simulated

# =============== Helper functions ===============#

def print_queue(cell_queue, cell_index):
    print('In queue (total=', len(cell_queue), ')', end=' ', sep='')
    for i in range(len(cell_queue)):
        if i == cell_index:
            print('|', cell_queue[i].cell_id, '<-',  cell_queue[i].parent_id ,'|', end=' ', sep='')
        else:
            print(cell_queue[i].cell_id, end='  ', sep='')
    print('')

def print_simulated(cells_simulated):
    print('Simulated:', end=' ')
    for i in range(len(cells_simulated)):
        print(cells_simulated[i].cell_id, end=' ', sep='')
    print('')

def get_ids(cell):
    id_list = []
    for i in range(len(cell)):
        id_list.append(cell[i].cell_id)
    return id_list

def get_current_leafs_idx(cell_queue, cells_simulated):
    sim_ids = get_ids(cells_simulated)
    idx = []
    for i in range(len(cell_queue)):
        if cell_queue[i].parent_id in sim_ids:
            idx.append(i)
    return idx
        


# ====================== BUILD DATA ====================== #
def mk_mising_dir(path_name):
    if not os.path.exists(path_name):
        os.mkdir(path_name)
    return path_name

def get_next_file_name(out_dir, force_index=None):
    sample = out_dir.split('/')[-1]
    print(force_index)

    if force_index !=None:
        new_dir = os.path.join(out_dir, sample+"_{:d}".format(force_index),'')
        new_file = os.path.join(new_dir, sample+"_{:d}".format(force_index)+".csv")
        print("force i=",force_index)
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)
        return new_dir, new_file

    for i in range(1000):
        new_dir = os.path.join(out_dir, sample+"_{:d}".format(i),'')
        new_file = os.path.join(new_dir, sample+"_{:d}".format(i)+".csv")
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)
            return new_dir, new_file


def write_param_file(filename, parameters, 
                    init_scale=0.5, 
                    bound_scales=[0.2, 5], 
                    fixed_p=[], bound_p=[]):
    with open(filename, "w") as fin:
        fin.write("# Generated config file for simulated data\n")
        for k, v in parameters.items():
            if k in fixed_p:
                fin.write("{:s} = {:.2E}\n".format(k, v))
            elif k in bound_p:
                fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format(k, v,
                                                                            v*init_scale, 
                                                                            v*bound_scales[0], 
                                                                            v*bound_scales[1]))
            else:            
                fin.write("{:s} = {:.2E}, {:.2E}\n".format(k, v, np.max([v*init_scale, 1e-11])))

def write_param_file_noise(filename, parameters, log_noise=0, rel_step=None, fixed=[]):
    with open(filename, "w") as fin:
        fin.write("# Generated config file for simulated data with noise: " + str(log_noise) + "\n")
        for param_name in parameters.keys():
            if param_name not in fixed:
                if log_noise>0:
                    init = np.exp(np.log(parameters[param_name]) + np.random.normal(0, scale=log_noise))
                else:
                    init = parameters[param_name]
                fin.write("{:s} = {:.2E}, {:.2E}\n".format(param_name, init, init*rel_step))
            else:
                init = parameters[param_name]
                fin.write("{:s} = {:.2E}\n".format(param_name, init))

    return filename

def write_csv_config(filename, segment=None, lane=None):
    with open(filename, "w") as fin:
        fin.write("# Generated config file for simulated data\n")
        fin.write("time_col = time_min  \n")
        fin.write("length_col = log_length_noise \n")
        fin.write("length_islog = true \n")
        fin.write("fp_col = gfp_noise \n")
        if lane!=None:
            fin.write("cell_tags = " + lane + ", cell_id \n")
            fin.write("parent_tags = " + lane + ", parent_id \n")
        else:
            fin.write("cell_tags = cell_id \n")
            fin.write("parent_tags = parent_id \n")

        fin.write("rescale_time = 1 \n")
        if segment!=None:
            fin.write("segment_col = " + segment + " \n")
    return filename


def build_data_set_fixed_dt(cells_simulated, var_x, var_g, dt, atol=1e-4):
    dataset = pd.DataFrame()
    for i in range(len(cells_simulated)):
        next_celldf = cells_simulated[i].to_df(1)
        next_celldf['log_length_noise'] = next_celldf['log_length'] + np.random.normal(loc=np.zeros_like( next_celldf['log_length']), scale=np.sqrt(var_x))
        next_celldf['gfp_noise'] = next_celldf['gfp'] + np.random.normal(loc=np.zeros_like( next_celldf['gfp']), 
                                                                        scale=np.sqrt(var_g))
        dataset = dataset.append(next_celldf)
    
    return dataset[ np.abs( (dataset['time_min'] + atol/2) % dt) < atol ] 


def build_data_set_scale_gfp_noise(cells_simulated, var_x, var_g, dt, atol=1e-4):

    dataset = pd.DataFrame()
    for i in range(len(cells_simulated)):
        next_celldf = cells_simulated[i].to_df(1)
        next_celldf['log_length_noise'] = next_celldf['log_length'] + np.random.normal(loc=np.zeros_like( next_celldf['log_length']), scale=np.sqrt(var_x))
        next_celldf['gfp_noise'] = next_celldf['gfp'] + np.random.normal(loc=np.zeros_like( next_celldf['gfp']), 
                                                                        scale=np.sqrt(var_g * np.abs(next_celldf['gfp'])))
        dataset = dataset.append(next_celldf)
    
    return dataset[ np.abs( (dataset['time_min'] + atol/2) % dt) < atol ] 


# =============== RUN COMMAND =============== #
def suggest_run_command(directory, filename, 
                        modes ="-s -m -p", 
                        out_dir=None, 
                        t=1e-1, 
                        param_file = "parameters.txt", 
                        add_flag='', 
                        binary = "../bin/gfp_gaussian" ):
    cmd = binary + " -c " + os.path.join(directory, "csv_config.txt") + \
        " -b " + os.path.join(directory, param_file) + \
        " -t " + str(t) +  " -i " + filename + " -l 0 " + ' ' + add_flag + ' '
    if out_dir==None:
        return cmd + modes
    else:
        return cmd + modes + " -o " + out_dir


# =============== PLOT =============== #
def plot_cells(cells, n_steps=1):
    _, axes = plt.subplots(2, figsize=(10,7))
    ax = axes.ravel()

    for j in range(len(cells[::n_steps])):
        cell = copy.deepcopy(cells[j])
        cell.time = np.array(cell.time)

        if n_steps>1:
            ax[0].set_title("log length (showing every {:d}th cell)".format(n_steps))
            ax[1].set_title("gfp (showing every {:d}th cell)".format(n_steps))

        else:
            ax[0].set_title("log length")
            ax[1].set_title("gfp")

        # ax[0].set_ylim([1.2, 2.2])
        
        if len(cells[::n_steps]) <20:
            ax[0].axvline(cell.time[-1], ls='--', color='tab:blue')
            ax[1].axvline(cell.time[-1], ls='--', color='tab:orange')

        if j ==0:
            ax[0].plot(cell.time, np.array(cell.log_length), label='log length', color='tab:blue')
            ax[1].plot(cell.time, np.array(cell.gfp), color='tab:orange', label='gfp')

        else:
            ax[0].plot(cell.time, np.array(cell.log_length), color='tab:blue')
            ax[1].plot(cell.time, np.array(cell.gfp), color='tab:orange')

    for j in range(2):
        ax[j].legend()
    plt.show()




# =================== Stand alone simulation =================== #

def get_final_files(directory):
    entries = os.listdir(directory)
    final_files = []
    for e in entries:
        if e.endswith("_final.csv"):
            final_files.append(os.path.join(directory,e))
    return final_files


def get_final_parameter_set(filename, return_type=str):

    parameter_names =  ["mean_lambda", 
                        "gamma_lambda",
                        "var_lambda",
                        "mean_q",
                        "gamma_q",
                        "var_q",
                        "beta",
                        "var_x",
                        "var_g",
                        "var_dx",
                        "var_dg"]

    final_set = []
    with open(filename) as fin:
        next(fin)
        for k in parameter_names:
            line = fin.readline()
            final_set.append(float(line.split(',')[-1].strip()))
    return dict(zip(parameter_names, final_set))

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Simulate cells with OU processes')

    parser.add_argument('--i',
                        dest='infile',
                        help='input file, final file from gfp_gaussian run or directory with final file(s) from gfp_gaussian run',
                        required=False,
                        default=None)
    parser.add_argument('--o',
                        dest="outdir",
                        help="output dir",
                        required=True)
    parser.add_argument('--l',
                        dest="lanes",
                        help="simulate lineages instead of a random tree of cells",
                        required=False,
                        default=1,
                        type=int)
    parser.add_argument('--n',
                        dest="n_cells",
                        help="number of cells per lane",
                        required=False,
                        default=250,
                        type=int)
    parser.add_argument('--m',
                        dest="div_mode",
                        help="division mode",
                        required=False,
                        default="sizer")
    parser.add_argument('--d',
                        dest="discard_cells",
                        help="discard first n cells",
                        required=False,
                        default=0, 
                        type=int)
    parser.add_argument('--g',
                        dest="genealogy",
                        help="genealogy: 't' = tree (default) or 's' = single lineage",
                        required=False,
                        default="t")
    parser.add_argument('--t',
                        dest="dt",
                        help="dt between measurment points",
                        required=False,
                        default=1e-3,
                        type=float)

    args = parser.parse_args()


    # ========== Simulation parameters ========== #

    dt = 1e-2
    dt_measument = 3 # in minutes

    overwrite_params = { 
                    "var_x": 0.0,
                    "var_g": 0.0,
                    "var_dx": 0.0,
                    "var_dg": 0.0}

    n_cells = args.n_cells # number of cells that will be simulated

    div_mode = args.div_mode

    division_log_length = 1+np.log(2)   # for sizer: division, when log_length hits division_log_length
    division_time = 60 - 1e-10          # for timer: division, when cell cycle time hits division_time
    division_addition = np.log(2)       # for adder: divsion, when division_addition in log_length was added in cell cycle


    if os.path.isfile(args.infile):
        files = [args.infile]
    else:
        files = get_final_files(args.infile)

    print("Found final files:", *files)

    for file in files:
        mk_mising_dir(args.outdir)
        outfile = os.path.join(args.outdir, file.split("/")[-1][:-len("_final.csv")] ) + "_simulation.csv"
        print("New simulation will be saved in:", outfile)
        parameter_set = get_final_parameter_set(file)

        for k in overwrite_params.keys():
            parameter_set[k] = overwrite_params[k]

        if args.genealogy[0] == "s":
            tree = False
        else:
            tree = True

        dataset = pd.DataFrame()
        for i in np.arange(args.lanes):
            cells_simulated = simulate_cells(dt, n_cells, parameter_set, div_mode,
                            division_log_length, 
                            division_time, 
                            division_addition, tree=tree)[args.discard_cells+1:]
            temp_dataset = build_data_set_fixed_dt(cells_simulated, parameter_set['var_x'], parameter_set['var_g'], dt_measument, atol=1e-4)
            temp_dataset['lane'] = i+1
            # add the last "lane" to the data frame
            dataset = dataset.append(temp_dataset)

        # finally save the data
        dataset.to_csv( outfile )

if __name__ == "__main__":
    main()