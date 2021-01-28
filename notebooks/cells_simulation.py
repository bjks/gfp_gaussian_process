
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
import subprocess
import os

# ========================================== #
# =============== SIMULATION =============== #
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

    def to_df(self, n=1):
        return pd.DataFrame({   "cell_id": ([self.cell_id]*len(self.time))[::n],
                                "time_min": self.time[::n],
                                "parent_id": ([self.parent_id]*len(self.time))[::n],
                                "log_length": self.log_length[::n], 
                                "gfp": self.gfp[::n],
                                "lt": self.lt[::n],
                                "qt": self.qt[::n]})

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

            if lt == None:    lambda0 = None
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



# =============== Simulation functions ===============#

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

def cell_divsion(cell, var_dx, var_dg, no_cells):

    cell1 = Cell(np.random.normal(loc=cell.log_length[-1] - np.log(2), scale=np.sqrt(var_dx)),
                    np.random.normal(loc=cell.gfp[-1]/2, scale=np.sqrt(var_dg)), cell.lt[-1], cell.qt[-1],
                    time0 = cell.time[-1],
                    cell_id = no_cells + 1, parent_id=cell.cell_id)

    cell2 = Cell(np.random.normal(loc=cell.log_length[-1] - np.log(2), scale=np.sqrt(var_dx)),
                    np.random.normal(loc=cell.gfp[-1]/2, scale=np.sqrt(var_dg)), cell.lt[-1], cell.qt[-1], 
                    time0 = cell.time[-1],
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
        if cell.log_length[-1] - cell.log_length[0] > division_addition:
            return True

# ============= FULL SIMULATION ============= # 

def simulate_cells(dt, n_cells, parameter_set, div_mode, division_log_length=None, division_time=None, division_addition=None):
    gfp0 = 3*parameter_set['mean_q']/parameter_set['mean_lambda']
    cell_queue = [Cell(division_log_length-np.log(2), gfp0, parameter_set['mean_lambda'], parameter_set['mean_q'])]
    cells_simulated = []
    no_cells = 0   # total number of cells (in queue and calculated)

    while len(cells_simulated) < n_cells:
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
                                                            no_cells)

                # remove the simulated cell from queue and add the new ones 
                cell_queue.pop(cell_index)
                cell_queue.append(cell1)
                cell_queue.append(cell2)
                break
            else:
                pass
        progress_bar_n = np.around(len(cells_simulated)/n_cells*20).astype(int)
        progress_bar = '='*progress_bar_n + ' '*(20-progress_bar_n)

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

def get_next_file_name(out_dir, no=None):
    sample = out_dir.split('/')[-1]
    if no != None:
        new_dir = os.path.join(out_dir, sample+"_{:d}".format(no),'')
        new_file = os.path.join(new_dir, sample+"_{:d}".format(no)+".csv")
        return new_dir, new_file

    for i in range(1000):
        new_dir = os.path.join(out_dir, sample+"_{:d}".format(i),'')
        new_file = os.path.join(new_dir, sample+"_{:d}".format(i)+".csv")
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)
            return new_dir, new_file


def write_param_file(filename, parameters, non_default={}):
    with open(filename, "w") as fin:
        fin.write("# Generated config file for simulated data\n")
        for k, v in parameters.items():
            if k in non_default:
                if non_default[k][0] == "fixed":
                    fin.write("{:s} = {:.2E}\n".format(k, v))
                elif non_default[k][0] == "free":
                    fin.write("{:s} = {:.2E}, {:.2E}\n".format(k, v, v))
                else:
                    fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format(k,   v*non_default[k][0],    
                                                                                    v*non_default[k][1], 
                                                                                    v*non_default[k][2], 
                                                                                    v*non_default[k][3] ))
            else:
                if v==0:
                    fin.write("{:s} = {:.2E}, {:.2E}, 0, {:.2E}\n".format(k, v, 1e-3, 1e15)) 
                else:
                    fin.write("{:s} = {:.2E}, {:.2E}, {:.2E}, {:.2E}\n".format(k, v, v*1e-2, v*0.3, v*11. ))

def write_csv_config(filename):
    with open(filename, "w") as fin:
        fin.write("# Generated config file for simulated data\n")
        fin.write("time_col = time_min  \n")
        fin.write("length_col = log_length_noise \n")
        fin.write("length_islog = true \n")
        fin.write("fp_col = gfp_noise \n")
        fin.write("parent_tags = parent_id \n")
        fin.write("cell_tags = cell_id \n")
        fin.write("rescale_time = 1 \n")

def build_data_set(cells_simulated, var_x, var_g, n):
    print("Every", n, "th datapoint is saved")
    dataset = pd.DataFrame()
    for i in range(len(cells_simulated)):
        next_celldf = cells_simulated[i].to_df(n)
        next_celldf['log_length_noise'] = next_celldf['log_length'] + np.random.normal(loc=np.zeros_like( next_celldf['log_length']), scale=np.sqrt(var_x))
        next_celldf['gfp_noise'] = next_celldf['gfp'] + np.random.normal(loc=np.zeros_like( next_celldf['gfp']), 
                                                                        scale=np.sqrt(var_g))
        dataset = dataset.append(next_celldf)

    return dataset


# =============== RUN COMMAND =============== #
def suggest_run_command(directory, filename, modes ="-s -m -p"):
    cmd = "../bin/gfp_gaussian -c " + os.path.join(directory, "csv_config.txt") +  " -b " + os.path.join(directory, "parameters.txt") + " -r 1e-3  -i " + filename + " -l 0 "
    return cmd + modes


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