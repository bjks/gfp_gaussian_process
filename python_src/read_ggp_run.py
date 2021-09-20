
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl

import os

plt.rcParams["font.family"] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = ['CMU Sans Serif']
plt.rcParams['mathtext.default'] = 'regular'
params = {'text.usetex': False, 'mathtext.fontset': 'cm'}
plt.rcParams.update(params)

SMALL_SIZE = 14
MEDIUM_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title


LATEX_LABEL = {'mean_lambda': r'$\bar{\lambda}$',
            'gamma_lambda': r'$\gamma_\lambda$',
            'var_lambda': r'$\sigma^2_\lambda$',
            'mean_q': r'$\bar{q}$',
            'gamma_q': r'$\gamma_q$',
            'var_q': r'$\sigma^2_q$',
            'beta': r'$\beta$',
            'var_x': r'$\sigma^2_x$',
            'var_g': r'$\sigma^2_g$',
            'var_dx': r'$\sigma^2_{dx}$',
            'var_dg': r'$\sigma^2_{dg}$'} 

class GGP_cell:
    def __init__(self, cell_id = 0, parent_id=-1):
        self.parent_id = parent_id
        self.cell_id = cell_id
        self.log_length = []
        self.gfp = []
        self.time = []

        self.mean_x = []
        self.mean_g = []
        self.mean_l = []
        self.mean_q = []

        self.cov_xx = []
        self.cov_gg = []
        self.cov_ll = []
        self.cov_qq = []


def df2ggp_cells(dataset, 
            time="time", 
            log_length="log_length", gfp="fp", 
            mean_x="mean_x", mean_g="mean_g", 
            mean_l="mean_l", mean_q="mean_q",
            cov_xx="cov_xx",
            cov_gg="cov_gg",
            cov_ll="cov_ll",
            cov_qq="cov_qq",
            cell_id="cell_id", 
            parent_id="parent_id"):
    """ 
    dataset (pandas data frame as read from csv file) to list of GGP_cell instances, m
    written for ggp output
    """
    cell_list = []
    last_cell = ""
    for _, row in dataset.iterrows(): 
        if row[cell_id] != last_cell:
            new_cell = GGP_cell(
                        cell_id=row[cell_id], 
                        parent_id=row[parent_id])
            cell_list.append(new_cell)

        cell_list[-1].log_length.append(row[log_length])
        cell_list[-1].gfp.append(row[gfp])
        cell_list[-1].time.append(row[time])

        cell_list[-1].mean_x.append(row[mean_x])
        cell_list[-1].mean_g.append(row[mean_g])
        cell_list[-1].mean_l.append(row[mean_l])
        cell_list[-1].mean_q.append(row[mean_q])

        cell_list[-1].cov_xx.append(row[cov_xx])
        cell_list[-1].cov_gg.append(row[cov_gg])
        cell_list[-1].cov_ll.append(row[cov_ll])
        cell_list[-1].cov_qq.append(row[cov_qq])

        last_cell = row[cell_id]
    return cell_list

                              

# ============== get file names ============== #

# directory: directory in which the input file is
# sample: filename of the input file

# File base: example_dir/example_sample/example_sample
def default_filebase(directory, sample, ext="_out"):
    return os.path.join(directory, sample+ext, sample) 

def get_filebase(directory, sample, out):
    if out == None:
        return os.path.join(directory, sample + "_out", sample)
    return os.path.join(directory, out, sample)

# Files
def get_param_code(paramter_settings):
    param_code = '_f'
    for i,k in enumerate(paramter_settings.keys()):
        if paramter_settings[k]=='free':
            param_code+=str(i)

    param_code += "_b"
    for i,k in enumerate(paramter_settings.keys()):
        if paramter_settings[k]=='bound':
            param_code+=str(i)
    return param_code


def get_data_file(directory, sample):
    return os.path.join(directory, sample) + '.csv' 


def get_minimization_file(filebase, paramter_settings):
    minimization_iter_file = filebase + get_param_code(paramter_settings) + "_iterations.csv"
    minimization_final_file = filebase + get_param_code(paramter_settings) + "_final.csv"
    return minimization_iter_file, minimization_final_file


def get_prediction_files(filebase, paramter_settings=None):
    if paramter_settings==None:
        return (filebase + '_prediction_forward.csv',
            filebase + '_prediction_backward.csv',
            filebase + '_prediction.csv')

    param_code = ''
    for p in paramter_settings:
        param_code += get_param_code(p) 
    return (filebase + param_code + '_prediction_forward.csv',
            filebase + param_code + '_prediction_backward.csv',
            filebase + param_code + '_prediction.csv')


def get_scan_files(filebase, paramter_settings):
    base = filebase + "_scan_"
    files = []
    for _,k in enumerate(paramter_settings.keys()):
        temp = base + k + ".csv"
        if os.path.isfile(temp):
            files.append(temp)
    return files


def get_all_filenames(filebase, paramter_settings):
    minimization_iter_file, minimization_final_file = get_minimization_file(filebase, paramter_settings)
    scans = get_scan_files(filebase, paramter_settings)
    f, b, p = get_prediction_files(filebase, [paramter_settings])
    return {"iter": minimization_iter_file, 
            "final": minimization_final_file,
            "scans": scans,
            "prediction": p,
            "backward": b,
            "forward": f}

def fixed_params():
    return {'mean_lambda': 'fixed',
            'gamma_lambda': 'fixed',
            'var_lambda': 'fixed',
            'mean_q': 'fixed',
            'gamma_q': 'fixed',
            'var_q':'fixed',
            'beta':'fixed',
            'var_x':'fixed',
            'var_g':'fixed',
            'var_dx':'fixed',
            'var_dg':'fixed'} 

# ============== READING ============== #
def read_params_config(filename):
    df = pd.read_csv(filename, nrows=11)
    return df

def get_params_config(df, param):
    return df.loc[df['name'] == param]
    
def read_1dscan(filename, l='log_likelihood'):
    tag = "scan_"
    base = filename.split("/")[-1]
    parameter = base[base.find(tag)+len(tag):-4]

    df = pd.read_csv(filename, skiprows=14)    
    return df[[parameter, l]], parameter

def read_iteration_process(filename):
    df = pd.read_csv(filename, skiprows=14)    
    return df[["iteration", 'log_likelihood']]

def read_correlation(filename, n=12):
    return pd.read_csv(filename, skiprows=n)

def read_minimization(filename, last_n=0):
    df = pd.read_csv(filename, skiprows=14)    
    if last_n != 0:
        return df.iloc[[-last_n]]
    return df

# ============== PLOTTING ============== #
# some basic plots for a ggp run
def plot_1dscans(filenames, plot_file=None, cols=3, width=14, l_col='log_likelihood'):
    """ plots the scan for all filenames as a grid """
    rows = np.ceil(len(filenames)/cols).astype(int)
    _, axes = plt.subplots(rows, cols, figsize=(width,0.7*width/cols*rows))
    # fig = plt.figure()
    for i, ax in enumerate(axes.ravel()):
        if i<len(filenames):
            scan, parameter = read_1dscan(filenames[i], l_col)

            param_range = scan.to_numpy()[:,0]
            ll = scan.to_numpy()[:,1]

            ax.plot(param_range,ll, label='scan')

            params_config = read_params_config(filenames[i])
            init = get_params_config(params_config, parameter)["init"].values[0]
            init_idx = np.searchsorted(param_range, init)
            ax.axvline(x=init, label='init', ls='--', color='tab:orange')
            ax.axhline(y=ll[init_idx], ls='--', color='tab:orange')

            x_scan_max, y_scan_max = param_range[np.nanargmax(ll)], np.nanmax(ll)

            ax.axvline(x=x_scan_max, label='ll max', ls='--', color='tab:green')
            ax.axhline(y=y_scan_max, ls='--', color='tab:green')

            ax.set_xlabel(parameter)
            ax.set_ylabel('log likelihood')
            ax.ticklabel_format(style='sci', scilimits=(0,1), useOffset=False)
            ax.legend()
        else:
            plt.delaxes(ax)
    plt.tight_layout()
    if plot_file != None:
        plt.savefig(plot_file + '_1dscans.pdf')
    plt.show()


def plot_minimization(filenames, plot_file=None, labels=None, log=None):
    """ plots log likelihood vs number of iterations, needs iterations file """
    _, ax = plt.subplots(figsize=(10,7))
    for i, filename in enumerate(filenames):
        mini_data = read_iteration_process(filename)
        iterations = mini_data.to_numpy()[:,0].astype(int)
        ll = mini_data.to_numpy()[:,1]
        ax.plot(iterations, ll)
        if labels != None:
            label = labels[i]+", iteration: {:d}, log likelihood: {:.2f}".format(iterations[-1], ll[-1])
        else:
            label = "iteration: {:d}, log likelihood: {:.2f}".format(iterations[-1], ll[-1])
        ax.scatter(iterations[-1], ll[-1], label=label)
    
    ax.legend()
    ax.set_xlabel("iterations")
    ax.set_ylabel("log likelihood")
    if log == 'x':
        plt.xscale('log')
    if log == True:
        plt.xscale('log')
        plt.yscale('log')

    if plot_file != None:
        plt.savefig(plot_file + '_minimization.pdf')
    plt.show()



def plot_errors(filename, plot_file=None):
    """ plot errors for all epsilons of the hessian matrix given in file, needs final file """

    final_params = pd.read_csv(filename, nrows=11)['final']
    errors_df = pd.read_csv(filename, skiprows=14)

    mpl.style.use('default')
    # mpl.style.use('seaborn')

    fig, ax = plt.subplots(figsize=(7,5))
    for i, para in enumerate(errors_df.columns[1:]):
        ax.set_xlabel('relative $h_r= h / \\theta_i$ ')
        ax.set_ylabel('relative error')

        ax.scatter(errors_df['epsilon'].astype(str), errors_df[para]/final_params[i], label=para)
        ax.plot(errors_df['epsilon'].astype(str), errors_df[para]/final_params[i])

    fig.legend(bbox_to_anchor=[1, 0.5], loc=10)
    if plot_file != None:
        plt.savefig(plot_file + '_errors.pdf')
    plt.show()  



def compare_init_final(filename, plot_file=None, except_param=[]):
    """ bar plot of relative deviation between init and final parameter value, needs final file  """
    init = read_params_config(filename)[["name", "init"]]
    init = init.set_index("name")

    minimized = read_minimization(filename, 1).transpose()
    reldev = [(init.loc[key].values[0]-minimized.loc[key].values[0])/init.loc[key].values[0] if init.loc[key].values[0]!=0 and key not in except_param else None for key in init.index]

    # Table
    comparison = pd.DataFrame({'parameter': init.index, 
                                'simulation': [init.loc[key].values[0] for key in init.index],
                                'minimization': [minimized.loc[key].values[0] for key in init.index],
                                'relative deviation': reldev})

    # plot
    colors = ['tab:orange' if 1 else 'tab:blue' for key in init.index]
    _, ax = plt.subplots(figsize=(7,5))
    plt.title(r"Relative deviation to initial value $(\Theta_i - \Theta_{i, opt})/\Theta_i $")

    reldev = [ri for ri in reldev if ri != None]
    x = np.arange(len(reldev))
    ax.bar(x, reldev, color=colors)

    for i, r in enumerate(reldev):
        if r<0:
            ax.text(x[i] -0.3, r -0.03 , "{:.3f}".format(r), va='center')
        else:
            ax.text(x[i] -0.3, r +0.03 , "{:.3f}".format(r), va='center')

    xticks = [key for key in init.index if key not in except_param]

    plt.xticks(np.arange(len(reldev)), xticks, rotation=45)
    # custom_legend = [Patch(facecolor='tab:orange', label='variances'),
    #                 Patch(facecolor='tab:blue', label='others')]
    # ax.legend(handles=custom_legend, bbox_to_anchor=(1,1), loc="upper left")
    ax.set_ylim([-1,1 ])

    plt.tight_layout()
    if plot_file != None:
        plt.savefig(plot_file + '_estim_params.pdf')

    plt.show()
    return comparison
        

# ==================================================== #
# ==================================================== #
# "Advanced" analysis
def plot_noisy_param_run(filenames, params_config, skip=0, cols=3, width=14):
    """ needs list of final files """
    no_params = len(params_config)
    rows = np.ceil(no_params/cols).astype(int)
    _, axes = plt.subplots(rows, cols, figsize=(width,0.7*width/cols*rows))
    ax = axes.ravel()
    for i, k in enumerate(params_config.items()):
        n = 0
        for minimization_final_file in filenames:
            if os.path.isfile(minimization_final_file):
                parameters_settings = read_params_config(minimization_final_file)
                final = get_params_config(parameters_settings, k[0])["final"].values[0]
                init = get_params_config(parameters_settings, k[0])["init"].values[0]

                err_bar = float(pd.read_csv(minimization_final_file, skiprows=14)[k[0]].iloc[2])
                ax[i].errorbar(n, final, yerr=err_bar, color='tab:blue',  fmt='o', ms=3, 
                            label="deviation prediction")

                # ax[i].scatter(n, final, color="tab:blue")
                ax[i].scatter(n, init, color="tab:orange")
            else:
                print(minimization_final_file, "not found!")
            n += 1

        # ax[i].ticklabel_format(style='sci', scilimits=(0,1), useOffset=False)
        ax[i].set_ylabel(k[0])
        ax[i].set_xlabel("run")

    for i in range(len(params_config), len(ax)):
        plt.delaxes(ax[i])
    plt.tight_layout()

    plt.show()



# ==================================================== #
# Prediction #
# ==================================================== #
def plot_predictions(filename, start=None, stop=None, step=None, 
                    time_unit=("min", 60), skip_row=13, xlim=[None, None]):
    """ needs a prediction file, start, stop, step refers to cells """
    _, axes = plt.subplots(4, 1, figsize=(8,10))
    ax = axes.ravel()

    data = pd.read_csv(filename, skiprows=skip_row)
    cells_data = df2ggp_cells(data)[start: stop: step]
     

    norm = mpl.colors.Normalize(vmin=-len(cells_data)/2, vmax=len(cells_data))
    if len(cells_data)==1:
        norm = mpl.colors.Normalize(vmin=-10*len(cells_data), vmax=len(cells_data))

    cmap_data = mpl.cm.ScalarMappable(cmap='Oranges', norm=norm)
    cmap_data.set_array([])

    cmap_prediction = mpl.cm.ScalarMappable(cmap='Blues', norm=norm)
    cmap_prediction.set_array([])

    s = 4
    lw = 1
    for i, cell in enumerate(cells_data):
        time = np.array(cell.time) / time_unit[1]
        data_color = cmap_data.to_rgba(i)
        prediction_color = cmap_prediction.to_rgba(i)

        ax[0].scatter(time, cell.log_length, color=data_color, s=s)
        ax[0].plot(time, cell.mean_x, color=prediction_color, lw=lw)
        ax[0].fill_between(time, cell.mean_x-np.sqrt(cell.cov_xx), cell.mean_x+np.sqrt(cell.cov_xx), 
                    color=prediction_color, alpha=0.4)

        ax[1].scatter(time, cell.gfp, color=data_color, s=1)
        ax[1].plot(time, cell.mean_g, color=prediction_color, lw=lw)
        ax[1].fill_between(time, cell.mean_g-np.sqrt(cell.cov_gg), cell.mean_g+np.sqrt(cell.cov_gg), 
                    color=prediction_color, alpha=0.4)
        
        ax[2].plot(time, cell.mean_l, color=prediction_color, lw=lw)
        ax[2].fill_between(time, cell.mean_l-np.sqrt(cell.cov_ll), cell.mean_l+np.sqrt(cell.cov_ll), 
                    color=prediction_color, alpha=0.4)

        ax[3].plot(time, cell.mean_q, color=prediction_color, lw=lw)
        ax[3].fill_between(time, cell.mean_q-np.sqrt(cell.cov_qq), cell.mean_q+np.sqrt(cell.cov_qq), 
                    color=prediction_color, alpha=0.4)

    ax[0].set_ylabel("log lentgh")   
    ax[1].set_ylabel("YFP content (a.u.)")   
    ax[2].set_ylabel(r"growth rate $\lambda$")   
    ax[3].set_ylabel(r"production rate $q$")   


    ax[0].set_xlabel("time ({:s})".format(time_unit[0]))  
    ax[1].set_xlabel("time ({:s})".format(time_unit[0]))  
    ax[2].set_xlabel("time ({:s})".format(time_unit[0]))  
    ax[3].set_xlabel("time ({:s})".format(time_unit[0]))  

    for i in range(4):
        ax[i].set_xlim(xlim)

    for i in range(len(ax)):
        # ax[i].legend()
        ax[i].grid(True)
    plt.show()

# ==================================================== #

def plot_raw_data(filename, start=None, stop=None, step=None, time_unit=("min", 60), skip_row=13, scatter=True):
    """ needs a prediction file, start, stop, step refers to cells """
    _, axes = plt.subplots(2, 1, figsize=(8,5))
    ax = axes.ravel()

    data = pd.read_csv(filename, skiprows=skip_row)
    cells_data = df2ggp_cells(data)[start: stop: step]
     

    norm = mpl.colors.Normalize(vmin=-len(cells_data)/2, vmax=len(cells_data))
    if len(cells_data)==1:
        norm = mpl.colors.Normalize(vmin=-10*len(cells_data), vmax=len(cells_data))
    cmap_data = mpl.cm.ScalarMappable(cmap='Oranges', norm=norm)
    cmap_data.set_array([])

    cmap_prediction = mpl.cm.ScalarMappable(cmap='Blues', norm=norm)
    cmap_prediction.set_array([])

    for i, cell in enumerate(cells_data):
        time = np.array(cell.time) / time_unit[1]
        data_color = cmap_data.to_rgba(i)

        if scatter:
            ax[0].scatter(time, cell.log_length, color=data_color, s=10)
            ax[1].scatter(time, cell.gfp, color=data_color, s=10)
        else:
            ax[0].plot(time, cell.log_length, color=data_color, lw=0.2)
            ax[1].plot(time, cell.gfp, color=data_color, lw=0.2)



    ax[0].set_ylabel("log lentgh")   
    ax[1].set_ylabel("YFP content (a.u.)")   

    ax[0].set_xlabel("time ({:s})".format(time_unit[0]))  
    ax[1].set_xlabel("time ({:s})".format(time_unit[0]))  


    for i in range(len(ax)):
        # ax[i].legend()
        ax[i].grid(True)
    plt.show()
