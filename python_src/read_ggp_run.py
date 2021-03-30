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



# ============== get file names ============== #
def default_filebase(directory, sample, ext="_out"):
    return os.path.join(directory, sample+ext, sample) 

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

def get_minimization_file(filebase, paramter_settings):
    minimization_iter_file = filebase + get_param_code(paramter_settings) + ".csv"
    minimization_final_file = filebase + get_param_code(paramter_settings) + "_final.csv"
    return minimization_iter_file, minimization_final_file

def get_scan_files(filebase, paramter_settings):
    base = filebase + "_scan_"
    files = []
    for i,k in enumerate(paramter_settings.keys()):
        temp = base + k + ".csv"
        if os.path.isfile(temp):
            files.append(temp)
    return files

def get_prediction_files(filebase):
    return (filebase + '_prediction_forward.csv',
            filebase + '_prediction_backward.csv',
            filebase + '_prediction.csv')

def get_data_file(directory, sample):
    return os.path.join(directory, sample) + '.csv' 


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

def read_minimization(filename, last_n=None):
    df = pd.read_csv(filename, skiprows=14)    
    if last_n != None:
        return df.iloc[[-last_n]]
    return df

# ============== PLOTTING ============== #
# some basic plots for a ggp run

def plot_1dscans(filenames, plot_file, cols=3, width=14, l_col='log_likelihood'):
    """ plots the scan for all filenames as a grid """
    rows = np.ceil(len(filenames)/cols).astype(int)
    fig, axes = plt.subplots(rows, cols, figsize=(width,0.7*width/cols*rows))
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
    """ plots log likelihood vs number of iterations """

    fig, ax = plt.subplots(figsize=(10,7))
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



def plot_errors(errors_df, final_params, plot_file):
    """ plot errors for all epsilons of the hessian matrix given in file """
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



def compare_init_final(filename, plot_file, except_param=[]):
    """ bar plot of relative deviation between init and final parameter value """
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
    fig, ax = plt.subplots(figsize=(7,5))
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
    custom_legend = [Patch(facecolor='tab:orange', label='variances'),
                    Patch(facecolor='tab:blue', label='others')]
    # ax.legend(handles=custom_legend, bbox_to_anchor=(1,1), loc="upper left")
    ax.set_ylim([-1,1 ])

    plt.tight_layout()
    if plot_file != None:
        plt.savefig(plot_file + '_estim_params.pdf')

    plt.show()
    return comparison


