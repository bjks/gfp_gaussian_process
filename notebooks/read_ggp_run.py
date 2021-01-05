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



# ============== FILES ============== #
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
    
def read_1dscan(filename):
    tag = "scan_"
    base = filename.split("/")[-1]
    parameter = base[base.find(tag)+len(tag):-4]

    df = pd.read_csv(filename, skiprows=14)    
    return df[[parameter, 'likelihood']], parameter

def read_iteration_process(filename):
    df = pd.read_csv(filename, skiprows=14)    
    return df[["iteration", 'likelihood']]

def read_minimization(filename, last_n=None):
    df = pd.read_csv(filename, skiprows=14)    
    if last_n != None:
        return df.iloc[[-last_n]]
    return df

# ============== PLOTTING ============== #

def plot_1dscans(filenames, plot_file, cols=3, width=14):
    """ plots the scan for all filenames as a grid """
    rows = np.ceil(len(filenames)/cols).astype(int)
    fig, axes = plt.subplots(rows, cols, figsize=(width,0.7*width/cols*rows))
    # fig = plt.figure()
    for i, ax in enumerate(axes.ravel()):
        if i<len(filenames):
            scan, parameter = read_1dscan(filenames[i])

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


def plot_minimization(filename, plot_file):
    """ plots log likelihood vs number of iterations """
    mini_data = read_iteration_process(filename)
    fig, ax = plt.subplots(figsize=(7,5))
    iterations = mini_data.to_numpy()[:,0].astype(int)
    ll = mini_data.to_numpy()[:,1]
    ax.plot(iterations, ll)
    ax.scatter(iterations[-1], ll[-1], label="iteration: {:d}, log likelihood: {:.2f}".format(iterations[-1], ll[-1]))
    ax.legend()
    ax.set_xlabel("iterations")
    ax.set_ylabel("log likelihood")
    if plot_file != None:
        plt.savefig(plot_file + '_minimization.pdf')
    plt.show()



def compare_predictions(predictions, labels, col, data, data_col, data_slice, title=None, ratio=None, plot_file=None):
    """compares predictions

    Args:
        predictions (list of pandas DataFrames): predictions
        labels (list of stings): labels of predictios for legend
        col (sting): column in predictions DataFrames that will be plotted
        data (pandas DataFrame): data that will be plotted as dots 
        data_col (string): column in data DataFrame that will be plotted
        data_slice (slice): show only those data points
        title: Defaults to None.
        ratio (tuple of ints, optional): indices of predictions of which the ratio will be plotted on second axis. Defaults to None.
        plot_file (string, optional): output file. Defaults to None.
    """
    cmap = plt.cm.tab10
    fig, ax = plt.subplots(figsize=(15,7))
    plots = []

    if title !=None:
        plt.title(title)
        ax.set_ylabel(title)

    time = data['time'][data_slice]
    for i, p in enumerate(predictions):
        if i == 0:
            ls = '-'
        else:
            ls ='--'
        plots.append(ax.plot(time, p[col][data_slice], ls, c=cmap(i), label=labels[i])[0])
    ax.set_xlabel("time")
    plots.append(ax.plot(time, data[data_col][data_slice], 'o', c=cmap(0), label="data")[0])

    if ratio!=None:
        ax2 = ax.twinx()
        plots.append(ax2.plot(time, 
                    predictions[ratio[0]][col][data_slice]/predictions[ratio[1]][col][data_slice],
                    color='grey', label='ratio')[0])
        ax2.axhline(y=1, ls='--', color='grey')
        ax2.set_ylabel("ratio")

    ax.legend(plots, [l.get_label() for l in plots])
    if plot_file != None:
        plt.savefig(plot_file + '_' + data_col + '_prediction.pdf')
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



def compare_init_final(filename, plot_file):
    """ bar plot of relative deviation between init and final parameter value """
    init = read_params_config(filename)[["name", "init"]]
    init = init.set_index("name")

    minimized = read_minimization(filename, 1).transpose()
    reldev = [(init.loc[key].values[0]-minimized.loc[key].values[0])/init.loc[key].values[0] if init.loc[key].values[0]!=0 else None for key in init.index]

    # Table
    comparison = pd.DataFrame({'parameter': init.index, 
                                'simulation': [init.loc[key].values[0] for key in init.index],
                                'minimization': [minimized.loc[key].values[0] for key in init.index],
                                'relative deviation': reldev})

    # plot
    colors = ['tab:orange' if 1 else 'tab:blue' for key in init.index]
    fig, ax = plt.subplots(figsize=(7,5))
    plt.title(r"Relative deviation to initial value $(\Theta_i - \Theta_{i, opt})/\Theta_i $")

    x = [xi for xi, ri in zip(np.arange(len(reldev)), reldev) if ri != None]
    reldev = [ri for ri in reldev if ri != None]
    ax.bar(x, reldev, color=colors)

    for i, r in enumerate(reldev):
        if r<0:
            ax.text(x[i] -0.3, r -0.03 , "{:.3f}".format(r), va='center')
        else:
            ax.text(x[i] -0.3, r +0.03 , "{:.3f}".format(r), va='center')

    plt.xticks(np.arange(len(reldev)), init.index[x], rotation=45)
    custom_legend = [Patch(facecolor='tab:orange', label='variances'),
                    Patch(facecolor='tab:blue', label='others')]
    # ax.legend(handles=custom_legend, bbox_to_anchor=(1,1), loc="upper left")
    ax.set_ylim([-1,1 ])

    plt.tight_layout()
    if plot_file != None:
        plt.savefig(plot_file + '_estim_params.pdf')

    plt.show()
    return comparison



def compare_time_series(data, cols, data2, cols2, data_slice, title=None):
    cmap = plt.cm.tab10
    fig, ax = plt.subplots(figsize=(15,7))

    if title !=None:
        plt.title(title)

    time = data['time_min'][data_slice]
    plots = []

    for i, col in enumerate(cols):
        d = data[col][data_slice]
        plots.append(ax.plot(time, d, 'o', c=cmap(i), label=col)[0])

    for i, col in enumerate(cols2):
        d = data2[col][data_slice]
        plots.append(ax.plot(time, d, '-', c=cmap(i), label=col)[0])

    ax.set_xlabel('time (min)')
    ax.legend(plots, [l.get_label() for l in plots])
    plt.show()
