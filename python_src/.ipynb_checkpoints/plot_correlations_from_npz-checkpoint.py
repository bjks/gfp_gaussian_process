from correlation_from_joint import *
from read_ggp_run import *
import re
import numpy as np 
import matplotlib.pyplot as plt
import argparse
import copy
import os
import matplotlib.colors as mcolors
from matplotlib import cm
import pandas as pd

import scipy.stats
from scipy.optimize import curve_fit
import scipy.optimize as optimize


SMALL_SIZE = 16
MEDIUM_SIZE = 18

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

marker_by_condition = {"glycerol": 'o', "glucose": '^', "glucoseaa": 's', "acetate": 'p'}

color_by_promoter = {'hi1': 'tab:blue', 'hi3': 'tab:orange', 'med2': 'tab:green', 'med3': 'tab:red', 
                    'rrnB': 'tab:purple', 'rpsB': 'tab:brown', 'rpmB': 'tab:pink', 'rplN': 'tab:olive'}

color_by_condition = {"acetate": 'tab:blue',
                    "glycerol": 'tab:orange', 
                    "glucose": 'tab:green', 
                    "glucoseaa": 'tab:red'}


labels = {"l(t+dt)": "$\lambda(t+dt)$", 
            "l(t)": "$\lambda(t)$", 
            "q(t+dt)": "$q(t+dt)$", 
            "q(t)": "$q(t)$"}

def flips_zn_zm(x):
    return np.roll(x, 4, axis=(1,2))

class Correlation_function:
    '''
    Represents a Correlation function

    Constructed from a list of correlation instances

    Attributes:
        dt (np.array(N)):  List of time lags
        n (np.array(N)): List of number of joints
        corr_naive (np.array((8,8,N))): List of naive corrlation matrices
        corr_naive (np.array((8,8,N))): List of naive correlation matrices
        corr_mle (np.array((8,8,N))): List of mle correlation matrices
        corr_mle_err (np.array((8,8,N))): List of errors of mle 
    '''

    mapping = {z:i for i, z in enumerate(["x(t+dt)", "g(t+dt)", "l(t+dt)", "q(t+dt)", "x(t)", "g(t)", "l(t)", "q(t)"])}

    def __init__(self, correlations):
        dt = []
        n = []

        corr_naive = []
        corr_mle = []
        corr_mle_err = []
        cov_mle = []
        cov_mle_err = []

        corr_concentration_naive =[]
        corr_concentration_mle = []
        corr_concentration_mle_err = []
        cov_concentration_mle = []
        cov_concentration_mle_err = []

        for c in correlations:
            dt.append(c.dt)
            n.append(c.n)

            # zz
            corr_naive.append(c.corr_naive)
            corr_mle.append(c.corr_mle)
            corr_mle_err.append(c.corr_mle_err)

            cov_mle.append(c.cov_mle)
            cov_mle_err.append(c.cov_mle_err)

            #cc
            corr_concentration_naive.append(c.corr_concentration_naive)
            corr_concentration_mle.append(c.corr_concentration_mle)
            corr_concentration_mle_err.append(c.corr_concentration_mle_err)

            cov_concentration_mle.append(c.cov_concentration_mle)
            cov_concentration_mle_err.append(c.cov_concentration_mle_err)

        # save as numpy arrays
        self.dt = np.array(dt).astype(float)
        self.n = np.array(n)

        
        # add negative dt direction       
        ## z_n
        self.dt = np.append(-self.dt[::-1], self.dt)
        self.n = np.append(self.n[::-1], self.n)
        
        corr_naive  = np.concatenate((flips_zn_zm(corr_naive[::-1]), corr_naive), axis=0)
                
        corr_mle  = np.concatenate((flips_zn_zm(corr_mle[::-1]), corr_mle), axis=0)
        corr_mle_err  = np.concatenate((flips_zn_zm(corr_mle_err[::-1]), corr_mle_err), axis=0)
      
        cov_mle  = np.concatenate((flips_zn_zm(cov_mle[::-1]), cov_mle), axis=0)
        cov_mle_err  = np.concatenate((flips_zn_zm(cov_mle_err[::-1]), cov_mle_err), axis=0)
        
        

        ## concetration
        corr_concentration_naive  = np.concatenate((flips_zn_zm(corr_concentration_naive[::-1]),
                                                    corr_concentration_naive), axis=0)
                
        corr_concentration_mle  = np.concatenate((flips_zn_zm(corr_concentration_mle[::-1]), 
                                                  corr_concentration_mle), axis=0)
        corr_concentration_mle_err  = np.concatenate((flips_zn_zm(corr_concentration_mle_err[::-1]), 
                                                      corr_concentration_mle_err), axis=0)
      
        cov_concentration_mle  = np.concatenate((flips_zn_zm(cov_concentration_mle[::-1]), 
                                                 cov_concentration_mle), axis=0)
        cov_concentration_mle_err  = np.concatenate((flips_zn_zm(cov_concentration_mle_err[::-1]), 
                                                     cov_concentration_mle_err), axis=0)
        
        
        ## Cor(z,z)
        self.corr_naive = np.stack(corr_naive, axis=2)
        
        self.corr_mle = np.stack(corr_mle, axis=2)
        self.corr_mle_err = np.stack(corr_mle_err, axis=2)
        
        self.cov_mle = np.stack(cov_mle, axis=2)
        self.cov_mle_err = np.stack(cov_mle_err, axis=2)
    
        
        # Corr(c,c)
        self.corr_concentration_naive = np.stack(corr_concentration_naive, axis=2)

        self.corr_concentration_mle = np.stack(corr_concentration_mle, axis=2)
        self.corr_concentration_mle_err = np.stack(corr_concentration_mle_err, axis=2)

        self.cov_concentration_mle = np.stack(cov_concentration_mle, axis=2)
        self.cov_concentration_mle_err = np.stack(cov_concentration_mle_err, axis=2)


    def filter_by_n(self, n=0):
        filter_n = np.where(self.n>=n, True, False)
        
        self.dt = self.dt[filter_n]
        self.n = self.n[filter_n]

        # zz
        self.corr_naive = self.corr_naive[:, :, filter_n]
       
        self.corr_mle = self.corr_mle[:, :, filter_n]
        self.corr_mle_err = self.corr_mle_err[:, :, filter_n]

        self.cov_mle = self.cov_mle[:, :, filter_n]
        self.cov_mle_err = self.cov_mle_err[:, :, filter_n]

        
        # cc
        self.corr_concentration_naive = self.corr_concentration_naive[:, :, filter_n]
        
        self.corr_concentration_mle = self.corr_concentration_mle[:, :, filter_n]
        self.corr_concentration_mle_err = self.corr_concentration_mle_err[:, :, filter_n]

        self.cov_concentration_mle = self.cov_concentration_mle[:, :, filter_n]
        self.cov_concentration_mle_err = self.cov_concentration_mle_err[:, :, filter_n]


    def get_corr_naive(self, zi, zj):
        if zi=="c(t+dt)" and zj=="c(t)":
            return self.corr_concentration_naive[0, 1]
        return self.corr_naive[self.mapping[zi], self.mapping[zj]]

    def get_corr_mle(self, zi, zj):
        if zi=="c(t+dt)" and zj=="c(t)":
            return self.corr_concentration_mle[0, 1]
        return self.corr_mle[self.mapping[zi], self.mapping[zj]]

    def get_corr_mle_err(self, zi, zj):
        if zi=="c(t+dt)" and zj=="c(t)":
            return self.corr_concentration_mle_err[0, 1]
        return self.corr_mle_err[self.mapping[zi], self.mapping[zj]]

    def get_cov_mle(self, zi, zj):
        if zi=="c(t+dt)" and zj=="c(t)":
            return self.cov_concentration_mle[0, 1]
        return self.cov_mle[self.mapping[zi], self.mapping[zj]]

    def get_cov_mle_err(self, zi, zj):
        if zi=="c(t+dt)" and zj=="c(t)":
            return self.cov_concentration_mle_err[0, 1]
        return self.cov_mle_err[self.mapping[zi], self.mapping[zj]]
    
    
    
### PLOTS

def plot_correlation_tiles(corr_func, 
                            gamma_lambda=None, 
                            gamma_q=None, 
                            plot_file=None, 
                            xlim=[0, None], ylim=[-0.4, 1], 
                            min_joint_number=0,
                            title=None,
                            cov=False):


    if cov:
        get_func=Correlation_function.get_cov_mle
        get_err_func=Correlation_function.get_cov_mle_err
    else:
        get_func=Correlation_function.get_corr_mle
        get_err_func=Correlation_function.get_corr_mle_err
        
        
    fig, _ = plt.subplots(figsize=(8,8))

    if title!=None:
        fig.suptitle(title.replace('_', ' '), fontweight="bold")
    
    # ax0 = plt.subplot(321) # the figure has 3 row, 2 columns, and this plot is the first plot. 
    # ax1 = plt.subplot(322, sharex=ax0, sharey=ax0)

    # ax2 = plt.subplot(323, sharex=ax0, sharey=ax0)
    # ax3 = plt.subplot(324, sharex=ax0, sharey=ax0)
    # ax4 = plt.subplot(313) 

    ax0 = plt.subplot2grid((3, 2), (0, 0)) # the figure has 3 row, 2 columns, and this plot is the first plot. 
    ax1 = plt.subplot2grid((3, 2), (0, 1), sharex=ax0, sharey=ax0)
    ax2 = plt.subplot2grid((3, 2), (1, 0), sharex=ax0, sharey=ax0)
    ax3 = plt.subplot2grid((3, 2), (1, 1), sharex=ax0, sharey=ax0)
    ax4 = plt.subplot2grid((3, 2), (2, 0), colspan=2) 
    fig.tight_layout(h_pad=4)

    for ax in [ax0, ax1, ax2, ax3, ax4]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # set the x-spine
        ax.spines['left'].set_position('zero')

        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()


    for ax in [ax0, ax1, ax2, ax3]:
        # set the y-spine
        ax.spines['bottom'].set_position('zero')

        # turn off the top spine/ticks
        ax.spines['top'].set_color('none')
        ax.xaxis.tick_bottom()
        ax.tick_params(axis='x', pad=14)
        # remove first tick
        ax.xaxis.get_major_ticks()[0].set_visible(False)


    cf = copy.deepcopy(corr_func)
    cf.filter_by_n(min_joint_number)
    ######### Actual plots #########
    # ll
    ax0.errorbar(cf.dt, 
                get_func(cf,"l(t+dt)", "l(t)"), 
                yerr=get_err_func(cf,"l(t+dt)", "l(t)"),
                color="tab:blue")

    ax0.set_title(r"$\lambda(t+dt),\lambda(t)$", y=1, pad=-14)
    if gamma_lambda!=None:
        ax0.plot(cf.dt, np.exp(-gamma_lambda*cf.dt), color='grey', ls='--')

    # qq
    ax1.errorbar(cf.dt, 
                get_func(cf,"q(t+dt)", "q(t)"), 
                yerr=get_err_func(cf,"q(t+dt)", "q(t)"),
                color="tab:blue")

    ax1.set_title(r"$q(t+dt),q(t)$", y=1, pad=-14)
    if gamma_q!=None:
        ax1.plot(cf.dt, np.exp(-gamma_q*cf.dt), color='grey', ls='--')

    #ql
    ax2.errorbar(cf.dt, 
                get_func(cf,"q(t+dt)", "l(t)"), 
                yerr=get_err_func(cf,"q(t+dt)", "l(t)"),
                color="tab:orange")
    ax2.set_title(r"$q(t+dt),\lambda(t)$", y=1, pad=-14)

    #lq
    ax3.errorbar(cf.dt, 
                get_func(cf,"l(t+dt)", "q(t)"), 
                yerr=get_err_func(cf,"l(t+dt)", "q(t)"),
                color="tab:orange")
    ax3.set_title(r"$\lambda(t+dt),q(t)$", y=1, pad=-14)
    

    ax4.plot(cf.dt, cf.n, color="darkgrey")
    ax4.set_title("Number of joints", y=1.0, pad=-14)
    ax2.set_xlabel(r"time lag $dt$ (min)")
    ax3.set_xlabel(r"time lag $dt$ (min)")
    ax4.set_xlabel(r"time lag $dt$ (min)")

    if xlim[1]!=None:
        ns = [n for i,n in enumerate(cf.n) if cf.dt[i]<xlim[1] and cf.dt[i]>xlim[0] ]
        ax4.set_ylim([np.min(ns)*0.99, np.max(ns) *1.01 ])

    for ax in [ax0, ax1, ax2, ax3, ax4]:
        ax.set_xlim(xlim)
    for ax in [ax0, ax1, ax2, ax3]:
        ax.set_ylim(ylim)

    fig.tight_layout(h_pad=4)
    if plot_file != None:
        print("Saved in", plot_file)
        fig.savefig(plot_file, dpi=300, facecolor="white")
    else:
        plt.show()
    
    
def plot_xy_correlation(corr_func, x, y, 
                        gamma=None, 
                        mean_lambda=None, 
                        plot_file=None, 
                        title=None, 
                        xlim=[0, None], 
                        ylim=[None, None], 
                        min_joint_number=0):

    get_func=Correlation_function.get_corr_mle
    get_err_func=Correlation_function.get_corr_mle_err
    get_naive_func=Correlation_function.get_corr_naive

        
    correlation_function = copy.deepcopy(corr_func)

    correlation_function.filter_by_n(min_joint_number)
    dts = correlation_function.dt
    rs = get_func(correlation_function, x,y)
    errs = get_err_func(correlation_function, x,y)
    corr_naive = get_naive_func(correlation_function, x,y)



    # =========== plot =========== #
    # prepare axes depending on whether exponential is plotted
    if gamma!=None:
        fig, _ = plt.subplots(figsize=(8,8))
        ax0 = plt.subplot2grid((2, 1), (0, 0)) 
        ax1 = plt.subplot2grid((2, 1), (1, 0), sharex=ax0)
        axs = [ax0, ax1]
    else:
        fig, _ = plt.subplots(figsize=(8,4))

        ax0 = plt.subplot2grid((1,1), (0, 0)) 
        axs = [ax0]

    for ax in axs:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # set the x-spine
#         ax.spines['left'].set_position('zero')

        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()

        
    # add top x axis
    if mean_lambda != None:
        secax0 = ax0.secondary_xaxis('top', functions=(lambda x:  x/(np.log(2)/mean_lambda), 
                                                        lambda x: x*(np.log(2)/mean_lambda )))
        secax0.set_xlabel(r'norm. time lag $dt/$(mean doubling time)')
        ax0.axvline(np.log(2)/mean_lambda, color="grey", ls='--')

    if title!=None:
        ax0.set_title(title.replace('_', ' '),fontweight="bold")
        
    dts = np.array(dts).astype(float)

    x_label = x.replace("l", "\lambda")
    y_label = y.replace("l", "\lambda")

    ax0.errorbar(dts, rs, yerr=errs, label=r"mle $\langle {:s}, {:s}\rangle$".format(x_label,y_label), lw=2, color="tab:blue")
    ax0.plot(dts, corr_naive, label=r"naive $\langle {:s}, {:s}\rangle$".format(x_label,y_label), lw=2, color="tab:green")

    ax0.set_ylabel("correlation")
    ax0.legend()
    
    if gamma!=None:
        ax0.plot(dts, np.exp(-dts*gamma), label="exponential", color="darkgrey")
        ax1.errorbar(dts, rs-np.exp(-dts*gamma), yerr=errs, label=r"mle $\langle {:s}, {:s}\rangle$".format(x_label,y_label), lw=2, color="tab:blue")
        ax1.plot(dts, corr_naive-np.exp(-dts*gamma), label=r"naive $\langle {:s}, {:s}\rangle$".format(x_label,y_label), lw=2, color="tab:green")

        ax1.plot(dts, dts*0, color="darkgrey")
        if mean_lambda != None:
            ax1.axvline(np.log(2)/mean_lambda, color="grey", ls='--')

        ax1.set_ylabel("deviation from exponential")
        ax1.legend()
        
        ax1.set_xlabel("time lag $dt$")

    else:
        ax0.plot(dts, dts*0, color="darkgrey")
        ax0.set_xlabel("time lag $dt$")

    for ax in axs:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    if plot_file != None:
        print("Saved in", plot_file)
        fig.savefig(plot_file, dpi=300, facecolor="white")
    else:
        plt.show()

    return


def plot_xy_correlation_list(correlation_function_list, x, y, 
                                plot_file=None, 
                                label=[], 
                                gamma=[], 
                                mean_lambda = [],
                                scale_t=False, 
                                normalize=[],
                                xlim=[None, None],
                                ylim=[None,None], 
                                ylabel=None,
                                err_style="bar", 
                                log=False, 
                                min_joint_number=0, 
                                fit=False, 
                                highlight_x0=False,
                                highlight_y0=False,
                                cov=False, 
                                 color_by="condition",
                                legend=True):

    if cov:
        get_func=Correlation_function.get_cov_mle
        get_err_func=Correlation_function.get_cov_mle_err
    else:
        get_func=Correlation_function.get_corr_mle
        get_err_func=Correlation_function.get_corr_mle_err
        
        
    
    color_norm = mcolors.Normalize(vmin=0, vmax=9.9)

    labels = copy.deepcopy(label)
    gammas = copy.deepcopy(gamma)
    mean_lambdas = copy.deepcopy(mean_lambda)

    # =========== figure =========== #
    if fit:
        fig, ax0 = plt.subplots(figsize=(8.3, 4.8))

    else:
        fig = plt.figure(figsize=(8.3, 4.8))
        ax0 = plt.axes()
    axs = [ax0]

    for ax in axs:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # set the x-spine
#         ax.spines['left'].set_position('zero')

        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()

    # "pad" the following lists with Nones
    labels += [None] * len(correlation_function_list)
    gammas = list(gammas)+ [None] * len(correlation_function_list)
    mean_lambdas = list(mean_lambdas) + [None] * len(correlation_function_list)
    normalize = list(normalize) + [None] * len(correlation_function_list)

    new_ylim = [np.inf,0]
    # =========== plot =========== #
    for i, cf in enumerate(correlation_function_list):
        condition = labels[i].split('_')[0]
      
        if color_by=="condition":
            color = color_by_condition[condition]
        else:
            color = cm.tab10(color_norm(i)) 
        
        correlation_function = copy.deepcopy(cf)
        correlation_function.filter_by_n(min_joint_number)
        label = labels[i].replace('_', ' ')
        label = label.replace('glucoseaa', 'glucose+aa')
        gamma = gammas[i]
        mean_lambda = mean_lambdas[i]
        norm = normalize[i]
        

        # =========== correlation =========== #
        dts = correlation_function.dt
        rs = get_func(correlation_function, x,y)
        errs = get_err_func(correlation_function, x,y)
        if norm != None:
            rs /= norm
            errs /= norm

        ### Scale dt by mean doubling time 
        if scale_t:
            dts /= (np.log(2.)/mean_lambda)
            if gamma!=None:
                gamma*=(np.log(2)/mean_lambda)
        else:
            if mean_lambda!=None:
                ax0.axvline(np.log(2)/mean_lambda, color=color, alpha=0.5)


#         if fit:
#             def exponentials(t, *p):
#                 # return a1*np.exp(-t*gamma) + a2*np.exp(-t*b2)
#                 return p[0] * np.exp(-dts*p[1])

#             # ax0.plot(dts, exponentials(dts, *popt), color=color, ls='--', alpha=0.6)
#             filter = np.where(dts>0.5, True, False) *  np.where(dts<1.5, True, False)
#             def my_error(p):
#                 yfit = exponentials(dts, *p)[filter]
#                 return (yfit-rs[filter])**2/(errs[filter])**2

#             popt, pcov = optimize.leastsq(my_error,x0=[1, gamma/2.])
#             ax0.plot(dts[filter], exponentials(dts, *popt)[filter], color=color, ls='--', alpha=1)
#             label += r" $\tau={:.2f}$".format(1/gamma)
#             label += r" $\tau^\prime={:.2f}$".format(1/popt[1])
#             # ax0.plot(dts[dts<0.5], np.exp(-gamma*dts)[dts<0.5], color=color, ls='--', alpha=1)

#         else:
        if gamma!=None:
            r0 = rs[dts==0][0]
            ax0.plot(dts, r0*np.exp(-dts*gamma), ls='--', color=color, alpha=1)

        if xlim!=[None,None]:
            dt_filter = (dts<=xlim[-1])*(dts>=xlim[0])
            dts=dts[dt_filter]
            rs=rs[dt_filter]
            errs=errs[dt_filter]
        
        ### correlation plot with \without error bars
        if err_style == "bar":
            ax0.errorbar(dts, rs, yerr=errs, lw=2, color=color, label=label)
        elif err_style == "fill":
            ax0.plot(dts, rs, lw=2, color=color, label=label)
            ax0.fill_between(dts, rs-errs, rs+errs, color=color, alpha=0.4)
        elif err_style == None:
            ax0.plot(dts, rs, lw=2, color=color, label=label)
        else:
            ax0.plot(dts, dts*0, color=color)  
            
    if highlight_x0:
        ax0.axhline(y=0, ls='--', color="grey")
    if highlight_y0:
        ax0.axvline(x=0, ls='--', color="grey")
    # ===== layout ===== #
    if ylim != [None, None]:
        new_ylim = ylim
    else:
        new_ylim = [None, None]
#         new_ylim[0] = np.min([np.min(rs), new_ylim[0]])
#         new_ylim[1] = np.max([np.max(rs), new_ylim[1]])

    if log:
        ax0.set_yscale('log')
        
#         yticks = np.ravel([np.array([1, 5])*(10.**i) for i in np.arange(-10,10)])
#         yticks = yticks[(yticks>=new_ylim[0]) * (yticks<=new_ylim[1])]
        
#         ax0.set_yticks(yticks)
#         ax0.set_yticklabels(yticks)
    ax0.set_ylim(new_ylim)
        
    if ylabel == None:
        ax0.set_ylabel(r"$\langle {:s}, {:s}\rangle$".format(x,y))
    else:
        ax0.set_ylabel(ylabel)
    
    ax0.set_xlim(xlim)
#     if x[0]==y[0]:
#         ax0.legend(ncol=2, bbox_to_anchor=(1,1), loc="upper right")
#     else:
#         ax0.legend(ncol=2, bbox_to_anchor=(1,0.), loc="lower right")
    if legend:
        ax0.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=4)
#         ax0.legend(ncol=2, loc="lower left")

    if scale_t:
        ax0.set_xlabel(r'time lag $dt/$(mean doubling time)')

    else:
        ax0.set_xlabel(r"time lag $dt$")

    fig.tight_layout()
    if plot_file != None:
        print("Saved in", plot_file)
        fig.savefig(plot_file, dpi=300, facecolor="white", bbox_inches="tight")
    else:
        plt.show()
    plt.close()

    
def mean_quantity(cells, quantity="c"):
    coll = np.array([])
    for cell in cells:
        if quantity == "c":
            coll = np.append(coll, cell.gfp/np.exp(cell.log_length))
        if quantity == "g":
            coll = np.append(coll, cell.lt)
        if quantity == "q":
            coll = np.append(coll, cell.qt)  
    return np.mean(coll)


def get_filtered_keys(indict, filter=None, keep="both"):
    out = []
    for k in indict.keys():
        if filter==None or filter in k.split('_'):
            if keep=="both":
                out.append(k)
            elif keep=="condition":
                out.append(k.split('_')[0])
            elif keep=="promoter":
                out.append(k.split('_')[1])
            else:
                print("ERROR")
                return None
    return out

def get_filtered_values(indict, filter=None):
    out = []
    for k in indict.keys():
        if filter==None or filter in k.split('_'):
            out.append(indict[k])
    return np.array(out)

def sort_dict(indict, sort_key):
    out = {}
    for sk in sort_key:
        for k in indict.keys():
            if sk in k:
                out[k] = indict[k]
    return out
