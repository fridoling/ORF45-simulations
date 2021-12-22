import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import pickle
import copy
import os
from expts import *
import re

params_file = "../../data/incell/params_incell.pickle"
with open(params_file, 'r') as f:
    model, popt, ens, cost = pickle.load(f)

times = np.hstack(([0], np.linspace(999, 1025, 101)))
t_iv = 1000

SUBtot = 0.01
trajs = {}
lower = {}
upper = {}
for net_id, net in nets_sub['INCELL_WT'].items():
    net.set_var_ic("SUBCtot", SUBtot)
    net.set_var_ic("SUBFtot", SUBtot)
    trajs[net_id] = net.integrate(times, popt)
    traj_set = Ensembles.ensemble_trajs(net, times, ens[::10])
    lower[net_id], upper[net_id] = Ensembles.traj_ensemble_quantiles(traj_set, (0.0, 1.0))

fontsize = 14
plot_vars = ["pSUBC", "pSUBF"]
for net_id in nets_sub["INCELL_WT"].keys():
    fig, ax = plt.subplots(1,1,figsize = (4,2.5))
    traj = trajs[net_id]
    time_traj = traj.timepoints[traj.timepoints>=999] - t_iv
    for var in plot_vars:
        var_traj = traj.get_var_traj(var)[traj.timepoints>=999]
        ax.plot(time_traj, var_traj/0.01, label = var, lw = 3)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_position(('data', -1))
    ax.set_xlim([-4,20])
    ax.set_xlabel("time (minutes)")
    ax.set_ylabel("pSUB/total")
    ax.set_ylim([-0.02,0.2])
    plt.tight_layout()    
    outfile = "../../res/Figure8B_substrates_"+net_id
    for ext in [".pdf", ".eps"]:
        plt.savefig(outfile+ext)
    plt.close("all")


for net_id in nets["INCELL_WT"].keys():
    fig, ax = plt.subplots(1,1,figsize = (4,2.5))
    traj = trajs[net_id]
    time_traj = traj.timepoints[traj.timepoints>=999] - t_iv
    for var in plot_vars:
        var_traj = traj.get_var_traj(var)[traj.timepoints>=999]
        ax.plot(time_traj, var_traj/0.01, label = var, lw = 3)
        lower_var = lower[net_id].get_var_traj(var)[lower[net_id].timepoints>=999]/0.01
        upper_var = upper[net_id].get_var_traj(var)[upper[net_id].timepoints>=999]/0.01
        time_traj_fill = lower[net_id].timepoints[lower[net_id].timepoints>=999] - t_iv
        ax.fill_between(time_traj_fill, lower_var, upper_var, alpha = 0.2)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_position(('data', -1))
    ax.set_xlim([-4,20])
    ax.set_xlabel("time (minutes)")
    ax.set_ylabel("pSUB/total")
    ax.set_ylim([-0.02,0.2])
    plt.tight_layout()    
    outfile = "../../res/Figure8B_substrates_with_ens_"+net_id
    for ext in [".pdf", ".svg"]:
        plt.savefig(outfile+ext)
