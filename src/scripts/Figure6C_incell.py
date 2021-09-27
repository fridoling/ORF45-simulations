import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import pickle
import copy
import os
from expts import *

params_file = "../../data/incell/params_incell.pickle"
params_file_invitro = "../../data/invitro/params_invitro.pickle"

if not os.path.exists(params_file):
    exp_ids = ['INCELL_WT', 'INCELL_RSK_KO']
    
    with open(params_file_invitro, 'r') as f:
        model, popt_invitro, ens_invitro, cost_invitro = pickle.load(f)
        
    pars_fixed = KeyedList([])

    pars_constrained = popt_invitro

    pars_free = KeyedList([
        ('KD_RP2', 1.0),
        ('kon_RP2', 1.0),
        ('kdp_R', 1.0),
        ('kp_K_bg', 1.0),
        ('kdp_K', 1.0),
        ('kp_K_egf', 1.0)
    ])

    model, popt, ens, cost = functions.fit_exps(exps,
                                 nets,
                                 pars_fixed,
                                 pars_constrained,
                                 pars_free,
                                 exp_ids=exp_ids,
                                 global_fit=True,
                                 global_it=10000,
                                 return_ens=True)

    with open(params_file, 'w') as f:
        pickle.dump((model, popt, ens, cost), f)
else: 
    with open(params_file, 'r') as f:
        model, popt, ens, cost = pickle.load(f)
## Save parameters for table
table_pars_incell = ['KD_ER', 'KD_EpR', 'KD_OR', 'KD_OpR', 'KD_EO', 'KD_pEO', 'KD_EK', 'KD_EP', 'KD_RP2', 'koff_ER', 'koff_OR', 
        'koff_EO', 'koff_pEO', 'koff_EK', 'koff_EP', 'koff_RP2', 'kp_E', 'kp_R', 'kp_K_bg', 'kp_K_egf', 'kdp_E', 'kdp_R', 'kdp_K', 'a', 'd']
net_basic.set_var_vals(popt)
table_dict_incell = { par: net_basic.get_var_val(par) for par in table_pars_incell }

with open("../../data/incell/table_pars_incell.pickle", "wb") as f:
    pickle.dump(table_dict_incell, f)

## Ensemble analysis:
out_vars = ["KD_RP2", "koff_RP2", "kp_K_bg", "kp_K_egf", "kdp_K", "kdp_R"]
fig, ax = plt.subplots(figsize = (len(out_vars)*1.0, 3))
functions.plot_ens(ens, net_basic, out_vars, params_opt = popt, step = 10, file = None, axis = ax, mode = "std")
ax.legend(bbox_to_anchor = (1,0.5), loc = "center left")
ax.set_title("Ensemble fit for incell experiments", loc = "left");
plt.tight_layout()
plt.savefig("../../res/ens_fit_incell.eps")


   
# plot
exp_ids = ['INCELL_WT', 'INCELL_RSK_KO']
colors = {'K': 'b', 'ORF': 'r'}
fontsize = 14
nf = {'pEtot': 1/0.7, 'pRtot': 1/2.0}

t_iv = 1000
times = np.hstack(([0], np.linspace(999, 1025, 101)))

params = popt
model.params.update(params)
model.CalculateForAllDataPoints(params)
sf = model.compute_scale_factors(1)[0]


for exp_id in exp_ids:
    data = exps[exp_id].GetData()
    vars = list(np.unique([var for key in data.keys() for var in data[key].keys()]))
    
    for var in vars:
        fig, ax = plt.subplots(1, 1, figsize=(5,3), sharey=True)
        sf_i = sf[exp_id][var]
        ymax = 0
        xmax = 0
        for net_key, net in nets[exp_id].items():
            net.set_var_ics(params)
            net.set_var_vals(params) 
            net_id = net.id
            t_data = np.array(exps[exp_id].GetData()[net_id][var].keys())
            var_data = np.array([d[0] for d in exps[exp_id].GetData()[net_id][var].values()])
            traj = net.integrate(times)
            t_traj = traj.timepoints[traj.timepoints>=999] - t_iv
            var_traj = traj.get_var_traj(var)[traj.timepoints>=999]
            ax.scatter(np.sort(t_data)-t_iv, var_data[np.argsort(t_data)]/sf_i*nf[var], color=colors[net_key], lw=1)
            ax.plot(t_traj, var_traj*nf[var], ls='-', color=colors[net_key], lw=2, alpha=0.8)
            xmax = np.max((xmax, np.max(data[net_id][var].keys())-t_iv))
            ymax = np.max((ymax, np.max(data[net_id][var].values())/sf_i*nf[var]))
        ax.set_xlim([-4, xmax+2])
        ax.set_xticks(np.arange(0,xmax+2,5, dtype=int))
        if exp_id=='INCELL_WT':
            ylim = [-0.1*ymax, ymax*1.5]
            ax.set_ylim(ylim)
            if var=='pEtot':
                ylim_for_ko = ylim
        else:
            ax.set_ylim(ylim_for_ko)
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_position(('data', -1))
        ax.tick_params(axis='x', labelsize= 14) 
        ax.tick_params(axis='y', labelsize= 14) 
        xlim = ax.get_xlim
        ax.set_xlim([-4,17])
        
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
#         ax.legend(frameon=False, bbox_to_anchor=(1,0.8), loc='upper left', fontsize=fontsize)
        plt.tight_layout()
        plt.savefig('../../res/Figure6C_incell_'+exp_id+'_'+var+'.pdf')
        plt.savefig('../../res/Figure6C_incell_'+exp_id+'_'+var+'.eps')
    plt.close("all")
    for var in vars:
        fig, ax = plt.subplots(1, 1, figsize=(5,3), sharey=True)
        sf_i = sf[exp_id][var]
        ymax = 0
        xmax = 0
        for net_key, net in nets[exp_id].items():
            net.set_var_ics(params)
            net.set_var_vals(params) 
            net_id = net.id
            t_data = np.array(exps[exp_id].GetData()[net_id][var].keys())
            var_data = np.array([d[0] for d in exps[exp_id].GetData()[net_id][var].values()])
            traj = net.integrate(times)
            t_traj = traj.timepoints[traj.timepoints>=999] - t_iv
            traj_set = Ensembles.ensemble_trajs(net, times, ens[::10])
            lower, upper = Ensembles.traj_ensemble_quantiles(traj_set, (0.0, 1.0))
            var_traj = traj.get_var_traj(var)[traj.timepoints>=999]
            ax.plot(t_traj, var_traj*nf[var], ls='-', color=colors[net_key], lw=2, alpha=0.8)
            lower_var = lower.get_var_traj(var)[lower.timepoints>=999]*nf[var]
            upper_var = upper.get_var_traj(var)[upper.timepoints>=999]*nf[var]
            time_traj_fill = lower.timepoints[lower.timepoints>=999] - t_iv
            ax.fill_between(time_traj_fill, lower_var, upper_var, alpha = 0.2)
            ax.scatter(np.sort(t_data)-t_iv, var_data[np.argsort(t_data)]/sf_i*nf[var], color=colors[net_key], lw=1)
            xmax = np.max((xmax, np.max(data[net_id][var].keys())-t_iv))
            ymax = np.max((ymax, np.max(data[net_id][var].values())/sf_i*nf[var]))
        ax.set_xlim([-4, xmax+2])
        ax.set_xticks(np.arange(0,xmax+2,5, dtype=int))
        if exp_id=='INCELL_WT':
            ylim = [-0.1*ymax, ymax*1.5]
            ax.set_ylim(ylim)
            if var=='pEtot':
                ylim_for_ko = ylim
        else:
            ax.set_ylim(ylim_for_ko)
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_position(('data', -1))
        ax.tick_params(axis='x', labelsize= 14) 
        ax.tick_params(axis='y', labelsize= 14) 
        xlim = ax.get_xlim
        ax.set_xlim([-4,17])
        
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
#         ax.legend(frameon=False, bbox_to_anchor=(1,0.8), loc='upper left', fontsize=fontsize)
        plt.tight_layout()
        plt.savefig('../../res/Figure6C_incell_ens_'+exp_id+'_'+var+'.pdf')
        plt.savefig('../../res/Figure6C_incell_ens_'+exp_id+'_'+var+'.eps')
