import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import pickle
import copy
import os
from expts import *
import SCfunctions as sc

params_file = "../../data/invitro/params_invitro.pickle"
global_it = 10000

if not os.path.exists(params_file):
    exp_ids = ['MKK1_ERK_RSK', 'PERK_PTASE']
    pars_fixed = KeyedList()

    with open('../../data/SPR/params_SPR_all.pickle', 'rb') as f:
        pars_SPR, sf_SPR = pickle.load(f)

    pars_constrained = KeyedList([
        ('KD_ER', 2.5),
        ('KD_OR', 1.2e-3),
        ('KD_OpR', 1.2e-3),
        ('KD_pEO', 0.8),
        ('KD_EK', 4.0),
        ('pE_factor', 100),
    ]) + pars_SPR

    pars_free = KeyedList([
        ('kp_R', 1.0),
        ('kp_E', 1.0),
        ('kon_EK', 1.0),
        ('KD_EP', 1.0),
        ('kon_EP', 1.0),
        ('kdp_E', 1.0)
    ])

    model, popt, ens, cost = sc.fit_exps(exps,
                                 nets,
                                 pars_fixed,
                                 pars_constrained,
                                 pars_free,
                                 exp_ids=exp_ids,
                                 global_fit=True,
                                 global_it=global_it,
                                 return_ens=True)
    
    with open(params_file, 'w') as f:
        pickle.dump((model, popt, ens, cost), f)
else: 
    with open(params_file, 'r') as f:
        model, popt, ens, cost = pickle.load(f)
    
# plot
exp_ids = ['MKK1_ERK_RSK', 'PERK_PTASE']
colors = {'K': 'b', '+ORF45': 'r', '+RSK': 'gray', '+RSK+ORF(1)': 'y', '+RSK+ORF(5)': 'r'}
fontsize = 14

for exp_id in exp_ids:
    params = popt
    model.params.update(params)
    model.CalculateForAllDataPoints(params)
    sf = model.compute_scale_factors(1)[0]
    
    data = exps[exp_id].GetData()
    vars = list(np.unique([var for key in data.keys() for var in data[key].keys()]))
    
    for var in vars:
        fig,ax = plt.subplots(1, 1, figsize=(5, 3))
        sf_i = sf[exp_id][var]
        for net_key, net in nets[exp_id].items():
            if net_key == '+RSK+ORF(1)':
                continue
            net.set_var_ics(params)
            net.set_var_vals(params) 
            net_id = net.id
            t_data = np.array(exps[exp_id].GetData()[net_id][var].keys())
            var_data = np.array([d[0] for d in exps[exp_id].GetData()[net_id][var].values()])
            traj = net.integrate([0, 999, np.max(t_data)*1.1])        
            t_traj = traj.timepoints[traj.timepoints>=999] - 1000
            var_traj = traj.get_var_traj(var)[traj.timepoints>=999]
            ax.scatter(np.sort(t_data)-1000, var_data[np.argsort(t_data)]/sf_i, color=colors[net_key])
            ax.plot(t_traj, var_traj, ls='-', color=colors[net_key], lw=2, alpha=0.8)
        ax.set_xlim([-4, 22])
        ax.set_xticks(np.arange(0,24,5))
        ymax = np.round(ax.get_ylim()[1])
        ystep = 0.25 if ymax<=1 else 0.5
        ax.set_yticks(np.arange(0, ymax+ystep, ystep))
        
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_position(('data', -1))

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticklabels(['{:d}'.format(x) for x in ax.get_xticks()], fontsize=14)
        fs = '{:.'+str(len(str(ystep)) - 2)+'f}'
        ax.set_yticklabels([fs.format(y) for y in ax.get_yticks()], fontsize=14)        
        plt.tight_layout()
        plt.savefig('../../res/Figure6B_invitro_'+exp_id+'_'+var+'.pdf')
        plt.savefig('../../res/Figure6B_invitro_'+exp_id+'_'+var+'.eps')
