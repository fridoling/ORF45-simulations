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
popt = OrderedDict()
with open(params_file, 'r') as f:
    model, popt['WT'], ens, cost = pickle.load(f)

net_K = nets['INCELL_WT']['K']
net_ORF = nets['INCELL_WT']['ORF']


times = [0, 999, 1025*1.1]
t_iv = 1000

y_labels = {'pRtot': u'pRSK ($\mu M$)', 'pEtot': u'ppERK ($\mu M$)'}

popt['pEO_EO'] = popt['WT'].copy()
popt['pEO_EO'].setByKey('KD_pEO', popt['WT'].getByKey('KD_pEO')*popt['WT'].getByKey('pE_factor'))
# popt['epo_eo'].setByKey('kon_EO', 0)
# popt['epo_eo'].setByKey('KD_EO', 0)
# popt['epo_eo'].setByKey('KD_pEO', 0)

popt['KD_OR_KD_OpR'] = popt['WT'].copy()
popt['KD_OR_KD_OpR'].setByKey('KD_OR', 10)
popt['KD_OR_KD_OpR'].setByKey('KD_OpR', 10)
# popt['KDor_KDorp'].setByKey('kon_OR', 0)
# popt['KDor_KDorp'].setByKey('KD_OR', 0)
# popt['KDor_KDorp'].setByKey('KD_OpR', 0)

# popt['ad1'] = popt['WT'].copy()
# popt['ad1'].setByKey('a', 1)
# popt['ad1'].setByKey('d', 1)

fontsize = 14

main_vars = ['ERK', 'RSK']
plot_vars = {'ERK': ['pEtot', 'pE_FSUB'], 'RSK': ['pRtot', 'pR_VFSUB']}
titles = {'ERK': {'pEtot': "total pERK", 'pE_FSUB': 'FSUB competent pERK'},
          'RSK': {'pRtot': "total pRSK/SUB competent pRSK", 'pR_VFSUB': 'VFSUB competent pRSK'}
         }
nf = {'RSK': 2.0, 'ERK': 0.7}
xlen = len(main_vars)
ylen = np.max([len(vars) for vars in plot_vars])

cols = {'WT': 'red', 'KD_OR_KD_OpR': 'green', 'ad1': 'purple', 'pEO_EO': 'magenta'}
# labels = {'WT': 'WT', 'KDor_KDorp': 'no RSK:ORF45', 'ad1': 'no closed complex', 'epo_eo': 'no ERK:ORF45'}
labels = {'WT': 'WT', 'KD_OR_KD_OpR': 'KD_OR = KD_OpR = $10\mu M$', 'ad1': 'a = d = 1', 'pEO_EO': 'KD_EO = KD_pEO'}
 
Otot_vec = np.logspace(-2, 2, 51) 
trajs = {}
for key, params in popt.items():
    trajs[key] = np.empty_like(Otot_vec, dtype = object)
    for k, Otot in zip(range(len(Otot_vec)), Otot_vec):
        net = net_ORF.copy()
        net.set_var_ic('Otot', Otot)
        trajs[key][k] = net.integrate([0,10000], params)

for j, main_var in zip(range(len(main_vars)), main_vars):
    for i, var in zip(range(len(plot_vars[main_var])), plot_vars[main_var]):
        fig, ax = plt.subplots(1,1,figsize = (4,3))
        if var is None:
            ax.axis('off')
            continue
        for key in popt.keys():
            sub_vec = np.array([traj.get_var_traj(var)[-1] for traj in trajs[key]])                
            ax.plot(Otot_vec, sub_vec/nf[main_var], lw=2, alpha=0.8, color = cols[key], label = labels[key])
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
#       #  This doesn't work with log scale:
#         ax.spines['left'].set_position(('data', 1e-3))        
        # Only show ticks on the left and bottom spines
        ax.set_xlabel("ORF45 concentration ($\mu M$)")
        ax.set_xscale('log')
        ax.set_ylim([0,1] if 'tot' in var else [0,0.1])
#         ax.axvline(x = 2, color = 'k', lw = 1, ls = '--')
        ax.set_ylabel("normalized concentration")
        plt.tight_layout()    
        plt.savefig('../../res/Figure6D_mutants_'+var+'.pdf')
        plt.savefig('../../res/Figure6D_mutants_'+var+'.eps')
