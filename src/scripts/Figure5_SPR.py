import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import pickle
import copy
import os
import functions


# load model
net_basic = IO.from_SBML_file('./models/model_SPR.xml','net_basic')

# load data
if not os.path.exists('../data/Figure5_SPR_data.pickle'):
    df_all = pd.read_excel('../data/SPR/SPR_data1.xlsx', sheet_name=3)
    df_all.columns = [
        'time',
        'RSK_ppERK_ORF_012',
        'RSK_ppERK_ORF_012_fit',
        'RSK_ppERK_ORF_111',
        'RSK_ppERK_ORF_111_fit']
    df_all = df_all[df_all['time']>0]
    t_iv = {}
    ind_max_012 = np.argmax(df_all['RSK_ppERK_ORF_012_fit'].diff()<0)
    ind_max_111 = np.argmax(df_all['RSK_ppERK_ORF_111_fit'].diff()<0)
    t_iv['RSK_ppERK_ORF_012'] = df_all['time'][ind_max_012]
    t_iv['RSK_ppERK_ORF_111'] = df_all['time'][ind_max_111]
    row_select = np.int(df_all.shape[0]/200)
    df_all = df_all.loc[::row_select]
    df_all[df_all<0] = 0
    with open('../data/Figure5_SPR_data.pickle', 'w') as f:
        pickle.dump((df_all, t_iv), f)
else:
    with open('../data/Figure5_SPR_data.pickle', 'r') as f:
        df_all, t_iv = pickle.load(f)

df = {}
df['RSK_ppERK_ORF_012'] = df_all[['time']]
df['RSK_ppERK_ORF_012'].insert(1, 'RSK_ppERK_ORF_012',
df_all['RSK_ppERK_ORF_012'].divide(np.max(df_all['RSK_ppERK_ORF_012_fit'])))

df['RSK_ppERK_ORF_111'] = df_all[['time']]
df['RSK_ppERK_ORF_111'].insert(1, 'RSK_ppERK_ORF_111', df_all['RSK_ppERK_ORF_111'].divide(np.max(df_all['RSK_ppERK_ORF_111_fit'])))


big_error = 1
small_error = 0.1

if not os.path.exists('../data/SPR/Figure5_SPR_params.pickle'):    
    KD_ER = 2.5
    KD_OR = 1.2e-3
    KD_pEO = 0.8

    with open('../data/SPR/params_SPR.pickle') as f:
        params_i = pickle.load(f)
        
    dyn_vars = [
        'ER',
        'ER_OR',
        'EO',
        'EO_OR',
        'EO_ER',
        'EO_ER_OR',    
        'pER',
        'OR',
        'pER_OR',
        'pEO',
        'pEO_OR',
        'pEO_pER', 
        'pEO_pER_OR'
    ]

    nonvars = list(set(net_basic.dynamicVars.keys()) - set(dyn_vars))
    nets = {}
    exps = {}
    models = {}
    data = {}
    params_opt = {}
    sf_opt = {}
    
    with open('../data/SPR/params_SPR_all.pickle', 'rb') as f:
        params_opt['RSK_ppERK_ORF_111'], sf_opt['RSK_ppERK_ORF_111'] = pickle.load(f)


    # Now with fixed a = d = 1
    params_i = KeyedList([
            ('kon_ER', params_opt['RSK_ppERK_ORF_111'].getByKey('kon_ER')),
            ('kon_OR', params_opt['RSK_ppERK_ORF_111'].getByKey('kon_OR')),
            ('kon_EO', params_opt['RSK_ppERK_ORF_111'].getByKey('kon_EO'))
            ])

    net_id = 'RSK_ppERK_ORF_111'
    ep = 1.11
    net = net_basic.copy(new_id=net_id)
    non_opt = list(set(net.optimizableVars.keys()) - set(params_i.keys()))    
    for key in non_opt:
        net.set_var_constant(key, is_constant=True)
        net.set_var_optimizable(key, is_optimizable=False)
    for var, val in params_i.items():
        net.set_var_optimizable(var, True) 
    for key in nonvars:
        net.set_var_constant(key, is_constant=True)
        net.set_var_optimizable(key, is_optimizable=False)
    net.set_var_ics(params_i)
    net.set_var_ic('a', 1)
    net.set_var_ic('d', 1)
    net.set_var_ic('pE', ep) 
    net.set_var_ic('O', 10)
    net.set_var_ic('Rtot', 1)
    net.set_var_ic('KD_OR', KD_OR)
    net.set_var_ic('KD_ER', KD_ER)   
    net.set_var_ic('KD_pEO', KD_pEO)       
    net.set_var_constant('kon_ER_factor', False)
    net.set_var_constant('kon_OR_factor', False)    
    net.set_var_constant('kon_EO_factor', False)        
    net.add_event(id='flow_off', trigger='gt(time, '+str(t_iv[net_id])+')',
                  event_assignments={'kon_ER_factor': 0, 'kon_OR_factor': 0, 'kon_EO_factor': 0, 'pEO': 0})
    nets[net_id] = net

    exp_id = 'RSK_ppERK_ORF_111'
    data[exp_id] = {
        exp_id: {
            'Rcomp': {t: (v, 0) for t,v in zip(df[exp_id]['time'], df[exp_id][exp_id])}
            }
        }
    iv_exp = [t_iv[exp_id]-2, t_iv[exp_id]+50]
    data[exp_id] = functions.set_errors(data[exp_id], error=big_error, iv=iv_exp, iv_error=small_error, absolute=True, min_error=0.1)
    exp = Experiment(name=exp_id, data=data[exp_id])
    exps[exp_id] = exp

    m_id = 'RSK_ppERK_ORF_111'
    exp_ids = [m_id]

    exp_set = [exps[exp_id] for exp_id in exp_ids]
    net_set = [nets[exp_id] for exp_id in exp_ids]

    models[m_id] = Model(exp_set, net_set)

    for var, val in params_i.items():        
        if var in ['a','d']:
            res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(1000)))
        else:
            res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(1.1)))
        models[m_id].AddResidual(res)
    m = models[m_id]
    params_opt[m_id+"_noC"] = Optimization.fmin_lm_log_params(m, params=params_i, maxiter=10, disp=False)
    sf_opt[m_id+"_noC"] = m.GetScaleFactors()
        
    with open('../data/Figure5_SPR_params.pickle', 'w') as f:
        pickle.dump((params_opt, sf_opt, nets, exps), f)
else:
    with open('../data/Figure5_SPR_params.pickle', 'r') as f:
        params_opt, sf_opt, nets, exps = pickle.load(f)
    data = {}
    for exp_id in ['RSK_ppERK_ORF_111']:
        data[exp_id] = {
            exp_id: {
                'Rcomp': {t: (v, 0) for t,v in zip(df[exp_id]['time'], df[exp_id][exp_id])}
                }
            }
        iv_exp = [t_iv[exp_id]-2, t_iv[exp_id]+50]
        data[exp_id] = functions.set_errors(data[exp_id], error=big_error, iv=iv_exp, iv_error=small_error, absolute=True, min_error=0.1)
        
        
# plot
params_ids = [
    'RSK_ppERK_ORF_111',
    'RSK_ppERK_ORF_111_noC']
m_ids = ['RSK_ppERK_ORF_111', 'RSK_ppERK_ORF_111']

colors = ['b', 'r']

fig,axes = plt.subplots(1, 1, figsize=(5,3), sharey=True)

for i, params_id, col in zip(range(len(m_ids)), params_ids, colors):
    m_id = m_ids[i]
    params = params_opt[params_id]
    ax = axes
    var = data[m_id][data[m_id].keys()[0]].keys()[0]
    for net_id in data[m_id].keys():
        sf_i = sf_opt[params_id][net_id][var]
        nets[m_id].set_var_ics(params)
        nets[m_id].set_var_vals(params) 
        if '_noC' in params_id:
            nets[m_id].set_var_ics(KeyedList([('a', 1), ('d', 1)]))
            nets[m_id].set_var_vals(KeyedList([('a', 1), ('d', 1)]))            
        t_data = np.array(exps[net_id].GetData()[net_id][var].keys())
        var_data = np.array([d[0] for d in exps[net_id].GetData()[net_id][var].values()])
        traj = nets[m_id].integrate([0, np.max(t_data)*1.1])        
    ax.plot(np.sort(t_data), var_data[np.argsort(t_data)], color='k', label='data', lw=1)
    ax.plot(traj.timepoints, traj.get_var_traj(var)*sf_i, ls='--', label=m_id, color=col, lw=2, alpha=0.8)    
#     ax.set_xlabel('time (minutes)')
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.set_xlim([-3,17])
    ax.set_ylim([-0.2,1.5])
    ax.spines['left'].set_position(('data', -1))
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
# ax.set_ylabel('normalized RU')
xticks = np.arange(0,20,5)
yticks = np.round(np.arange(0,1.6,0.2), 2)
ax.set_xticks(xticks)
ax.set_xticklabels(xticks, fontsize=14)
ax.set_yticks(yticks)
ax.set_yticklabels(yticks, fontsize=14)
plt.tight_layout()
plt.savefig('../res/Figure5_SPR.pdf')
plt.savefig('../res/Figure5_SPR.eps')
