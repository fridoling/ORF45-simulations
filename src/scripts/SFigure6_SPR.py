## Fit and plot SPR experiments

# import required modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *
import re
from collections import OrderedDict
import pickle
import copy
import pandas as pd
import os
import functions

# should we plot labels on individual plots?
plot_labels = False
       
net_basic = IO.from_SBML_file('../models/model_SPR.xml','net_basic')
net_basic_ppERK_ORF = IO.from_SBML_file('../models/model_SPR.xml','net_ppERK_ORF')

data = {}
models = {}

with open('../../data/SPR/SPR_data.pickle', 'r') as f:
    df, t_iv = pickle.load(f)

if os.path.exists("../../data/SPR/params_SPR.pickle"):
    with open("../../data/SPR/params_SPR.pickle", 'r') as f:
        params_opt, sf_opt, data, nets, exps = pickle.load(f)
else:
    nets = {}
    exps = {}
    params_init = {}
    params_opt = {}
    sf_opt = {}
    
    big_error = 1
    small_error = 0.1
    
    exp_ids_dict = {}
    
    ### RSK ppERK
    m_id = 'RSK_ppERK'
    
    dyn_vars = [
         'pER',
    ]
    
    nonvars = list(set(net_basic.dynamicVars.keys()) - set(dyn_vars))
        
    params_init[m_id] = KeyedList([
             ('kon_ER', 228),
        ])
    
    for net_id, pE in zip(['RSK_ppERK_012', 'RSK_ppERK_111'], [0.12, 1.11]):
        net = net_basic.copy(new_id=net_id)    
        non_opt = list(set(net_basic.optimizableVars.keys()) - set(params_init[m_id].keys()))    
        for key in non_opt:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)
        for var, val in params_init[m_id].items():
            net.set_var_optimizable(var, True) 
        for key in nonvars:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)                
        net.set_var_ic('pE', pE)
        net.set_var_ic('O', 0)    
        net.set_var_ic('Rtot', 1)
        net.set_var_ic('kon_EO', 0)
        net.set_var_ic('KD_OR', 0)
        net.set_var_ic('KD_ER', 2.5) 
        net.set_var_ic('KD_EO', 0)
        net.set_var_ic('KD_pEO', 0)       
        net.set_var_constant('kon_ER_factor', False)
        net.add_event(id='kon_ER_off', trigger='gt(time, '+str(t_iv[net_id])+')', event_assignments={'kon_ER_factor': 0})
        nets[net_id] = net
        
    exp_ids_dict['RSK_ppERK'] = []
    for exp_id in ['RSK_ppERK_012', 'RSK_ppERK_111']:
        df_exp = df[exp_id]
        data[exp_id] = {
            exp_id: {
                'pER': {t: (v, 0) for t,v in zip(df_exp['time'], df_exp[exp_id])}
                }
            }
        iv_exp = [0, t_iv[exp_id]+50]
        data[exp_id] = functions.set_errors(data[exp_id], error=big_error, iv=iv_exp, iv_error=small_error, absolute=True, min_error=0.1)
        exp = Experiment(name=exp_id, data=data[exp_id])
        exps[exp_id] = exp
        exp_ids_dict['RSK_ppERK'].append(exp_id)
    
    ### RSK ORF
    m_id = 'RSK_ORF'
    
    dyn_vars = [
        'OR',
    ]
    
    nonvars = list(set(net_basic.dynamicVars.keys()) - set(dyn_vars))
    
    params_init[m_id] = KeyedList([
         ('kon_OR', 1),
        ])
    
    for net_id, O in zip(['RSK_ORF_05', 'RSK_ORF_135'], [5e-4, 1.35e-2]):
        net = net_basic.copy(new_id=net_id)
        non_opt = list(set(net_basic.optimizableVars.keys()) - set(params_init[m_id].keys()))        
        for key in non_opt:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)
        for var, val in params_init[m_id].items():
            net.set_var_optimizable(var, True) 
        for key in nonvars:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)                
        net.set_var_ic('pE', 0.0)
        net.set_var_ic('Rtot', 1)
        net.set_var_ic('O', O)
        net.set_var_ic('KD_OR', 1.2e-3)    
        net.set_var_constant('kon_OR_factor', False)
        net.add_event(id='kon_OR_off', trigger='gt(time, '+str(t_iv[net_id])+')', event_assignments={'kon_OR_factor': 0})
        nets[net_id] = net
    
    exp_ids_dict[m_id] = []
    for exp_id in ['RSK_ORF_05', 'RSK_ORF_135']:
        df_exp = df[exp_id]    
        data[exp_id] = {
            exp_id: {
                'OR': {t: (v, 0) for t,v in zip(df_exp['time'], df_exp[exp_id])}
                }
            }
        iv_exp = [0, t_iv[exp_id]+50]
        data[exp_id] = functions.set_errors(data[exp_id], error=big_error, iv=iv_exp, iv_error=small_error, absolute=True, min_error=0.1)
        exp = Experiment(name=exp_id, data=data[exp_id])
        exps[exp_id] = exp
        exp_ids_dict[m_id].append(exp_id)
    
    ### ppERK ORF
    m_id = 'ppERK_ORF'
    
    dyn_vars = [
        'pEO',
    ]
    
    nonvars = list(set(net_basic_ppERK_ORF.dynamicVars.keys()) - set(dyn_vars))
    
    params_init[m_id] = KeyedList([
         ('kon_EO', 1),
        ])
    
    pE = 1.66
    net = net_basic_ppERK_ORF.copy(new_id=m_id)
    non_opt = list(set(net_basic_ppERK_ORF.optimizableVars.keys()) - set(params_init[m_id].keys()))        
    for key in non_opt:
        net.set_var_constant(key, is_constant=True)
        net.set_var_optimizable(key, is_optimizable=False)
    for var, val in params_init[m_id].items():
        net.set_var_optimizable(var, True) 
    for key in nonvars:
        net.set_var_constant(key, is_constant=True)
        net.set_var_optimizable(key, is_optimizable=False)
    net.set_var_ic('pE', pE)
    net.set_var_ic('O', 1.0)
    net.set_var_ic('Rtot', 0.0)
    net.set_var_ic('KD_pEO', 0.8)    
    net.set_var_constant('kon_EO_factor', False)
    net.add_event(id='kon_EO_off', trigger='gt(time, '+str(t_iv[m_id])+')', event_assignments={'kon_EO_factor': 0})
    nets[m_id] = net
    
    exp_ids_dict[m_id] = []
    exp_id = m_id
    df_exp = df[exp_id] 
    data[exp_id] = {
        exp_id: {
            'pEO': {t: (v, 0) for t,v in zip(df_exp['time'], df_exp[exp_id])}
            }
        }
    iv_exp = [0, t_iv[exp_id]+50]
    data[exp_id] = functions.set_errors(data[exp_id], error=big_error, iv=iv_exp, iv_error=small_error, absolute=True, min_error=0.1)
    exp = Experiment(name=exp_id, data=data[exp_id])
    exps[exp_id] = exp
    exp_ids_dict[m_id].append(exp_id)
    
    binaries = [
        'RSK_ppERK',    
        'RSK_ORF',    
        'ppERK_ORF'
    ]
    
    for m_id in binaries:
        print "fitting "+m_id+"..."
        exp_ids = exp_ids_dict[m_id]
    
        exp_set = [exps[exp_id] for exp_id in exp_ids]
        net_set = [nets[exp_id] for exp_id in exp_ids]
    
        models[m_id] = Model(exp_set, net_set)
    
        for var, val in params_init[m_id].items():
            res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(10)))
            models[m_id].AddResidual(res)
        m = models[m_id]
        params_opt[m_id] = Optimization.fmin_lm_log_params(m, params=params_init[m_id], maxiter=100, disp=True)
        sf_opt[m_id] = m.GetScaleFactors()        
    
    
    ### RSK ppERK ORF with parameters from binary experiments
    print "Fitting RSK ppERK ORF...\n"
    
    m_id = 'RSK_ppERK_ORF'
    
    kon_OR = params_opt['RSK_ORF'].getByKey('kon_OR')
    kon_ER = params_opt['RSK_ppERK'].getByKey('kon_ER')
    kon_EO = params_opt['ppERK_ORF'].getByKey('kon_EO')
    
    params_binary = KeyedList([
            ('kon_ER', kon_ER),
            ('kon_OR', kon_OR),
            ('kon_EO', kon_EO),
            ])
    
    params_init[m_id] = KeyedList([
            ('a', 0.007),
             ('d', 0.0003)
            ])
    
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
    
    for net_id, pE in zip(['RSK_ppERK_ORF_012', 'RSK_ppERK_ORF_111'], [0.12, 1.11]):
        net = net_basic.copy(new_id=net_id)
        non_opt = list(set(net.optimizableVars.keys()) - set(params_init[m_id].keys()))    
        for key in non_opt:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)
        for var, val in params_init[m_id].items():
            net.set_var_optimizable(var, True) 
        for key in nonvars:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)
        net.set_var_ics(params_binary)
        net.set_var_ics(params_init[m_id])
        net.set_var_ic('E', 0.0)
        net.set_var_ic('pE', pE)    
        net.set_var_ic('O', 10)
        net.set_var_ic('Rtot', 1)
        net.set_var_ic('KD_OR', 1.2e-3)
        net.set_var_ic('KD_ER', 2.5)   
        net.set_var_ic('KD_pEO', 0.8) 
        net.set_var_ic('KD_EO', 100)        
        net.set_var_constant('kon_ER_factor', False)
        net.set_var_constant('kon_OR_factor', False)    
        net.set_var_constant('kon_EO_factor', False)        
        net.add_event(id='flow_off', trigger='gt(time, '+str(t_iv[net_id])+')', event_assignments={'kon_ER_factor': 0, 'kon_OR_factor': 0, 'kon_EO_factor': 0, 'pEO': 0})
        nets[net_id] = net
        
    exp_ids_dict[m_id] = []
    for exp_id in ['RSK_ppERK_ORF_012', 'RSK_ppERK_ORF_111']:
        df_exp = df[exp_id]
        data[exp_id] = {
            exp_id: {
                'Rcomp': {t: (v, 0) for t,v in zip(df_exp['time'], df_exp[exp_id])}
                }
            }
        iv_exp = [t_iv[exp_id]-2, t_iv[exp_id]+50]
        data[exp_id] = functions.set_errors(data[exp_id], error=big_error, iv=iv_exp, iv_error=small_error, absolute=True, min_error=0.1)
        exp = Experiment(name=exp_id, data=data[exp_id])
        exps[exp_id] = exp
        exp_ids_dict[m_id].append(exp_id)
    
    exp_ids = exp_ids_dict[m_id]

    exp_set = [exps[exp_id] for exp_id in exp_ids]
    net_set = [nets[exp_id] for exp_id in exp_ids]

    models[m_id] = Model(exp_set, net_set)

    for var, val in params_init[m_id].items():      
        if var in ['a', 'd']:
            res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(1000)))
        else:
            res = Residuals.PriorInLog(var+'_prior', var, np.log(val), np.log(np.sqrt(1.1)))
        models[m_id].AddResidual(res)
    m = models[m_id]
    params_opt[m_id] = Optimization.fmin_lm_log_params(m, params=params_init[m_id], maxiter=10, disp=True)
    sf_opt[m_id] = m.GetScaleFactors()        
    with open("../../data/SPR/params_SPR.pickle", 'w') as f:
         pickle.dump((params_opt, sf_opt, data, nets, exps), f)

params_SPR_all = params_opt['RSK_ORF'] + \
    params_opt['RSK_ppERK'] + \
    params_opt['ppERK_ORF'] + \
    KeyedList([(p, params_opt['RSK_ppERK_ORF'].getByKey(p)) for p in ['a', 'd']])
with open("../../data/SPR/params_SPR_all.pickle", 'w') as f:
    pickle.dump((params_SPR_all, sf_opt['RSK_ppERK_ORF']), f)

## Save parameters for table
table_pars_SPR = ['koff_ER', 'koff_OR', 'koff_pEO', 'a', 'd']
net = nets['RSK_ppERK_ORF_111']
net.set_var_vals(params_SPR_all)
table_dict_SPR = { par: net.get_var_val(par) for par in table_pars_SPR }
with open("../../res/table_pars_SPR.pickle", "wb") as f:
    pickle.dump(table_dict_SPR, f)


## Plot
print "Plotting...\n"
off_rates = {
        'RSK_ppERK': ['koff_ER'],
        'RSK_ORF': ['koff_OR'],
        'ppERK_ORF': ['koff_pEO'],
        'RSK_ppERK_ORF': ['koff_ER', 'koff_OR', 'koff_pEO', 'a', 'd'],
      }


titles = [
    u'RSK with ppERK = $0.12\mu M$',
    u'RSK with ppERK = $1.11\mu M$',
    u'RSK with ORF45 = $0.5\mu M$',
    u'RSK with ORF45 = $13.5\mu M$',
    u'ORF45 with ppERK = $1.66\mu M$',
    '',
    u'RSK with ORF45 = $10\mu M$ and ppERK = $0.12\mu M$',
    u'RSK with ORF45 = $10\mu M$ and ppERK = $1.11\mu M$',
]


m_ids = [
        'RSK_ppERK',
        'RSK_ORF',
        'ppERK_ORF',
        'RSK_ppERK_ORF',
        ]
exp_ids_dict = {}
for m_id in m_ids:
    exp_ids_dict[m_id] = []
    for exp_id in exps.keys():
        if re.match('^'+m_id, exp_id) and re.match(m_id+'_[A-Z]', exp_id) is None:
            exp_ids_dict[m_id].append(exp_id)


for m_id in m_ids:
    exp_ids = exp_ids_dict[m_id]
    for exp_id in exp_ids:
        fig,ax = plt.subplots(1, 1, figsize=(5, 3))
        params = params_opt[m_id]
        var = data[exp_id][data[exp_id].keys()[0]].keys()[0]
        for net_id in data[exp_id].keys():
            sf_i = sf_opt[m_id][net_id][var]
            nets[net_id].set_var_ics(params)
            nets[net_id].set_var_vals(params) 
            t_data = np.array(exps[net_id].GetData()[net_id][var].keys())
            var_data = np.array([d[0] for d in exps[net_id].GetData()[net_id][var].values()])
            traj = nets[net_id].integrate([0, np.max(t_data)*1.1])        
        ax.scatter(t_data, var_data, color='k', s=10, marker='o', label='data')
        ax.plot(traj.timepoints, traj.get_var_traj(var)*sf_i, ls=':', label=m_id, color='r')    
        # Hide the right and top spines
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_position(('data', -1))
        
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.set_xlim([-4,25])
        s = ''
        if plot_labels:
            for net_id in data[exp_id].keys():
                for var in off_rates[m_id]:
                    val = nets[net_id].get_var_val(var)
                    s=s+var+' = '+'{:.2g}'.format(val)+'\n'
            ax.text(0.8, 0.9, s, transform=ax.transAxes, verticalalignment='top')
        plt.tight_layout()
        plt.savefig('../../res/SFigure6_SPR_'+exp_id+'.pdf')
        plt.savefig('../../res/SFigure6_SPR_'+exp_id+'.eps')
