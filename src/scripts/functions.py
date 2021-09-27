## Functions used by other scripts

# import required modules
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt
from SloppyCell.ReactionNetworks import *

def set_errors(data, error=0.1, iv=[0,0], iv_error=0, min_error=0.1, absolute=False):
    data_out = copy.deepcopy(data)
    for net,vars in data.items():
        for var,times in vars.items():
            for t,tup in times.items():
                if absolute:
                    if t>iv[0] and t<=iv[1]:
                        data_out[net][var][t] = (tup[0], iv_error)
                    else:
                        data_out[net][var][t] = (tup[0], error)
                else:
                    if t>iv[0] and t<=iv[1]:                    
                        data_out[net][var][t] = (tup[0], np.max([iv_error*tup[0],min_error]))
                    else:
                        data_out[net][var][t] = (tup[0], np.max([error*tup[0],min_error]))                        
    return data_out

def subtract_bg(data, bg=0, norm=True, renormalize=True):
    data_out = copy.deepcopy(data)
    for net, vars in data.items():
        for var,times in vars.items():
            for t,tup in times.items():
                if renormalize:
                    data_out[net][var][t] = ((tup[0]-bg)/(1-bg), tup[1])
                else:
                    data_out[net][var][t] = (tup[0]-bg, tup[1])
    return data_out   

def order_params(params, m):
   params_ordered = KeyedList()
   for p in m.params.keys():
       val = params.getByKey(p)
       params_ordered.setByKey(p, val)
   return params_ordered
 
def plot_ens(ens, net, out_vars, params_opt = None, step = 10, file = None, axis = None, mode = "std", normalize = True):
    ens = ens[::step]
    ens_out = np.empty((len(ens), len(out_vars)))
    for i in range(len(ens)):
        ens_i = ens[i]
        net.set_var_vals(ens_i)
        ens_out_i = [net.get_var_val(var) for var in out_vars]
        ens_out[i,] = ens_out_i
    ens_mean = np.median(ens_out, axis = 0)
    if normalize:
        ens_out = ens_out/ens_mean
        ens_mean = np.log10(ens_mean)
        ens_mean_plot = np.zeros_like(ens_mean)
    else:
        ens_mean_plot = np.log10(ens_mean)
    if mode=="std":
        ens_error = np.std(np.log10(ens_out) - ens_mean, axis = 0)
    elif mode=="quantile":
        ens_error = np.abs(np.quantile(np.log10(ens_out) - ens_mean, [0.05, 0.95], axis = 0))
    else:
        raise ValueError("mode must either be 'std' or 'quantile'")
    if params_opt is not None:
        net.set_var_vals(params_opt)
        vars_opt = np.array([net.get_var_val(var) for var in out_vars])
        if normalize:
            vars_opt = np.log10(vars_opt) - ens_mean
        else:
            vars_opt = np.log10(vars_opt)
    if axis is None:
        fig, axis = plt.subplots(figsize = (len(out_vars)*0.5,3))
    if not normalize:
        axis.scatter(range(len(out_vars)), ens_mean_plot, color = "k", label = "median")
    if params_opt is not None:
        axis.scatter(range(len(out_vars)), vars_opt, color = "r", marker = "x", label = "best")  
    axis.errorbar(range(len(out_vars)), ens_mean_plot, ens_error, ls = "none", capsize = 3, color = "k")
    axis.set_xticks(range(len(out_vars)));
    axis.set_xticklabels(out_vars, rotation = 45);
    if mode=="std":
        ymin = np.floor(np.min(ens_mean_plot-ens_error))
        ymax = np.ceil(np.max(ens_mean_plot+ens_error))
    else:
        ymin = np.floor(np.min(ens_mean_plot-ens_error[0,:]))
        ymax = np.ceil(np.max(ens_mean_plot+ens_error[1,:]))
    axis.set_yticks(np.arange(ymin,ymax+1));
    axis.set_yticklabels([r"$10^{"+str(int(n))+"}$" for n in np.arange(ymin, ymax+1)]);
    plt.tight_layout()
    if file is not None:
        plt.savefig(file)
        
def plot_cor_matrix(ens, net, out_vars, step = 10, file = None):
    ens = ens[::step]
    ens_out = np.empty((len(ens), len(out_vars)))
    for i in range(len(ens)):
        ens_i = ens[i]
        net.set_var_vals(ens_i)
        ens_out_i = [net.get_var_val(var) for var in out_vars]
        ens_out[i,] = np.log10(ens_out_i)
    fig, ax = plt.subplots()
    cax = ax.matshow(np.corrcoef(ens_out.transpose()), cmap = "RdYlBu", vmin = -1, vmax = 1)
    ax.set_xticks(range(len(out_vars)));
    ax.set_xticklabels(out_vars, rotation = 90)
    ax.set_yticks(range(len(out_vars)));
    ax.set_yticklabels(out_vars)
    fig.colorbar(cax)

def fit_exps(exps, nets, params_fixed, params_constrained, params_free, global_fit = True, exp_ids=None, local_it=20, global_it=10000, return_ens=True):
    """ 
    Fit a set of experiments in a SloppyCell model.

    This function generates a SloppyCell model from a set of experiments
    and nets and performs local and global optimization given sets of fixed,
    constrained, and free parameters.
  
    Parameters:
    -----------
    exps : list
        A list of SloppyCell experiment objects.
    nets : list
        A list of SloppyCell network objects
    params_fixed : KeyedList
        A keyed list of fixed parameters
    params_constrained : KeyedList
        A keyed list of constrained parameters
    params_free : KeyedList
        A keyed list of free parameters
    global_fit : bool, optional
        When true, perform ensemble analysis.
    exp_ids : list, optional
        A list of ids to restrict the fit to a subset of experiments.
    local_it : int, optional
        Number of iterations for the local fit.
    global_it : int
        Number of iterations for the ensemble analysis.
    return_ens : bool
        When true, return the parameter ensemble.
  
    Returns:
    --------
    out : tuple of Model and KeyedList and, optionally, two lists
        If `return_ens` is false (the default case), return the SloppyCell 
        model and the optimal parameter set.

        Otherwise return also the parameter ensemble and the corresponding
        costs.
  
    """
    if exp_ids is None:
        exp_ids = exps.keys()

    exp_set = [exps[exp_id] for exp_id in exp_ids]
    net_set = [item for sublist in [nets[key].values() for key in exp_ids] for item in sublist]

    keys_opt = list(set(params_constrained.keys()).union(set(params_free.keys())))
    
    for net in net_set:
        keys_non_opt = list(set(net.optimizableVars.keys()) - set(keys_opt))    
        for key in keys_non_opt:
            net.set_var_constant(key, is_constant=True)
            net.set_var_optimizable(key, is_optimizable=False)
        for key in keys_opt:
            net.set_var_optimizable(key, True) 
        net.set_var_ics(params_fixed)
        net.set_var_ics(params_constrained)
        net.set_var_ics(params_free)

    m = Model(exp_set, net_set)

    for key, val in params_free.items():
        res = Residuals.PriorInLog(key+'_prior', key, np.log(val), np.log(np.sqrt(10)))
        m.AddResidual(res)

    for key, val in params_constrained.items():
        res = Residuals.PriorInLog(key+'_prior', key, np.log(val), np.log(np.sqrt(1.1)))
        m.AddResidual(res)

    print "Performing local optimization ..."
    params_opt = Optimization.fmin_lm_log_params(m, params=m.params, maxiter=local_it, disp=False)
    print " done.\n"


    if global_fit:
        print "Performing ensemble analysis ..."
        Network.full_speed()
        gs_start = np.Inf
        gs_end = 0
        params = params_opt
    
        while np.abs(gs_start-gs_end)>1:
            gs_start = gs_end
            j = m.jacobian_log_params_sens(np.log(params))
            jtj = np.dot(np.transpose(j), j)
            ens, gs, r, sf = Ensembles.ensemble_log_params(m, params, jtj, steps=np.round(global_it/10), temperature=1, save_scalefactors=True)
            gs_end = np.min(gs)
            params = ens[np.argmin(gs)]

        j = m.jacobian_log_params_sens(np.log(params))
        jtj = np.dot(np.transpose(j), j)
        ens, gs, r, sf = Ensembles.ensemble_log_params(m, params, jtj, steps=global_it, temperature=1, save_scalefactors=True) 
        params_opt = ens[np.argmin(gs)]
        print " done.\n"

    m.params.update(params_opt)
    m.CalculateForAllDataPoints(params_opt);

    if return_ens:
        if not global_fit:
            ens = []
            gs = []
        return m, params_opt, ens, gs
    else:
        return m, params_opt


