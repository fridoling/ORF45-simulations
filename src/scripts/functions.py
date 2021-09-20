## Functions used by other scripts

# import required modules
import numpy as np
import copy
import pandas as pd

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
 
