from SloppyCell.ReactionNetworks import *
from collections import OrderedDict
import copy
import numpy as np
import functions

### load basic networks
net_basic = IO.from_SBML_file('../models/model_invitro_incell.xml','net_basic')
net_basic_sub = IO.from_SBML_file('../models/model_with_substrates.xml','net_basic_sub')

nets = OrderedDict()
nets_sub = OrderedDict()
exps = OrderedDict()
data = OrderedDict()
nets['basic'] = net_basic
nets_sub['basic'] = net_basic_sub

### set total amounts constant
for key in ['Ptot', 'P2tot', 'Etot', 'Rtot', 'Ktot', 'Otot']:
    net_basic.set_var_constant(key, is_constant=True)
    net_basic.set_var_optimizable(key, is_optimizable=False)

for key in ['Ptot', 'P2tot', 'Etot', 'Rtot', 'Ktot', 'Otot', 'SUBFtot', 'SUBCtot']:
    net_basic_sub.set_var_constant(key, is_constant=True)
    net_basic_sub.set_var_optimizable(key, is_optimizable=False)

### IN VITRO EXPERIMENTS
### Experiment 1: PERK-RSK
# define net in which phosphorylation is switched on after t = 1000
nets['PERK_RSK'] = OrderedDict()
nets['PERK_RSK']['K'] = net_basic.copy(new_id='PERK_RSK')
nets['PERK_RSK']['K'].set_var_ic('Etot', 0.01)
nets['PERK_RSK']['K'].set_var_ic('pE', 0.01)
nets['PERK_RSK']['K'].set_var_ic('Ptot', 0)
nets['PERK_RSK']['K'].set_var_ic('P2tot', 0)
nets['PERK_RSK']['K'].set_var_ic('Rtot', 2)
nets['PERK_RSK']['K'].set_var_ic('Ktot', 0)
nets['PERK_RSK']['K'].set_var_ic('Otot', 0)
nets['PERK_RSK']['K'].set_var_ic('kp_R_on', 0)
nets['PERK_RSK']['K'].set_var_constant('kp_R_on', False)
nets['PERK_RSK']['K'].add_event(id='switch_kp_R_on', trigger='gt(time, 1000)', event_assignments={'kp_R_on': 1})
nets['PERK_RSK']['K'].set_name('K')

# define net with same switch, but without O
nets['PERK_RSK']['+ORF45'] = nets['PERK_RSK']['K'].copy(new_id='PERK_RSK_ORF')
nets['PERK_RSK']['+ORF45'].set_var_ic('Otot', 50)
nets['PERK_RSK']['+ORF45'].set_name('+ORF45')

data['PERK_RSK'] = {
    'PERK_RSK': {
        'pRtot': {
            1000: (0, 0.1),
            1001: (18.74, 1),
            1002: (27.71, 1),
            1003: (37.69, 1),
            1004: (44.51, 1),
            1006: (60.14, 1),
            1008.5: (75.24, 1),
            1012: (75.96, 1)
        }
    },
    'PERK_RSK_ORF': {
        'pRtot': {
            1000: (0, 0.1),
            1001: (2.9, 1),
            1002: (5.24, 1),
            1003: (9.17, 1),
            1004: (10.61, 1),
            1006: (13.33, 1),
            1008.5: (17.32, 1),
            1012: (19.32, 1)
        }
    }    
}
data['PERK_RSK'] = functions.set_errors(data['PERK_RSK'], error=0.05, absolute=False)
exp = Experiment(name='PERK_RSK', data=data['PERK_RSK'])
#exp.set_sf_prior([('pRtot')], 'gaussian in log sf', (np.log(38), np.log(2)))
exp.set_fixed_sf({'pRtot': 38})
exps['PERK_RSK'] = exp


### Experiment 2: MKK1-ERK-RSK

# define net in which phosphorylation is switched on after t = 1000
nets['MKK1_ERK_RSK'] = OrderedDict()
nets['MKK1_ERK_RSK']['K'] = net_basic.copy(new_id='MKK1_ERK_RSK')
nets['MKK1_ERK_RSK']['K'].set_var_ic('Etot', 2)
nets['MKK1_ERK_RSK']['K'].set_var_ic('pE', 0)
nets['MKK1_ERK_RSK']['K'].set_var_ic('Ptot', 0)
nets['MKK1_ERK_RSK']['K'].set_var_ic('P2tot', 0)
nets['MKK1_ERK_RSK']['K'].set_var_ic('Ktot', 0.25)
nets['MKK1_ERK_RSK']['K'].set_var_ic('Rtot', 2)
nets['MKK1_ERK_RSK']['K'].set_var_ic('Otot', 0)
nets['MKK1_ERK_RSK']['K'].set_var_ic('kp_R_on', 0)
nets['MKK1_ERK_RSK']['K'].set_var_constant('kp_R_on', False)
nets['MKK1_ERK_RSK']['K'].add_event(id='switch_kp_R_on', trigger='gt(time, 1000)', event_assignments={'kp_R_on': 1})
nets['MKK1_ERK_RSK']['K'].set_var_ic('kp_E_on', 0)
nets['MKK1_ERK_RSK']['K'].set_var_constant('kp_E_on', False)
nets['MKK1_ERK_RSK']['K'].add_event(id='switch_kp_E_on', trigger='gt(time, 1000)', event_assignments={'kp_E_on': 1})
nets['MKK1_ERK_RSK']['K'].set_name('K')

# define net with same switch, but without O
nets['MKK1_ERK_RSK']['+ORF45'] = nets['MKK1_ERK_RSK']['K'].copy(new_id='MKK1_ERK_RSK_ORF')
nets['MKK1_ERK_RSK']['+ORF45'].set_var_ic('Otot', 50)
nets['MKK1_ERK_RSK']['+ORF45'].set_name('+ORF45')

### data for MKK1-ERK-RSK experiment
data['MKK1_ERK_RSK'] = {
    'MKK1_ERK_RSK': {
        'pRtot': {
            1000: (0, 0.1),
            1001: (0.04, 0.1),
            1002: (0.22, 0.1),
            1003.5: (0.83, 0.1),
            1005: (2.11, 1),
            1007: (4.56, 1),
            1010: (7.59, 1),
            1015: (8.79, 1)
        }
    },
    'MKK1_ERK_RSK_ORF': {
        'pRtot': {
            1000: (0, 0.1),
            1001: (0.01, 0.1),
            1002: (0.02, 0.1),
            1003.5: (0.03, 0.1),
            1005: (0.06, 0.1),
            1007: (0.12, 0.1),
            1010: (0.21, 0.1),
            1015: (0.43, 0.1)
        }
    }    
}
data['MKK1_ERK_RSK'] = functions.set_errors(data['MKK1_ERK_RSK'], error=0.05, absolute=False)
exp = Experiment(name='MKK1_ERK_RSK', data=data['MKK1_ERK_RSK'])
exp.set_fixed_sf({'pRtot': 4.5})
#exp.set_sf_prior([('pRtot')], 'gaussian in log sf', (np.log(4.5), np.log(1.1)))
exps['MKK1_ERK_RSK'] = exp


### Experiment 3: PERK+PTASE

### nets for PERK+PTASE experiment
# define net in which phosphatase is switched on after t = 1000
nets['PERK_PTASE'] = OrderedDict()
nets['PERK_PTASE']['K'] = net_basic.copy(new_id='PERK_PTASE')
nets['PERK_PTASE']['K'].set_var_ic('Etot', 1)
nets['PERK_PTASE']['K'].set_var_ic('pE', 1)
nets['PERK_PTASE']['K'].set_var_ic('Ktot', 0)
nets['PERK_PTASE']['K'].set_var_ic('Rtot', 0)
nets['PERK_PTASE']['K'].set_var_ic('Otot', 0)
nets['PERK_PTASE']['K'].set_var_ic('Ptot', 0)
nets['PERK_PTASE']['K'].set_var_ic('P2tot', 0)
nets['PERK_PTASE']['K'].set_var_ic('kp_R_on', 0)
nets['PERK_PTASE']['K'].set_var_ic('kdp_E_on', 1)
nets['PERK_PTASE']['K'].set_var_constant('Ptot', False)
nets['PERK_PTASE']['K'].add_event(id='switch_kdp_E_on', trigger='gt(time, 1000)', event_assignments={'Ptot': 1})
nets['PERK_PTASE']['K'].set_name('K')

# define net with same switch, but with RSK = 1
nets['PERK_PTASE']['+RSK'] = nets['PERK_PTASE']['K'].copy(new_id='PERK_PTASE_RSK')
nets['PERK_PTASE']['+RSK'].set_var_ic('Rtot', 1)
nets['PERK_PTASE']['+RSK'].set_name('+RSK')

# define net with same switch, but with RSK = 1 and o = 1
nets['PERK_PTASE']['+RSK+ORF(1)'] = nets['PERK_PTASE']['+RSK'].copy(new_id='PERK_PTASE_ORF1')
nets['PERK_PTASE']['+RSK+ORF(1)'].set_var_ic('Otot', 1)
nets['PERK_PTASE']['+RSK+ORF(1)'].set_name('+RSK+ORF(1)')

# define net with same switch, but with RSK = 1 and o = 5
nets['PERK_PTASE']['+RSK+ORF(5)'] = nets['PERK_PTASE']['+RSK'].copy(new_id='PERK_PTASE_ORF5')
nets['PERK_PTASE']['+RSK+ORF(5)'].set_var_ic('Otot', 5)
nets['PERK_PTASE']['+RSK+ORF(5)'].set_name('+RSK+ORF(5)')


### data for PERK+PTASE experiment
data['PERK_PTASE'] = {
    'PERK_PTASE': {
        'pEtot': {
            1000: (1, 0.1),
            1001: (0.3, 0.1),
            1003: (0.25, 0.1),
            1010: (0.21, 0.1),
            1020: (0.18, 0.1)
        }
    },
    'PERK_PTASE_RSK': {
        'pEtot': {            
            1000: (1, 0.1),
            1001: (0.46, 0.1),
            1003: (0.22, 0.1),
            1010: (0.17, 0.1),
            1020: (0.14, 0.1)
        }
    },
    'PERK_PTASE_ORF1': {
        'pEtot': {
            1000: (1, 0.1),
            1001: (0.91, 0.1),
            1003: (0.87, 0.1),
            1010: (0.41, 0.1),
            1020: (0.25, 0.1)
        }
    },    
    'PERK_PTASE_ORF5': {
        'pEtot': {
            1000: (1, 0.1),
            1001: (1.07, 0.1),
            1003: (0.89, 0.1),
            1010: (0.85, 0.1),
            1020: (0.59, 0.1)
        }
    }
}

data['PERK_PTASE'] = functions.subtract_bg(data['PERK_PTASE'], bg=0.14)
data['PERK_PTASE'] = functions.set_errors(data['PERK_PTASE'], error=0.05, min_error=0.1, absolute=False)
exp = Experiment(name='PERK_PTASE', data=data['PERK_PTASE'])
exp.set_fixed_sf({'pEtot': 1})
exps['PERK_PTASE'] = exp



#### IN CELL EXPERIMENTS
net_incell = net_basic.copy(new_id='net_incell')
net_incell_sub = net_basic_sub.copy(new_id='net_incell_sub')
for net in [net_incell, net_incell_sub]:
    for key in net_basic.optimizableVars.keys():
        net.set_var_optimizable(key, is_optimizable=True)
    net.set_var_optimizable('kon_RP2', is_optimizable=True)
    net.set_var_optimizable('KD_RP2', is_optimizable=True)
    net.set_var_optimizable('kdp_R', is_optimizable=True)  
    net.set_var_optimizable('kp_K_bg', is_optimizable=True)
    net.set_var_optimizable('kp_K_egf', is_optimizable=True)
    net.set_var_optimizable('kdp_K', is_optimizable=True)

nets['INCELL_WT'] = OrderedDict()
nets['INCELL_WT']['K'] = net_incell.copy(new_id='WT')
nets['INCELL_WT']['K'].set_var_ic('Ktot', 1.2)
nets['INCELL_WT']['K'].set_var_ic('K', 1.2)
nets['INCELL_WT']['K'].set_var_ic('Etot', 0.7)
nets['INCELL_WT']['K'].set_var_ic('Rtot', 2.0)
nets['INCELL_WT']['K'].set_var_ic('Rtot_factor', 1)
nets['INCELL_WT']['K'].set_var_ic('Otot_factor', 0)
nets['INCELL_WT']['K'].set_var_ic('Ptot', 1)
nets['INCELL_WT']['K'].set_var_ic('P2tot', 1)
nets['INCELL_WT']['K'].set_var_ic('kp_R_on', 1)
nets['INCELL_WT']['K'].set_var_ic('kp_K_on', 0)
nets['INCELL_WT']['K'].set_var_constant('kp_K_on', False)
nets['INCELL_WT']['K'].add_event(id='switch_kin_on', trigger='gt(time, 1000)', event_assignments={'kp_K_on': 1})
nets['INCELL_WT']['K'].set_name('K')


nets['INCELL_WT']['ORF'] = nets['INCELL_WT']['K'].copy(new_id='WT_ORF')
nets['INCELL_WT']['ORF'].set_var_ic('Otot_factor', 1)
nets['INCELL_WT']['ORF'].set_var_ic('Otot', 1.0)
nets['INCELL_WT']['ORF'].set_name('+ORF')

nets_sub['INCELL_WT'] = OrderedDict()
nets_sub['INCELL_WT']['K'] = net_incell_sub.copy(new_id='WT')
nets_sub['INCELL_WT']['K'].set_var_ic('Ktot', 1.2)
nets_sub['INCELL_WT']['K'].set_var_ic('K', 1.2)
nets_sub['INCELL_WT']['K'].set_var_ic('Etot', 0.7)
nets_sub['INCELL_WT']['K'].set_var_ic('Rtot', 2.0)
nets_sub['INCELL_WT']['K'].set_var_ic('Rtot_factor', 1)
nets_sub['INCELL_WT']['K'].set_var_ic('Otot_factor', 0)
nets_sub['INCELL_WT']['K'].set_var_ic('Ptot', 1)
nets_sub['INCELL_WT']['K'].set_var_ic('P2tot', 1)
nets_sub['INCELL_WT']['K'].set_var_ic('kp_R_on', 1)
nets_sub['INCELL_WT']['K'].set_var_ic('kp_K_on', 0)
nets_sub['INCELL_WT']['K'].set_var_constant('kp_K_on', False)
nets_sub['INCELL_WT']['K'].add_event(id='switch_kin_on', trigger='gt(time, 1000)', event_assignments={'kp_K_on': 1})
nets_sub['INCELL_WT']['K'].set_name('K')

nets_sub['INCELL_WT']['ORF'] = nets_sub['INCELL_WT']['K'].copy(new_id='WT_ORF')
nets_sub['INCELL_WT']['ORF'].set_var_ic('Otot_factor', 1)
nets_sub['INCELL_WT']['ORF'].set_var_ic('Otot', 1.0)
nets_sub['INCELL_WT']['ORF'].set_name('+ORF')


nets['INCELL_RSK_KO'] = OrderedDict()
nets['INCELL_RSK_KO']['K'] = nets['INCELL_WT']['K'].copy(new_id='RSK_KO')
nets['INCELL_RSK_KO']['K'].set_var_ic('Rtot_factor', 0)
nets['INCELL_RSK_KO']['K'].set_name('K')

nets['INCELL_RSK_KO']['ORF'] = nets['INCELL_WT']['ORF'].copy(new_id='RSK_KO_ORF')
nets['INCELL_RSK_KO']['ORF'].set_var_ic('Rtot_factor', 0)
nets['INCELL_RSK_KO']['ORF'].set_name('+ORF')

### data for PERK+PTASE experiment
data['INCELL_WT'] = {
    'WT': {
        'pEtot': {
            1000: (0.11, 0),
            1001.5: (0.10, 0),
            1003: (0.13, 0),
            1005: (0.21, 0),
            1007: (0.31, 0),
            1010: (0.56, 0),
             1015: (0.96, 0),
#             1020: (1.00, 0),
#             1030: (0.67, 0),
#             1040: (0.58, 0),
#             1050: (0.54, 0),
#             1060: (0.52, 0)            
        },
        'pRtot': {
            1000: (0.02, 0),
            1001.5: (0.21, 0),
            1003: (0.20, 0),
            1005: (0.28, 0),
            1007: (0.27, 0),
            1010:(0.37, 0),
            1015: (0.87, 0),
#             1020: (1.00, 0),
#             1030: (0.83, 0),
#             1040: (0.65, 0),
#             1050: (0.54, 0),
#             1060: (0.27, 0)
        }        
    },
    'WT_ORF': {
        'pEtot': {            
            990: (0.26, 0),
            1000: (0.26, 0),
            1001.5: (0.35, 0),
            1003: (0.62, 0),
            1005: (1.12, 0),
            1007: (1.55, 0),
            1010:(1.69, 0),
            1015: (1.92, 0),
#             1020: (1.67, 0),
#             1030: (1.35, 0),
#             1040: (1.43, 0),
#            1050: (1.45, 0),
#             1060: (1.16, 0)
        },
        'pRtot': {            
            990: (0.65,0),            
            1000: (0.65,0),
            1001.5: (0.73,0),
            1003: (0.76,0),
            1005: (0.84,0),
            1007: (0.91,0),
            1010: (1.08,0),
            1015: (1.73,0),
#             1020: (1.58,0),
#             1030: (1.55,0),
#             1040: (1.49,0),
#             1050: (1.47,0),
#             1060: (1.03,0)          
        }        
    }
}

data['INCELL_WT'] = functions.set_errors(data['INCELL_WT'], error=0.05, min_error=0.01, absolute=True)
exp = Experiment(name='INCELL_WT', data=data['INCELL_WT'])
exp.set_fixed_sf({'pEtot': 9.14})
exps['INCELL_WT'] = exp


### data for PERK+PTASE experiment
data['INCELL_RSK_KO'] = {
    'RSK_KO': {
        'pEtot': {
            1000.0: (0.17,0),
            1002: (0.41,0),
            1003: (0.76,0),
            1005: (0.95,0),
            1007: (0.97,0),
            1010: (0.87,0),
            1015: (0.81,0),
#             1020: (0.68,0),
#             1030: (0.57,0),
#             1040: (0.54,0),
#             1050: (0.53,0),
#             1060: (0.41,0)
        }
    },    
    'RSK_KO_ORF': {
        'pEtot': {
            1000.0: (0.17,0),
            1002: (0.33,0),
            1003: (0.53,0),
            1005: (0.68,0),
            1007: (0.71,0),
            1010: (0.71,0),
            1015: (0.63,0),
#             1020: (0.56,0),
#             1030: (0.55,0),
#             1040: (0.47,0),
#             1050: (0.48,0),
#             1060: (0.35,0)
        }
    }
}

data['INCELL_RSK_KO'] = functions.set_errors(data['INCELL_RSK_KO'], error=0.05, min_error=0.01, absolute=True)
exp = Experiment(name='INCELL_RSK_KO', data=data['INCELL_RSK_KO'])
exp.set_fixed_sf({'pEtot': 9.14})
exps['INCELL_RSK_KO'] = exp
