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


if not os.path.exists('../data/SPR/SPR_data.pickle'):
    print "Importing data...\n"
    df = {}
    t_iv = {}
    df['RSK_ppERK'] = pd.read_excel('../data/SPR/SPR_data1.xlsx', sheet_name=0)
    df['RSK_ppERK'].columns = ['time', 'RSK_ppERK_012', 'RSK_ppERK_012_fit', 'RSK_ppERK_111', 'RSK_ppERK_111_fit']
    df['RSK_ppERK'] = df['RSK_ppERK'][df['RSK_ppERK']['time']>0]
    ind_max_012 = (df['RSK_ppERK']['RSK_ppERK_012_fit'].diff()<0).idxmax()
    ind_max_111 = (df['RSK_ppERK']['RSK_ppERK_111_fit'].diff()<0).idxmax()
    t_iv['RSK_ppERK_012'] = df['RSK_ppERK']['time'][ind_max_012]
    t_iv['RSK_ppERK_111'] = df['RSK_ppERK']['time'][ind_max_111]
    row_select = np.int(df['RSK_ppERK'].shape[0]/200)
    df['RSK_ppERK'] = df['RSK_ppERK'][::row_select]
    
    
    df['RSK_ppERK_012'] = df['RSK_ppERK'][['time']]
    df['RSK_ppERK_012'].insert(1, 'RSK_ppERK_012', df['RSK_ppERK']['RSK_ppERK_012'].divide(np.max(df['RSK_ppERK']['RSK_ppERK_012_fit'])))
    
    df['RSK_ppERK_111'] = df['RSK_ppERK'][['time']]
    df['RSK_ppERK_111'].insert(1, 'RSK_ppERK_111', df['RSK_ppERK']['RSK_ppERK_111'].divide(np.max(df['RSK_ppERK']['RSK_ppERK_111_fit'])))
    
    
    df['RSK_ORF'] = pd.read_excel('../data/SPR/SPR_data1.xlsx', sheet_name=1, usecols=[0,1,3,6,8])
    df['RSK_ORF'].columns = ['time', 'RSK_ORF_05', 'RSK_ORF_05_fit', 'RSK_ORF_135', 'RSK_ORF_135_fit']
    df['RSK_ORF'] = df['RSK_ORF'][df['RSK_ORF']['time']>0]
    ind_max_05 = (df['RSK_ORF']['RSK_ORF_05_fit'].diff()<0).idxmax()
    ind_max_135 = (df['RSK_ORF']['RSK_ORF_135_fit'].diff()<0).idxmax()
    t_iv['RSK_ORF_05'] = df['RSK_ORF']['time'][ind_max_05]
    t_iv['RSK_ORF_135'] = df['RSK_ORF']['time'][ind_max_135]
    row_select = np.int(df['RSK_ORF'].shape[0]/200)
    df['RSK_ORF'] = df['RSK_ORF'].loc[::row_select]
    df['RSK_ORF'][df['RSK_ORF']<0] = 0
    
    
    df['RSK_ORF_05'] = df['RSK_ORF'][['time']]
    df['RSK_ORF_05'].insert(1, 'RSK_ORF_05', df['RSK_ORF']['RSK_ORF_05'].divide(np.max(df['RSK_ORF']['RSK_ORF_05_fit'])))
    
    df['RSK_ORF_135'] = df['RSK_ORF'][['time']]
    df['RSK_ORF_135'].insert(1, 'RSK_ORF_135', df['RSK_ORF']['RSK_ORF_135'].divide(np.max(df['RSK_ORF']['RSK_ORF_135_fit'])))
    
    
    df['RSK_ppERK_ORF'] = pd.read_excel('../data/SPR/SPR_data1.xlsx', sheet_name=3)
    df['RSK_ppERK_ORF'].columns = ['time', 'RSK_ppERK_ORF_012', 'RSK_ppERK_ORF_012_fit', 'RSK_ppERK_ORF_111', 'RSK_ppERK_ORF_111_fit']
    df['RSK_ppERK_ORF'] = df['RSK_ppERK_ORF'][df['RSK_ppERK_ORF']['time']>0]
    ind_max_012 = np.argmax(df['RSK_ppERK_ORF']['RSK_ppERK_ORF_012_fit'].diff()<0)
    ind_max_111 = np.argmax(df['RSK_ppERK_ORF']['RSK_ppERK_ORF_111_fit'].diff()<0)
    t_iv['RSK_ppERK_ORF_012'] = df['RSK_ppERK_ORF']['time'][ind_max_012]
    t_iv['RSK_ppERK_ORF_111'] = df['RSK_ppERK_ORF']['time'][ind_max_111]
    row_select = np.int(df['RSK_ppERK_ORF'].shape[0]/200)
    df['RSK_ppERK_ORF'] = df['RSK_ppERK_ORF'].loc[::row_select]
    df['RSK_ppERK_ORF'][df['RSK_ppERK_ORF']<0] = 0
    
    
    df['RSK_ppERK_ORF_012'] = df['RSK_ppERK_ORF'][['time']]
    df['RSK_ppERK_ORF_012'].insert(1, 'RSK_ppERK_ORF_012', df['RSK_ppERK_ORF']['RSK_ppERK_ORF_012'].divide(np.max(df['RSK_ppERK_ORF']['RSK_ppERK_ORF_012_fit'])))
    
    df['RSK_ppERK_ORF_111'] = df['RSK_ppERK_ORF'][['time']]
    df['RSK_ppERK_ORF_111'].insert(1, 'RSK_ppERK_ORF_111', df['RSK_ppERK_ORF']['RSK_ppERK_ORF_111'].divide(np.max(df['RSK_ppERK_ORF']['RSK_ppERK_ORF_111_fit'])))
    
    
    df_raw = pd.read_excel('../data/SPR/SPR_data2.xlsx', sheet_name=0, usecols=[0,1,3])
    df_raw.columns = ['time', 'ppERK_ORF', 'ppERK_ORF_fit']
    df_raw = df_raw[df_raw['time']>0]
    ind_max = np.argmax(df_raw['ppERK_ORF'].diff())
    t_iv['ppERK_ORF'] = df_raw['time'][ind_max]/60.0
    
    df['ppERK_ORF'] = df_raw[['time']].divide(60)
    df['ppERK_ORF'].insert(1, 'ppERK_ORF', df_raw['ppERK_ORF'].divide(np.max(df_raw['ppERK_ORF'])))
    row_select = np.int(df['ppERK_ORF'].shape[0]/200)
    df['ppERK_ORF'] = df['ppERK_ORF'].loc[::row_select]
    df['ppERK_ORF'][df['ppERK_ORF']<0] = 0

    with open('../data/SPR/SPR_data.pickle', 'w') as f:
        pickle.dump((df, t_iv), f)
else:
    print "Data has already been processed."
