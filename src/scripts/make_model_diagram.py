import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import networkx as nx
import os
import pickle
import pandas as pd


## Diagram of reaction scheme
fig, ax = plt.subplots(1,1, figsize=(3,2))
ax.clear()
ax.axis('off');
DG = nx.DiGraph()
dx = 1
dy = 0.5
node_list = ['x', 'y', 'xy']
node_colors = ['lightgray',
               'lightgray',
               'lightgray']
node_labels = {
    'x': 'x',
    'y': 'y',
    'xy': 'xy',
}
edge_labels = {
    ('x', 'xy'): 'KD_xy',
    ('y', 'xy'): 'KD_xy'
}


nodes = {
    'x': (-dx,dy),
    'y': (-dx,-dy),
    'xy': (0,0)
}

edges1 = [('x', 'xy'), ('y', 'xy')]
edges_r = [('xy', 'x'), ('xy', 'y')]


node_sizes = [500, 500, 1000]
node_sizes_f = [500, 500]
nx.draw_networkx_nodes(DG, pos=nodes, nodelist = node_list, ax=ax, edgecolors='darkgrey', node_color=node_colors, node_size = node_sizes, node_shape='o', alpha=1, zorder=1);
nx.draw_networkx_labels(DG, pos=nodes, labels=node_labels, ax=ax, font_size=8,  font_family='monospace');
nx.draw_networkx_edges(DG, pos=nodes, node_size=1000, edgelist=edges1, alpha=1, ax=ax, zorder=2, width=1.5);
nx.draw_networkx_edges(DG, pos=nodes, node_size=500, edgelist=edges_r, alpha=1, ax=ax, zorder=2, width=1.5);
nx.draw_networkx_edge_labels(DG, pos=nodes, edge_labels=edge_labels, font_size=8,  font_family='monospace', label_pos=0.6);

plt.tight_layout(pad = 0.0)
plt.savefig('../../res/diagram_reaction_scheme.pdf')
plt.savefig('../../res/diagram_reaction_scheme.eps')


## Diagrams of catalytic reactions
dx = 1
dy = 0.5
node_list = ['sub', 'enz1', 'subenz', 'enz2', 'prod']
node_colors = {
        "EK": ['mistyrose', 'yellow', 'orange', 'yellow', 'lightyellow'],
        "EP": ['lightyellow', 'cyan', 'dodgerblue', 'cyan', 'mistyrose'],
        "pER": ['mistyrose', 'lightyellow', 'lightyellow', 'lightyellow', 'lavender'],
        "RP2": ['lavender', 'cyan', 'dodgerblue', 'cyan', 'mistyrose']
        }

node_labels = {
        "EK": {'sub': 'E', 'enz1': 'K', 'subenz': 'EK', 'prod': 'pE', 'enz2': 'K'},
        "EP": {'sub': 'pE', 'enz1': 'P', 'subenz': 'pEP', 'prod': 'E', 'enz2': 'P'},
        "pER": {'sub': 'R', 'enz1': 'pE', 'subenz': 'pER', 'prod': 'pR', 'enz2': 'pE'},
        "RP2": {'sub': 'pR', 'enz1': 'P2', 'subenz': 'pRP2', 'prod': 'R', 'enz2': 'P2'}
}
edge_labels = {
        "EK": {('sub', 'subenz'): 'KD_EK', ('enz1', 'subenz'): 'KD_EK', ('subenz', 'enz2'): 'kp_E', ('subenz', 'prod'): 'kp_E'},
        "EP": {('sub', 'subenz'): 'KD_EP', ('enz1', 'subenz'): 'KD_EP', ('subenz', 'enz2'): 'kdp_E', ('subenz', 'prod'): 'kdp_E'},
        "pER": {('sub', 'subenz'): 'KD_ER', ('enz1', 'subenz'): 'KD_ER', ('subenz', 'enz2'): 'kp_R', ('subenz', 'prod'): 'kp_R'},
        "RP2": {('sub', 'subenz'): 'KD_RP2', ('enz1', 'subenz'): 'KD_RP2', ('subenz', 'enz2'): 'kdp_R', ('subenz', 'prod'): 'kdp_R'}
}


nodes = {
    'sub': (-dx,-dy),
    'enz1': (-dx,dy),
    'subenz': (0,0),
    'prod': (dx,-dy),
    'enz2': (dx,dy),    
}

edges1 = [('sub', 'subenz'), ('enz1', 'subenz')]
edges2 = [('subenz', 'prod'), ('subenz', 'enz2')]
edges_r = [('subenz', 'sub'), ('subenz', 'enz1')]


node_sizes = [500, 500, 1000, 500, 500]
node_sizes_f = [500, 500, 1000, 500, 500]

for key in ['EK', 'EP', 'pER', 'RP2']:
    fig, ax = plt.subplots(1,1, figsize=(4,2))
    ax.clear()
    ax.axis('off');
    DG = nx.DiGraph()
    
    nx.draw_networkx_nodes(DG, pos=nodes, nodelist = node_list, ax=ax, edgecolors='darkgrey', node_color=node_colors[key], node_size = node_sizes, node_shape='o', alpha=1, zorder=1);
    nx.draw_networkx_labels(DG, pos=nodes, labels=node_labels[key], ax=ax, font_size=8,  font_family='monospace');
    nx.draw_networkx_edges(DG, pos=nodes, node_size=1000, edgelist=edges1, alpha=1, ax=ax, zorder=2, width=1.5);
    nx.draw_networkx_edges(DG, pos=nodes, node_size=500, edgelist=edges2, alpha=1, ax=ax, zorder=2, width=1.5);
    nx.draw_networkx_edges(DG, pos=nodes, node_size=500, edgelist=edges_r, alpha=1, ax=ax, zorder=2, width=1.5);
    nx.draw_networkx_edge_labels(DG, pos=nodes, edge_labels=edge_labels[key], font_size=8,  font_family='monospace');
    
    plt.tight_layout(pad = 0.0)
    plt.savefig('../../res/diagram_'+key+'.pdf')
    plt.savefig('../../res/diagram_'+key+'.eps')



## Diagram of full model
nodes = {}
node_pos = {}
nodes[0] = ['O']
nodes[1] = ['R', 'pR', 'E', 'pE']
nodes[2] = ['OpR', 'pEO', 'OR', 'EO', 'EpR', 'pEpR', 'pER', 'ER']
nodes[3] = ['EO_OR', 'ER_OR', 'EO_ER', 'EO_EpR', 'EO_OpR', 'EpR_OpR', 'pEpR_OpR', 'pEO_OpR', 'pEO_pEpR', 'pEO_pER', 'pEO_OR', 'pER_OR']
nodes[4] = ['EO_EpR_OpR', 'pEO_pEpR_OpR', 'pEO_pER_OR', 'EO_ER_OR']
d = [1, 2.8, 5, 8]
node_pos[0] = {'O': (0,0)}
node_pos[1] = {'R': (-d[0],0),
            'E': (0,d[0]),
            'pR': (d[0],0),
            'pE': (0,-d[0])
          }
node_pos[2] = {'EO': (0,d[1]),
            'EpR': (d[1]*np.sin(np.pi/4.0), d[1]*np.cos(np.pi/4.0)),
            'OpR': (d[1],0),
            'pEpR': (d[1]*np.sin(3*np.pi/4.0), d[1]*np.cos(3*np.pi/4.0)),
            'pEO': (0,-d[1]),
            'pER': (d[1]*np.sin(5*np.pi/4.0), d[1]*np.cos(5*np.pi/4.0)),
            'OR': (-d[1],0),
            'ER': (d[1]*np.sin(7*np.pi/4.0), d[1]*np.cos(7*np.pi/4.0))
            }
node_pos[3] = {'EO_EpR': (d[2]*np.sin(np.pi/12.0), d[2]*np.cos(np.pi/12.0)),
             'EO_OpR': (d[2]*np.sin(3*np.pi/12.0), d[2]*np.cos(3*np.pi/12.0)),
             'EpR_OpR': (d[2]*np.sin(5*np.pi/12.0), d[2]*np.cos(5*np.pi/12.0)),
             'pEpR_OpR': (d[2]*np.sin(7*np.pi/12.0), d[2]*np.cos(7*np.pi/12.0)),
             'pEO_OpR': (d[2]*np.sin(9*np.pi/12.0), d[2]*np.cos(9*np.pi/12.0)),
             'pEO_pEpR': (d[2]*np.sin(11*np.pi/12.0), d[2]*np.cos(11*np.pi/12.0)),
             'pEO_pER': (d[2]*np.sin(13*np.pi/12.0), d[2]*np.cos(13*np.pi/12.0)),
             'pEO_OR': (d[2]*np.sin(15*np.pi/12.0), d[2]*np.cos(15*np.pi/12.0)),
             'pER_OR': (d[2]*np.sin(17*np.pi/12.0), d[2]*np.cos(17*np.pi/12.0)),
             'ER_OR': (d[2]*np.sin(19*np.pi/12.0), d[2]*np.cos(19*np.pi/12.0)),
             'EO_OR': (d[2]*np.sin(21*np.pi/12.0), d[2]*np.cos(21*np.pi/12.0)),
             'EO_ER': (d[2]*np.sin(23*np.pi/12.0), d[2]*np.cos(23*np.pi/12.0)),
            }
node_pos[4] = {'EO_EpR_OpR': (d[3]*np.sin(np.pi/4.0), d[3]*np.cos(np.pi/4.0)),
             'pEO_pEpR_OpR': (d[3]*np.sin(3*np.pi/4.0), d[3]*np.cos(3*np.pi/4.0)),
             'pEO_pER_OR': (d[3]*np.sin(5*np.pi/4.0), d[3]*np.cos(5*np.pi/4.0)),
             'EO_ER_OR': (d[3]*np.sin(7*np.pi/4.0), d[3]*np.cos(7*np.pi/4.0)),
}

edge_labels = {}
edge_labels[0] = {
    ('O', 'EO'): 'KD_EO',
    ('O', 'OpR'): 'KD_OpR',
    ('O', 'pEO'): 'KD_pEO',
    ('O', 'OR'): 'KD_OR'
}

edge_labels[1] = {
    ('E', 'ER'): 'KD_ER',
    ('E', 'EpR'): 'KD_EpR',
    ('pE', 'pER'): 'KD_ER',
    ('pE', 'pEpR'): 'KD_EpR',
    ('R', 'ER'): 'KD_ER',
    ('pR', 'EpR'): 'KD_EpR',
    ('R', 'pER'): 'KD_ER',
    ('pR', 'pEpR'): 'KD_EpR'
}
edge_labels[2] = {
    ('EO', 'EO_OR'): 'KD_OR',
    ('EO', 'EO_OpR'): 'KD_OpR',
    ('OR', 'EO_OR'): 'KD_EO',
    ('OR', 'pEO_OR'): 'KD_pEO',
    ('pEO', 'pEO_OR'): 'KD_OR',
    ('pEO', 'pEO_OpR'): 'KD_OpR',
    ('OpR', 'pEO_OpR'): 'KD_pEO',
    ('OpR', 'EO_OpR'): 'KD_EO'
}
edge_labels[3] = {
    ('EO', 'EO_ER'): 'KD_ER',
    ('EO', 'EO_EpR'): 'KD_EpR',
    ('OR', 'ER_OR'): 'KD_ER',
    ('OR', 'pER_OR'): 'KD_ER',
    ('pEO', 'pEO_pER'): 'KD_ER',
    ('pEO', 'pEO_pEpR'): 'KD_EpR',
    ('OpR', 'pEpR_OpR'): 'KD_EpR',
    ('OpR', 'EpR_OpR'): 'KD_EpR'
}
edge_labels[4] = {
    ('ER', 'EO_ER'): 'KD_EO',
    ('ER', 'ER_OR'): 'KD_OR',
    ('pER', 'pER_OR'): 'KD_OR',
    ('pER', 'pEO_pER'): 'KD_pEO',
    ('pEpR', 'pEO_pEpR'): 'KD_pEO',
    ('pEpR', 'pEpR_OpR'): 'KD_OpR',
    ('EpR', 'EpR_OpR'): 'KD_OpR',
    ('EpR', 'EO_EpR'): 'KD_EO',
}
edge_labels[5] = {
    ('EO_ER', 'EO_ER_OR'): 'KD_OR',
    ('EO_OR', 'EO_ER_OR'): 'KD_ER',
    ('ER_OR', 'EO_ER_OR'): 'KD_EO',
    ('pER_OR', 'pEO_pER_OR'): 'KD_pEO',
    ('pEO_OR', 'pEO_pER_OR'): 'KD_ER',
    ('pEO_pER', 'pEO_pER_OR'): 'KD_OR',
    ('pEO_pEpR', 'pEO_pEpR_OpR'): 'KD_OpR',
    ('pEO_OpR', 'pEO_pEpR_OpR'): 'KD_EpR',
    ('pEpR_OpR', 'pEO_pEpR_OpR'): 'KD_pEO',
    ('EpR_OpR', 'EO_EpR_OpR'): 'KD_EO',
    ('EO_OpR', 'EO_EpR_OpR'): 'KD_EpR',
    ('EO_EpR', 'EO_EpR_OpR'): 'KD_OpR'
}

edge_labels[6] = {
    ('O', 'EO_ER'): 'KD_EO',
    ('O', 'EO_EpR'): 'KD_EO',
    ('O', 'ER_OR'): 'KD_OR',
    ('O', 'pER_OR'): 'KD_OR',
    ('O', 'pEO_pER'): 'KD_pEO',
    ('O', 'pEO_pEpR'): 'KD_pEO',
    ('O', 'pEpR_OpR'): 'KD_OpR',
    ('O', 'EpR_OpR'): 'KD_OpR'
}

edge_labels[7] = {
    ('E', 'EO_OR'): 'KD_EO',
    ('E', 'EO_OpR'): 'KD_EO',
    ('E', 'EO_OR'): 'KD_EO',
    ('pE', 'pEO_OR'): 'KD_pEO',
    ('pE', 'pEO_OR'): 'KD_pEO',
    ('pE', 'pEO_OpR'): 'KD_pEO',
    ('pE', 'pEO_OpR'): 'KD_pEO',
    ('E', 'EO_OpR'): 'KD_EO',

    ('R', 'EO_OR'): 'KD_OR',
    ('pR', 'EO_OpR'): 'KD_OpR',
    ('R', 'EO_OR'): 'KD_OR',
    ('R', 'pEO_OR'): 'KD_OR',
    ('R', 'pEO_OR'): 'KD_OR',
    ('pR', 'pEO_OpR'): 'KD_OpR',
    ('pR', 'pEO_OpR'): 'KD_OpR',
    ('pR', 'EO_OpR'): 'KD_OpR',

    ('R', 'EO_ER'): 'KD_ER',
    ('pR', 'EO_EpR'): 'KD_EpR',
    ('R', 'pEO_pER'): 'KD_ER',
    ('pR', 'pEO_pEpR'): 'KD_EpR',

    ('E', 'ER_OR'): 'KD_ER',
    ('pE', 'pER_OR'): 'KD_ER',
    ('pE', 'pEpR_OpR'): 'KD_EpR',
    ('E', 'EpR_OpR'): 'KD_EpR',

}



edges = {}
edges_r = {}
for i in range(len(edge_labels)):
    edges[i] = edge_labels[i].keys()
    edges_r[i] = [(b,a) for (a,b) in edge_labels[i].keys()]
node_pos_all = {}
for i in range(5):
    node_pos_all.update(node_pos[i])

cols = {}
cols[0] = 'lightgrey'
for i in np.arange(1,5):
    cols[i] = np.empty(len(node_pos[i]), dtype=object)
    for j in range(len(node_pos[i])):
        if 'pR' in nodes[i][j] and not 'pE' in nodes[i][j]:
            cols[i][j] = 'lavender'
        if not 'pR' in nodes[i][j] and 'pE' in nodes[i][j]:
            cols[i][j] = 'lightyellow'
        if not 'pR' in nodes[i][j] and not 'pE' in nodes[i][j]:
            cols[i][j] = 'mistyrose'
        if 'pR' in nodes[i][j] and 'pE' in nodes[i][j]:
            cols[i][j] = 'aquamarine'

node_labels = {}
for i in range(5):
    node_labels[i] = {k: k for k in node_pos[i].keys()}

node_sizes_n = [500, 500, 1000, 1500, 3500]
node_sizes_e = [1000, 1000, 1500, 1500, 1500, 3500, 1500, 1500]
node_sizes_r = [500, 500, 1000, 1000, 1000, 1500, 500, 500]
label_pos = [0.375, 0.5, 0.35, 0.55, 0.45, 0.6, 0.3, 0.3]
node_shapes = ['o', 'o', 'o', 'o', 'o']

edge_cols = {}
cmap = matplotlib.cm.get_cmap('Dark2')
KD_cols = {}
KDs = ['KD_ER', 'KD_EpR', 'KD_OR', 'KD_OpR', 'KD_EO', 'KD_pEO']
for i, KD in zip(range(len(KDs)), KDs):
    KD_cols[KD] = cmap(i)
for i in range(len(edges)):
    edge_cols[i] = []
    for edge in edges[i]:
        KD = edge_labels[i][edge]
        edge_cols[i].append(KD_cols[KD])

from matplotlib.font_manager import FontProperties

fig, ax = plt.subplots(1,1, figsize=(7,7))
ax.clear()
ax.axis('off');
DG = nx.DiGraph()

for j in range(len(node_labels)):
    nx.draw_networkx_nodes(DG, pos=node_pos[j], nodelist=nodes[j], ax=ax, edgecolors='darkgrey', node_color=cols[j], node_size = node_sizes_n[j], node_shape=node_shapes[j], alpha=1, zorder=1);
    nx.draw_networkx_labels(DG, pos=node_pos[j], labels=node_labels[j], ax=ax, font_size=8, font_family='monospace');

for i in range(len(edge_labels)):
    nx.draw_networkx_edges(DG, pos=node_pos_all, node_size=node_sizes_e[i], edgelist=edges[i], alpha=1, ax=ax, zorder=2, width=1.5, edge_color=edge_cols[i]);
    nx.draw_networkx_edges(DG, pos=node_pos_all, node_size=node_sizes_r[i], edgelist=edges_r[i], alpha=1, ax=ax, zorder=2, width=1.5, edge_color=edge_cols[i]);
    nx.draw_networkx_edge_labels(DG, pos=node_pos_all, label_pos=label_pos[i], edge_labels=edge_labels[i], font_size=5,  font_family='monospace');

plt.tight_layout(pad = 0.0)
plt.savefig('../../res/diagram_all.pdf')
plt.savefig('../../res/diagram_all.eps')
plt.close('all')
