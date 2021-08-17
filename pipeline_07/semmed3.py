import os
import pickle
import pandas as pd
import seaborn as sns
import re
from collections import defaultdict
from tqdm import tqdm


def drop_edges(df, drop_edges):
    #idx = df.query('PREDICATE in @drop_edges').index
    idx = df.query('PREDICATE not in @drop_edges').index
    df.drop(idx, inplace=True)

edges = pd.read_csv("edges4.csv", index_col=0)
nodes = pd.read_csv("nodes_blm.csv", index_col=0)
mappings= pd.read_csv("semmed_mapping.csv",index_col=False)
print(mappings)

support_predicate=[]
for item in mappings['predicate']:
    support_predicate.append(item)

#remove edges with no nodes
print(len(edges))
edges = edges[edges.SUBJECT_CUI.isin(nodes.ID) & edges.OBJECT_CUI.isin(nodes.ID)]
print(len(edges))

#vc = edges.PREDICATE.value_counts()
#print(vc)

print("before drop:",len(edges))
remove_edges = ['compared_with', 'higher_than', 'lower_than', 'different_from', 'different_than', 'same_as',
               'OCCURS_IN', 'PROCESS_OF', 'DIAGNOSES', 'METHOD_OF', 'USES',
               'AUGMENTS', 'ADMINISTERED_TO', 'COMPLICATES']


### drop edge/predicate not in yaml file
#drop_edges(edges, remove_edges)
drop_edges(edges, support_predicate)
print("after drop:",len(edges))


### remove nodes with no edges
print(len(nodes))
nodes = nodes[nodes.ID.isin(set(list(edges.SUBJECT_CUI) + list(edges.OBJECT_CUI)))]
print(len(nodes))


edges.to_csv('edges_filtered.csv', index=False)
nodes.to_csv("nodes_filtered.csv", index=False)
