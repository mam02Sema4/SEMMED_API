import os
import pandas as pd


edges = pd.read_csv("edges4.csv", index_col=0)
nodes = pd.read_csv("nodes_umls.csv", index_col=0)


#remove edges with no nodes
print("total edge count:",len(edges))
edges = edges[edges.SUBJECT_CUI.isin(nodes.ID) & edges.OBJECT_CUI.isin(nodes.ID)]
print("Remove edges with no node info in nodes_umls.csv:",len(edges))


### remove nodes with no related edges
print("total node count:",len(nodes))
nodes = nodes[nodes.ID.isin(set(list(edges.SUBJECT_CUI) + list(edges.OBJECT_CUI)))]
print("drop nodes with no related edges:",len(nodes))

edges.to_csv('edges_filtered_new.csv', index=False)
nodes.to_csv("nodes_filtered_new.csv", index=False)
