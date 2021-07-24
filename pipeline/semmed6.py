import os
import pickle
import pandas as pd


nodes = pd.read_csv("nodes_xref.csv")

nodes.ID = "UMLS:" + nodes.ID
nodes['category:STRING'] = nodes.blm_category
nodes['id:STRING'] = nodes.ID
nodes.umls_type_label = nodes.umls_type_label.str.replace("|", ";")
nodes.umls_type = nodes.umls_type.str.replace("|", ";")

nodes.rename(columns = {'ID': ':ID', 
                        'LABEL': 'name:STRING', 
                        'blm_category': ':LABEL',
                        'umls_type_label': 'umls_type_label:STRING[]',
                        'umls_type': 'umls_type:STRING[]',
                        'xrefs': 'xrefs:STRING[]'}, inplace=True)

nodes.to_csv("nodes_neo4j.csv", index=False)

edges = pd.read_csv('edges_biolink.csv')
edges['START_ID'] = "UMLS:" + edges.SUBJECT_CUI
edges['END_ID'] = "UMLS:" + edges.OBJECT_CUI

edges['n_pmids'] = edges.PMID.str.count(";")+1


edges.bl_pred = edges.bl_pred.str.lower()
edges.rename(columns = {'START_ID': ':START_ID', 'END_ID': ':END_ID', 
                        'bl_pred': ':TYPE', 'NEG': 'negated',
                       'PMID': 'pmids'}, inplace=True)

edges['is_defined_by'] = "semmeddb"
edges['relation'] = "semmeddb:" + edges[":TYPE"].str.lower()
edges['provided_by'] = "semmeddb_sulab"

del edges['SUBJECT_CUI']
del edges['OBJECT_CUI']

edges.to_csv("edges_neo4j.csv", index=False)

print(edges.head())
