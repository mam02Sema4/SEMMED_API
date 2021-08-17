import os
import pickle
import pandas as pd
import seaborn as sns
import re
from collections import defaultdict, Counter
from tqdm import tqdm


def is_allowed_edge(domain, pred, rnge):
    allowed_domain, allowed_range = allowed_domain_range[pred] if pred in allowed_domain_range else (None, None)
    return ((domain in allowed_domain if allowed_domain else True) and 
              (rnge in allowed_range if allowed_range else True))

semmed_mappings_file = pd.read_csv("semmed_mapping.csv",index_col=False)
edges = pd.read_csv('edges_filtered.csv')
nodes = pd.read_csv("nodes_filtered.csv")
node_category = dict(zip(nodes.ID, nodes.blm_category))

edges['bl_domain'] = edges.SUBJECT_CUI.apply(node_category.get)
edges['bl_pred'] = edges.PREDICATE.apply(lambda p: p.lower())
edges['bl_range'] = edges.OBJECT_CUI.apply(node_category.get)
edges['bl_type'] = edges['bl_domain'] + "." + edges['bl_pred'] + "." + edges['bl_range']


edges.rename(columns={'PREDICATE': 'SEMMED_PRED'}, inplace=True)

len(set(edges.bl_pred))

allowed_domain_range = {
    'causes': (None, {'biologicalprocessoractivity', 'diseaseorphenotypicfeature'}),
    'location_of': ({'grossanatomicalstructure', 'anatomicalentity', 'cellularcomponent', 'cell'}, None),
    'treats': (None, {'diseaseorphenotypicfeature'}),
    'predisposes': (None, {'diseaseorphenotypicfeature'}),
    'prevents': (None, {'diseaseorphenotypicfeature'}),
}


d = {x: is_allowed_edge(*x.split(".")) for x in set(edges.bl_type)}
allowed_edges = {k for k,v in d.items() if v}

###check ineligible edges
for edge in edges.bl_type:
    if edge not in allowed_edges:
       print ("Not In:",edge)
       print("-------------")

print(len(edges))
idx = edges.bl_type.isin(allowed_edges)
edges = edges[idx]
print(len(edges))

### rename semmed predicate into biolink predicate (only exact mapping now)
#print(edges.bl_pred.value_counts())

for index, row in semmed_mappings_file.iterrows():
    if ((row['biolink_key'].lower()!=row['predicate'].lower()) ):
      #if (row['mapping_type']=='exact_mappings'):
        edges.loc[edges.bl_pred == row['predicate'], "bl_pred"] = row['biolink_key']
        #print("Rename:",row['predicate'],"Into:",row['biolink_key'])
        #print("======")


edges.loc[lambda df:(df['bl_pred'] == "related_to") & (df['bl_domain'] == "Gene") & 
      (df['bl_range'] == "DiseaseOrPhenotypicFeature"), "bl_pred"] = 'gene_associated_with_condition'

nodes = nodes[nodes.ID.isin(set(list(edges['SUBJECT_CUI']) + list(edges['OBJECT_CUI'])))]

edges['bl_type'] = edges['bl_domain'] + "." + edges['bl_pred'] + "." + edges['bl_range']

del edges['bl_type']
del edges['bl_domain']
del edges['bl_range']

edges.to_csv("edges_biolink.csv", index=None)
nodes.to_csv("nodes_biolink.csv", index=None)
