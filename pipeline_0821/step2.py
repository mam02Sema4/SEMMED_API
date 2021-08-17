import os
import pandas as pd
import shelve
from itertools import chain
from collections import defaultdict, Counter
import requests


names = pd.read_csv("https://metamap.nlm.nih.gov/Docs/SemanticTypes_2018AB.txt", sep="|", 
                   names=['abv', 'ID', 'label'])
print(names.head())

type_label = dict(zip(names.ID, names.label))

df = pd.read_csv("MRSTY.RRF.gz", sep="|", names=['ID', 'TYPE','a', 'b', 'c', 'd'], index_col=False, usecols=['ID', 'TYPE'], dtype=str)

id_type = df.groupby("ID").TYPE.aggregate(set).to_dict()
id_type_label = {k:{type_label.get(x) for x in v} for k,v in id_type.items()}

nodes = pd.read_csv("nodes1.csv", index_col=0)
print("Read file - node count:",len(nodes))

nodes['umls_type'] = nodes.ID.map(lambda x: id_type.get(x))
nodes['umls_type_label'] = nodes.umls_type.map(lambda x:{type_label.get(y) for y in x} if x else None)

nodes.dropna(subset=['umls_type'], inplace=True)
print("Drop NA (UMLS type):",len(nodes))

nodes.umls_type_label = nodes.umls_type_label.apply("|".join)
nodes.umls_type = nodes.umls_type.apply("|".join)
nodes.to_csv("nodes_umls.csv")
