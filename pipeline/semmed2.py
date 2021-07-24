import os
import pickle
import pandas as pd
import seaborn as sns
import shelve
from itertools import chain
import re
from collections import defaultdict, Counter
from tqdm import tqdm
import requests
import json
from numpy import nan

groups = pd.read_csv("https://metamap.nlm.nih.gov/Docs/SemGroups_2018.txt", sep="|",
                    names=['abv', 'group', 'id', 'label'])
print(groups.head())


names = pd.read_csv("https://metamap.nlm.nih.gov/Docs/SemanticTypes_2018AB.txt", sep="|", 
                   names=['abv', 'ID', 'label'])
print(names.head())

type_label = dict(zip(names.ID, names.label))

df = pd.read_csv("MRSTY.RRF.gz", sep="|", names=['ID', 'TYPE','a', 'b', 'c', 'd'], index_col=False, usecols=['ID', 'TYPE'], dtype=str)

id_type = df.groupby("ID").TYPE.aggregate(set).to_dict()
id_type_label = {k:{type_label.get(x) for x in v} for k,v in id_type.items()}

nodes = pd.read_csv("nodes1.csv", index_col=0)
nodes['umls_type'] = nodes.ID.map(lambda x: id_type.get(x))
nodes['umls_type_label'] = nodes.umls_type.map(lambda x:{type_label.get(y) for y in x} if x else None)

nodes.dropna(subset=['umls_type'], inplace=True)

blm_to_umls = json.load(open("blm_to_umls_nodes.json"))
blm_to_umls = {k:set(v) for k,v in blm_to_umls.items()}

umls_to_blm_check = defaultdict(set)
umls_to_blm = dict()
for k,vv in blm_to_umls.items():
    for v in vv:
        umls_to_blm_check[v.lower()].add(k.lower())
        umls_to_blm[v.lower()] = k.lower()
assert set(len(x) for x in umls_to_blm_check.values()) == {1}

nodes['blm_category'] = nodes.umls_type_label.map(lambda x: {umls_to_blm.get(xx.lower()) for xx in x})
nodes.blm_category = nodes.blm_category.map(lambda v: {x for x in v if x})
nodes.blm_category = nodes.blm_category.map(lambda v: v if v else nan)
nodes.dropna(subset=['blm_category'], inplace=True)
print(nodes.head())

nodes[nodes.blm_category.map(len)>1].blm_category.map(frozenset).value_counts()
nodes2 = nodes[nodes.blm_category.map(len)>1]
Counter(nodes2[nodes2.umls_type_label.map(len)>1].umls_type_label.map(frozenset)).most_common(10)

nodes.loc[nodes.blm_category == {'Protein', 'ChemicalSubstance'}, "blm_category"] = {"Protein"}
nodes.loc[nodes.blm_category == {'Cell', 'ChemicalSubstance'}, "blm_category"] = {"Cell"}
nodes.loc[nodes.blm_category == {'GenomicEntity', 'ChemicalSubstance'}, "blm_category"] = {"GenomicEntity"}
nodes.blm_category = nodes.blm_category.map(lambda x:list(x)[0])

nodes.umls_type_label = nodes.umls_type_label.apply("|".join)
nodes.umls_type = nodes.umls_type.apply("|".join)
nodes.to_csv("nodes_blm.csv")
