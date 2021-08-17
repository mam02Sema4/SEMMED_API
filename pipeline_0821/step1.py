import os
import pickle

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
from collections import defaultdict
import itertools
import numpy as np
from numpy import NAN

from sqlalchemy import create_engine
import pandas as pd

db_connection_str = 'mysql+pymysql://acct:pwd@localhost/pubmed'
db_connection = create_engine(db_connection_str,pool_recycle=1)


#dfl = []  

# Create empty dataframe
#dfs = pd.DataFrame()  
#query="SELECT SUBJECT_CUI,PREDICATE,OBJECT_CUI,PMID FROM pubmed.PREDICATION"

# Start Chunking
#for chunk in pd.read_sql(query, con=db_connection,chunksize=1000000):

    # Start Appending Data Chunks from SQL Result set into List
#    dfl.append(chunk)
#    print("New Chunk")
#    print("---------")

# Start appending data from list to dataframe
#dfs = pd.concat(dfl, ignore_index=True)
#dfs.to_csv('SEMMED_result.csv', index=False)
#print("Done!")

cols = ['SUBJECT_CUI','PREDICATE','OBJECT_CUI','PMID']
gb_cols = ['SUBJECT_CUI','PREDICATE','OBJECT_CUI']
sem_df = pd.DataFrame(columns=cols)
df_iter = pd.read_csv('SEMMED_result.csv', dtype=str, chunksize=10000000)
for chunk in tqdm(df_iter, total=10):
    chunk.PMID = chunk.PMID.astype("str")
    c = chunk.groupby(gb_cols).PMID.agg(";".join).reset_index()
    sem_df = sem_df.append(c)

sem_df = sem_df.groupby(gb_cols).PMID.agg(";".join).reset_index()

print(sem_df.head())
sem_df.to_csv("edges1.csv")

multi_start = sem_df['SUBJECT_CUI'].str.contains('|', regex=False)
multi_end = sem_df['OBJECT_CUI'].str.contains('|', regex=False)
pipe_lines = sem_df[multi_start | multi_end].copy()
sem_df = sem_df[~multi_start & ~multi_end]
print('Rows with multiple subjects or objects {:,}'.format(len(pipe_lines)))
print('Rows with only 1 subject AND only 1 object {:,}'.format(len(sem_df)))

pipe_lines.SUBJECT_CUI = pipe_lines.SUBJECT_CUI.str.split('|')
pipe_lines.OBJECT_CUI = pipe_lines.OBJECT_CUI.str.split('|')

# do the combinations
lines = []
for row in tqdm(pipe_lines.itertuples(), total=len(pipe_lines)):
    a = [row.SUBJECT_CUI, row.OBJECT_CUI]
    c = list(itertools.product(*a))
    lines.extend([{'SUBJECT_CUI':x[0], 'PREDICATE':row.PREDICATE, 'OBJECT_CUI':x[1], 'PMID': row.PMID} for x in c])
expanded_df = pd.DataFrame(lines)

sem_df = sem_df.append(expanded_df, sort=True)
sem_df = sem_df.groupby(gb_cols).PMID.agg(";".join).reset_index()


sem_df.replace("", NAN, inplace=True)
print("before drop NAN:",len(sem_df))

# Then drop the rows with NaN
sem_df.dropna(subset=['SUBJECT_CUI','OBJECT_CUI'], inplace=True)
print("after drop NAN:",len(sem_df))

sem_df.to_csv("edges2.csv")

names = list("abcdefghijklmn")
iter_csv = pd.read_csv("MRSAT.RRF.gz", delimiter="|", names=names, index_col=None, chunksize=1000000)
chunks = []
umls_entrez = dict()
for chunk in tqdm(iter_csv, total=67668372/1000000):
    chunk.fillna(method='ffill', inplace=True)
    chunk = chunk[chunk.i == "ENTREZGENE_ID"]
    d = dict(zip(chunk.a, chunk.k))
    umls_entrez.update(d)

entrez_umls = {v:k for k,v in umls_entrez.items()}

sem_df.SUBJECT_CUI = sem_df.SUBJECT_CUI.map(lambda x:entrez_umls[x] if x in entrez_umls else x)
sem_df.OBJECT_CUI = sem_df.OBJECT_CUI.map(lambda x:entrez_umls[x] if x in entrez_umls else x)

noncdf = sem_df[~sem_df.SUBJECT_CUI.str.startswith("C")]
print("Not start with C:",len(noncdf))

# dump everything that doesn't starts with a "C"
sem_df = sem_df[sem_df.SUBJECT_CUI.str.startswith("C")]
sem_df = sem_df[sem_df.OBJECT_CUI.str.startswith("C")]
print("Start with C:",len(sem_df))
sem_df.to_csv("edges3.csv")


idx = sem_df["PREDICATE"].str.startswith("NEG_")
sem_df['NEG'] = False
sem_df.loc[idx, 'NEG'] = True
sem_df.loc[idx, 'PREDICATE'] = sem_df[idx].PREDICATE.str.replace("NEG_", "")
print(sem_df[sem_df.NEG].head())

sem_df.to_csv("edges4.csv")

conso = pd.read_csv("MRCONSO_ENG.RRF.gz", delimiter="|", index_col=None, names = list("abcdefghijklmnopqrs"))
conso = conso[(conso['c'] == "P") & (conso['e'] == "PF")]

node_label = dict(zip(conso.a, conso.o))
print("node label length:",len(node_label))
nodes = set(sem_df.SUBJECT_CUI) | set(sem_df.OBJECT_CUI)
print("node length:",len(nodes))

nodes = pd.DataFrame({"ID":x, "LABEL": node_label.get(x)} for x in nodes)
nodes.to_csv("nodes_na.csv")
nodes = nodes.dropna()
print("After dropping no labels:",len(nodes))

nodes.to_csv("nodes1.csv")

