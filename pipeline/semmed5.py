import sys
import os
import pickle

import pandas as pd
from numpy import nan
import seaborn as sns
import shelve
import re
from collections import defaultdict, Counter
from tqdm import tqdm
from tqdm.notebook import tqdm as tqdm_notebook
from itertools import chain
from more_itertools import chunked
from collections import Counter
from pprint import pprint
import requests
from pyquery import PyQuery as pq
from wikidataintegrator import wdi_helpers, wdi_core, wdi_login


uri_to_curie = lambda s: s.split("/")[-1].replace("_", ":")
nodes = pd.read_csv("nodes_biolink.csv",index_col=0)

names = "CUI,LAT,TS,LUI,STT,SUI,ISPREF,AUI,SAUI,SCUI,SDUI,SAB,TTY,CODE,STR,SRL,SUPPRESS,CVF,X".split(",")
umls = pd.read_csv("MRCONSO_ENG.RRF.gz", delimiter="|", names=names, index_col=None)
# only get CUIs in our list of nodes
umls = umls[umls.CUI.isin(nodes.index)]

umls['xref'] = umls.SAB + ":" + umls.CODE.map(str)
# easy fix to HGNC prefix duplication between SAB and CODE
umls.xref = umls.xref.str.replace("HGNC:HGNC:", "HGNC:")
# fix this MSH MESH nonsense
umls.xref = umls.xref.str.replace("MSH:", "MESH:")
# NCI_FDA is UNII
umls.xref = umls.xref.str.replace("NCI_FDA:", "UNII:")

XREF = dict(umls.groupby("CUI")['xref'].apply(set))
XREF = defaultdict(set, XREF)

# what xrefs are on chemicals?
chem_umls = nodes[nodes.blm_category == "chemicalsubstance"].index ###lower case
xref_chem = {k:v for k,v in XREF.items() if k in chem_umls}
#print(len(chem_umls))
c = Counter(list(chain(*[list(map(lambda x:x.split(":",1)[0], y)) for y in xref_chem.values()])))
#pprint(c.most_common(25))


pd.set_option("display.width", 120)

URL = "http://id.nlm.nih.gov/mesh/sparql"
PREFIX = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
"""

def sparql_query(query):
    params = {'query': PREFIX + query, 'format': 'JSON', 'limit': 1000, 'offset': 0}
    r = requests.get(URL, params=params)
    res = [{k: v['value'] for k, v in x.items()} for x in r.json()['results']['bindings']]
    t = tqdm()
    while True:
        t.update(1)
        params['offset'] += 1000
        r = requests.get(URL, params=params).json()['results']['bindings']
        if not r:
            break
        res.extend([{k: v['value'] for k, v in x.items()} for x in r])
    df = pd.DataFrame(res)
    return df


query = """
SELECT distinct ?mesh ?meshLabel ?r ?rr
FROM <http://id.nlm.nih.gov/mesh> WHERE {
  ?mesh meshv:active 1 .
  ?mesh meshv:preferredMappedTo ?p .
  ?p meshv:treeNumber ?treeNum .
  FILTER(STRSTARTS(STR(?treeNum), "http://id.nlm.nih.gov/mesh/D")) .
  ?mesh rdfs:label ?meshLabel .
  ?mesh meshv:preferredConcept [meshv:registryNumber ?r] .
  #OPTIONAL {?mesh meshv:preferredConcept [meshv:relatedRegistryNumber ?rr]}
}
"""
df = sparql_query(query)

df.r = df.r.replace("0", nan)
df.dropna(subset=["r"], inplace=True)
df = df[~df.r.str.startswith("EC ")]
df.mesh = df.mesh.str.replace("http://id.nlm.nih.gov/mesh/", "")
df.set_index("mesh", inplace=True)

df.to_csv("mesh_xrefs.csv")

mesh_xrefs = pd.read_csv("mesh_xrefs.csv",index_col=0)

mesh_xrefs.r = mesh_xrefs.r.apply(lambda x: "CAS:" + x if "-" in x else "UNII:" + x)
mesh_xrefs = mesh_xrefs.groupby("mesh").r.apply(set).to_dict()
mesh_xrefs = {"MESH:"+k:v for k,v in mesh_xrefs.items()}

for k,v in xref_chem.items():
    for vv in list(v):
        if vv in mesh_xrefs:
            v.update(mesh_xrefs[vv])

unii_df = pd.read_csv("UNII_Records_12Jun2021.txt", dtype=str, sep='\t', low_memory=False)
unii_df.dropna(subset=['INCHIKEY'], inplace=True)

n=0
for k,v in tqdm_notebook(xref_chem.items()):
    for vv in list(v):
        if vv.startswith("UNII:"):
            xref = vv.replace("UNII:", "")
            s = unii_df.query("UNII == @xref").INCHIKEY
            if not s.empty:
                n+=1
                v.add("INCHIKEY:" + list(s)[0])

xref_inchi = {k:v for k,v in xref_chem.items() if any(vv.startswith("INCHIKEY:") for vv in v)}
xref_inchi = {k:[vv for vv in v if vv.startswith("INCHIKEY:")][0].replace("INCHIKEY:", "") for k,v in xref_inchi.items()}
print(len(xref_inchi))


url = "https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_structures__standard_inchi_key__in={}&format=json&limit=100"
for chunk in tqdm(chunked(xref_inchi.items(), 100), total=len(xref_inchi)/100):
    chunk = dict(chunk)
    chunk = {v:k for k,v in chunk.items()}
    inchis = ",".join(chunk)
    mols = requests.get(url.format(inchis)).json()['molecules']
    for m in mols:
        chembl = m['molecule_chembl_id']
        inchi = m['molecule_structures']['standard_inchi_key']
        XREF[chunk[inchi]].add("CHEMBL:" + chembl)

print(len({k:v for k,v in XREF.items() if any(vv.startswith("CHEMBL:") for vv in v)}))

with open("xrefs.shelve", 'wb') as f:
    pickle.dump(XREF, f)

df = pd.read_csv("uberon.csv")
df = df[df.xref.str.startswith("UMLS:")]
df.xref = df.xref.str.replace("UMLS:", "")
df.item = df.item.apply(uri_to_curie)

s = df.groupby("xref")['item'].apply(set)
for umls, x in dict(s).items():
    XREF[umls].update(x)

print(XREF['C1272528'])
print("---------")

df = pd.read_csv("doid.csv")
df.dropna(inplace=True)
df = df[df.xref.str.startswith("UMLS_CUI:")]
df.xref = df.xref.str.replace("UMLS_CUI:", "")
df.item = df.item.apply(uri_to_curie)
print(df.head())


s = df.groupby("xref")['item'].apply(set)
for umls, x in dict(s).items():
    XREF[umls].update(x)

print(XREF['C0263518'])
print("---------")

names = list("abcdefghijklmn")
iter_csv = pd.read_csv("MRSAT.RRF.gz", delimiter="|", names=names, index_col=None, chunksize=1000000)
chunks = []
umls_uniprot = dict()
for chunk in tqdm(iter_csv, total=67668372/1000000):
    chunk.fillna(method='ffill', inplace=True)
    chunk = chunk[chunk.i == "SWISS_PROT"]
    d = dict(zip(chunk.a, chunk.k))
    umls_uniprot.update(d)


print(len(umls_uniprot))

for umls, uniprot in umls_uniprot.items():
    XREF[umls].add("UNIPROT:" + uniprot)

print(XREF['C0215993'])
print("---------")
with open("xrefs.shelve", 'wb') as f:
    pickle.dump(XREF, f)

nodes['xrefs'] = nodes.index.map(lambda x: ";".join(XREF.get(x,list())))

print(nodes.head())
nodes.to_csv("nodes_xref.csv")
