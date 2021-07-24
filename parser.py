from collections import defaultdict
import csv
import os
import sys
import pandas as pd
import yaml
import json

maxInt = sys.maxsize

#DOC_TYPE = "DiseaseOrPhenotypicFeature"
#SEMMED_TYPE = "diseaseorphenotypicfeature"

#DOC_TYPE = "ChemicalSubstance"
#SEMMED_TYPE = "chemicalsubstance"

#DOC_TYPE = "PhenotypicFeature"
#SEMMED_TYPE = "phenotypicfeature"

DOC_TYPE = "Gene"
SEMMED_TYPE = "gene"

#DOC_TYPE = "AnatomicalEntity"
#SEMMED_TYPE = "anatomicalentity"

#DOC_TYPE = "BiologicalProcess"
#SEMMED_TYPE = "biologicalprocessoractivity"

JSON_FILE_NAME=SEMMED_TYPE+".json"
print(JSON_FILE_NAME)

SEMMED_SEMANTIC_TYPE_MAPPING = {
    "chemicalsubstance": "ChemicalSubstance",
    "activityandbehavior": None,
    "anatomicalentity": "AnatomicalEntity",
    "biologicalentity": None,
    "biologicalprocessoractivity": "BiologicalProcess",
    "cell": "Cell",
    "cellularcomponent": "CellularComponent",
    "diseaseorphenotypicfeature": "DiseaseOrPhenotypicFeature",
    "gene": "Gene",
    "genomicentity": None,
    "grossanatomicalstructure": None,
    "phenotypicfeature": "PhenotypicFeature",
    "protein": "Gene"
}

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

def parse_biolink_yaml(yaml_file):
   print("Starting parsing YAML file:",yaml_file)
   with open(yaml_file) as file:
    documents = yaml.full_load(file)
    inverse_pair={}
    SEMMED_PRED_MAPPING={}
    semmed_list=[]
    for item, doc in documents.items():
        if (item=='slots'):

         for doc_key in doc:
           biolink_key=doc_key.replace(" ","_")
           if ('inverse' in doc[doc_key]):
               inverse=doc[doc_key]['inverse'].strip()
               inverse_pair[doc_key]=inverse

           for item_key in doc[doc_key].keys():
            if ('mappings' in item_key): ### no matter what kind of mappings. E.g., exact, narrow...
             mappings=doc[doc_key][item_key]
             if (mappings!=None):
                for mapping in mappings:
                  if ('SEMMEDDB' in mapping):
                    semmed_list.append(["biolink:"+biolink_key,item_key,mapping.replace("SEMMEDDB:","")])
    
    for rec in semmed_list:
      check=rec[0].replace("biolink:"," ").replace("_"," ").strip()
      semmed_type=rec[2]
      found=0
      for key,value in inverse_pair.items():
       if (value.strip()==check):
          found=1
          SEMMED_PRED_MAPPING[semmed_type]={'self':check.replace(' ','_'),'reverse':key.replace(' ','_')}
       elif (key.strip()==check): ###key check for part of only 07/2021
          found=1
          SEMMED_PRED_MAPPING[semmed_type]={'self':check.replace(' ','_'),'reverse':value.replace(' ','_')}
      if (found==0):
         SEMMED_PRED_MAPPING[semmed_type]={'self':check.replace(' ','_'),'reverse':check.replace(' ','_')}
    print("Finished!")
    return SEMMED_PRED_MAPPING

def load_data(data_folder):
    def construct_rec(sub_umls, obj_umls, line, reverse=False):

        if sub_umls not in gene_related:

            gene_related[sub_umls] = {
                '_id': sub_umls[5:],
                'umls': sub_umls[5:],
                'name': id_type_mapping[sub_umls]['name'],
                '@type': DOC_TYPE
            }
        semmed_pred=line[0].lower()
        if not reverse:
             pred = SEMMED_PRED_MAPPING[semmed_pred]['self']
        else:
             pred = SEMMED_PRED_MAPPING[semmed_pred]['reverse']

        semantic_type = SEMMED_SEMANTIC_TYPE_MAPPING[id_type_mapping[obj_umls]['type']]
        if semantic_type:
            if pred not in gene_related[sub_umls]:
                gene_related[sub_umls][pred] = {}
            pmid=str(line[1])
            assoc = sub_umls + pred + pmid + obj_umls
           
            if assoc not in unique_assocs:
                unique_assocs.add(assoc)
                if obj_umls not in gene_related[sub_umls][pred]:
                    gene_related[sub_umls][pred][obj_umls] = {
                        "pmid": set(),
                        'umls': obj_umls[5:],
                        'name': id_type_mapping[obj_umls]['name'],
                        '@type': semantic_type
                    }
                gene_related[sub_umls][pred][obj_umls]["pmid"] = gene_related[sub_umls
                                                                              ][pred][obj_umls]["pmid"] | set(line[1].split(';'))


    nodes_path = os.path.join(data_folder, "nodes_neo4j.csv")
    edges_path = os.path.join(data_folder, "edges_neo4j.csv")
    SEMMED_PRED_MAPPING=parse_biolink_yaml(os.path.join(data_folder, "biolink-model.yaml"))
   
    group_by_semmantic_dict = defaultdict(list)
    id_type_mapping = {}
    unique_assocs = set()
    with open(nodes_path) as f:
        csv_reader = csv.reader(f, delimiter=',')
        next(csv_reader)
        count=0
        for _item in csv_reader:
            sem_type=_item[-2]
            umls_id=_item[-1]
            name=_item[1] 
            if (sem_type=='protein'):
               sem_type='gene'
 
            group_by_semmantic_dict[sem_type].append(umls_id) ###put UMLS ids in the same SEM type
            id_type_mapping[umls_id] = {'type':sem_type, 'name': name}

 
    gene_related = {}
    with open(edges_path) as f:
        csv_reader = csv.reader(f, delimiter=',')
        next(csv_reader)
        count=0
        
        for _item in csv_reader:
             subject_id=_item[4]
             object_id=_item[5]
             count+=1
             print("Current:",count)
             if subject_id in group_by_semmantic_dict[SEMMED_TYPE]: ###check if UMLS id in one SEM type
                construct_rec(subject_id, object_id, _item)
                           
             if object_id in group_by_semmantic_dict[SEMMED_TYPE]:
                construct_rec(object_id, subject_id, _item, reverse=True)

    with open(JSON_FILE_NAME, 'a') as outfile:
     doc=[]         
     for v in gene_related.values():
        for m, n in v.items():
          if m not in ["_id", "umls", "name", "@type"]:
                tmp = []
                for item in n.values():
                    item['pmid'] = list(item['pmid'])
                    tmp.append(item)
                v[m] = tmp
        #print(v) #yield v
        doc.append(v)
     json.dump(doc,outfile,indent=4)

load_data("./")
