import yaml
import json
import pandas as pd


def parse_biolink_yaml(result,output_file):
   example_file=output_file.replace("_new","")
   print("Read example file:",example_file)  ### Read existing file and replace variables

   with open(example_file) as file:
    documents = yaml.full_load(file)
    for item, doc in documents.items():
      if (item=='paths'):
        for key in doc.keys():
            if (key=="/query"):
               for key1 in doc[key].keys():
                   if (key1=='post'):
                      doc[key][key1]['x-bte-kgs-operations'].clear()
                      for item in result:
                          ref_key="$ref"
                          doc[key][key1]['x-bte-kgs-operations'].append({ref_key: '#/components/x-bte-kgs-operations/'+item})
      elif (item=='components'):
           for key in doc.keys():
               if (key=='x-bte-kgs-operations'):
                   doc[key].clear()
                   for item in result:
                       item_list=item.split("-") 
                       value=[{'inputSeparator': ',', 'inputs': [{'id': 'UMLS', 'semantic':item_list[0]}], 'method': 'post', 'source': 'SEMMED', 'outputs': [{'id': 'UMLS', 'semantic': item_list[2]}], 'parameters': {'fields': item_list[1] }, 'path': '/query', 'requestBody': {'body': {'q': '{inputs[0]}', 'scopes': 'umls'}}, 'supportBatch': True, 'response_mapping': {'$ref': '#/components/x-bte-response-mapping/'+item}, 'predicate': item_list[1]}]  
                       doc[key][item]=value
               elif(key=='x-bte-response-mapping'):
                   doc[key].clear()
                   for item in result:
                       predicate=item.split("-")[1]
                       pred_umls=predicate+".umls"
                       pred_pmid=predicate+".pmid"
                       value={'umls': pred_umls, 'pmid': pred_pmid}
                       doc[key][item]=value

   with open(output_file, 'w') as file: ### produce new file
       yaml.dump(documents, file, sort_keys=False)
   print("New produced: "+output_file+" is ready for review!")
   print("------")

def parse_json(json_source):
    f = open("../"+json_source,)
    data = json.load(f)

    all_comb=[]
    for i in data:
      source=None
      target=None
      predicate=None
      
      for key in i.keys():
        if (key !='_id' and key!='umls' and key!='name'):
           if (type(i[key]) is not list):
              source=i[key]
           else:
              predicate=key
              target=i[key][-1]['@type']
              edge=source+"-"+predicate+"-"+target
              if (edge not in all_comb):
                 all_comb.append(edge)
    return all_comb

def main():
    #json_source="semmed_gene.json"
    json_source_list=["semmed_anatomy.json","semmed_biological_process.json","semmed_chemical.json","semmed_disease.json","semmed_phenotype.json"]
    for source_file in json_source_list:
      output_file=source_file.replace(".json","_new.yaml")
      result=parse_json(source_file)
      parse_biolink_yaml(result,output_file) 

if __name__ == "__main__":
    main()
