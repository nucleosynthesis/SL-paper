import yaml

## Simple functions to return list from Yaml table for Simplified Likelihoods

def yaml_table_to_pylist(yaml_file, table_index):
  ## The first argument is the path to the .yaml table file and the table index 
  ## is the column (0 is the first dependant variable column) which is to be read 
  with open(yaml_file) as fp: config = yaml.load(fp)
  values = config["dependent_variables"][table_index]["values"]
  data   = [values[k]["value"] for k in range(len(values))]
  return data

def yaml_multi_table_to_pylist(yaml_file):
  ## The first argument is the path to the .yaml table file 
  with open(yaml_file) as fp: config = yaml.load(fp)
  ret_list = []
  indep_vals = config["dependent_variables"][0]["values"]
  ret_list.extend([indep_vals[k]["value"] for k in range(len(indep_vals))])
  for r in range(1,len(indep_vals)):
   values = config["dependent_variables"][r]["values"]
   data   = [values[k]["value"] for k in range(len(values))]
   ret_list.extend(data)
  return ret_list
