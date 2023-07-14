import sys, os, json
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')
# add the module's parent directory to the system path
sys.path.append(module_dir)

from loadFile import getSeqEntry
from fastaResFun import fastaCreation

if __name__ == "__main__" :

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	data_dict = config_dict["data"]
	params_dict = config_dict["parameters"]

	# if params_dict["step"] == "fasta":
	# 	data_dict = getSeqEntry(data_dict)
	
	data_dict["sptree"] = params_dict["sptree"]
	data_dict = fastaCreation(data_dict, params_dict["remote"], params_dict["APIKey"], params_dict["step"], params_dict["duplication"],sys.argv[2])

	config_dict["data"] = data_dict
	config_dict["parameters"] = params_dict
	
	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out, indent="")
		