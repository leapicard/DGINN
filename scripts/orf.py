import sys, os, json
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')
# add the module's parent directory to the system path
sys.path.append(module_dir)

from analyseFun import getORFs
from loadFile import spTreeCheck


if __name__ == "__main__" :

	with open(sys.argv[1], 'r') as config_in:
		config_dict = json.load(config_in)

	parameters = config_dict["parameters"]
	data = config_dict["data"]
	
	# if parameters["step"] == "orf":
	# 	data = orfEntry(data)
		
	spTreeCheck(data, "orf", parameters["duplication"])

	ORFile = getORFs(data["seqFile"],
                     data["queryName"], 
                     sys.argv[2])
	data["ORFs"] = ORFile 

	config_dict["parameters"] = parameters
	config_dict["data"] = data
    
	with open(sys.argv[1],'w') as config_out:
		json.dump(config_dict, config_out, indent="")