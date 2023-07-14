
import sys, os, json
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')

# add the module's parent directory to the system path
sys.path.append(module_dir)

# import the module using its relative path
from loadFile import baseNameInit
from blastFun import blastTreatment, parseBlast, setGenAttr
from initConfig import paramDef, initLogger



if __name__ == "__main__" :
    with open(sys.argv[1], 'r') as json_in :
        json_dict = json.load(json_in)

    parameters = json_dict["parameters"]
    data = json_dict["data"]

    data["baseName"] = baseNameInit(data["baseName"], data["queryFile"], data["aln"])
        
    data["blastRes"] = blastTreatment(data["queryFile"], 
                            data["o"],
                            data["baseName"],
                            data["db"], 
                            parameters["evalue"], 
                            parameters["percID"], 
                            parameters["mincov"], 
                            parameters["APIKey"], 
                            parameters["remote"], 
                            parameters["entryQuery"])

    data["lBlastRes"] = parseBlast(data["blastRes"])
    data,params = setGenAttr(data,parameters)
    json_dict["parameters"] = params
    json_dict["data"] = data

    with open(sys.argv[1], 'w') as json_out :
        json.dump(json_dict, json_out, indent="")