import sys, os, json
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')
# add the module's parent directory to the system path
sys.path.append(module_dir)

from loadFile import pspEntry
from posselFun import pspAnalysis, pspAnalysis2


if __name__ == "__main__" :	
    with open(sys.argv[1], 'r') as config_in:
        config_dict = json.load(config_in)

    parameters = config_dict["parameters"]
    data = config_dict["data"]
    
    data, dAlTree = pspEntry(data, parameters)

    listArgsPosSel =  []
    filename = data["o"]+"possel_"+data["baseName"]
    fAT = open(filename+"_files_list.txt", "w")

    for aln, tree in data["dAlTree"].items():
        if len(data["dAlTree"]): 
            fAT.write(aln+"\n"+tree+"\n")
            outDir = pspAnalysis2(data, parameters, aln, tree)

