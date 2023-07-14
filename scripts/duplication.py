import sys, os, json, shutil, logging
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')
# add the module's parent directory to the system path
sys.path.append(module_dir)
from analyseFun import checkPhyMLTree, checkPhyMLTree2
from loadFile import spTreeCheck, duplPSEntry
from treeFun import treeTreatment, treeTreatment2


if __name__ == "__main__":

    with open(sys.argv[1], 'r') as config_in:
        config_dict = json.load(config_in)

    parameters = config_dict["parameters"]
    data = config_dict["data"]

    if parameters["step"] == "option":
        if (parameters["duplication"] and not parameters["positiveSelection"] and not parameters["recombination"]) or \
            (parameters["duplication"] and parameters["positiveSelection"] and not parameters["recombination"]):
        
            data, data["dAlTree"] = duplPSEntry(data)

            dTree = data["dAlTree"].pop(data["aln"])
            spTreeCheck(data, "duplication", parameters["duplication"])
            data["dAlTree"][data["aln"]] = dTree

        dAlTree_after_check = checkPhyMLTree2(data, data["dAlTree"], parameters["nbspecies"], parameters["LBopt"])
        dAlTree_after_treatment = treeTreatment2(data, parameters, dAlTree_after_check, sys.argv[2])

        data["dAlTree"] = dAlTree_after_check
        data["dAlTree"].update(dAlTree_after_treatment)

        config_dict["parameters"] = parameters
        config_dict["data"] = data

        with open(sys.argv[1],'w') as config_out:
            json.dump(config_dict, config_out, indent="")

    else :
        if (parameters["duplication"] and not parameters["positiveSelection"] and not parameters["recombination"]) or \
            (parameters["duplication"] and parameters["positiveSelection"] and not parameters["recombination"]):

            dTree = data["dAlTree"].pop(data["aln"])
            spTreeCheck(data, "orf", parameters["duplication"])
            data["dAlTree"][data["aln"]] = dTree

            dAlTree_after_check = checkPhyMLTree(data, data["dAlTree"], parameters["nbspecies"], parameters["LBopt"])
            dAlTree_after_treatment = treeTreatment(data, dAlTree_after_check, parameters["nbspecies"], parameters["phymlOpt"], sys.argv[2])
                
            data["dAlTree"] = dAlTree_after_check
            data["dAlTree"].update(dAlTree_after_treatment)

            config_dict["parameters"] = parameters
            config_dict["data"] = data

            with open(sys.argv[1],'w') as config_out:
                json.dump(config_dict, config_out, indent="")



