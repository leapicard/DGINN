
import sys, os, json
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')
# add the module's parent directory to the system path
sys.path.append(module_dir)

from clusterIsoFun import isoformAln

if __name__ == "__main__" :

        with open(sys.argv[1], 'r') as config_in :
                config_dict = json.load(config_in)
                
        aln = sys.argv[3]
        config_dict["data"]["aln"] = isoformAln(aln, sys.argv[2])

        with open(sys.argv[1], 'w') as config_out :
                json.dump(config_dict, config_out, indent="")
