

import sys, os, json
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir, '..', 'lib')
# add the module's parent directory to the system path
sys.path.append(module_dir)
from covAlnFun import covAln



if __name__ == "__main__" :
	"""
	arg1 = mafft output : fasta alignment file path
	arg2 = covAln output : fasta alignment file path
	arg3 = data dictionnary : pickle file containing dictionnary of dataObject path, before mafft update
	arg4 = data dictionnary : pickle file containing dictionnary of dataObject path, after mafft update
	"""

	with open(sys.argv[3], 'r') as config_in:
		config_dict = json.load(config_in)
	
	print(f"\nUpdating config after mafft alignment\n")

	config_dict["data"]["aln"] = sys.argv[1]
				
	print(f"\nRunning coverage check on file {sys.argv[1]}\n")

	fasCov, nbOut = covAln(aln = config_dict["data"]["aln"], 						
							      cov = 50, 
							      queryName = config_dict["data"]["queryName"],		 
							      o = sys.argv[2])

	with open(sys.argv[3],'w') as config_out:
		json.dump(config_dict, config_out, indent="")