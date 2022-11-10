import json, sys

"""
arg1 = output from prank aligner : fasta alignment file path
arg2 = data dictionnary : pickle file containing data dictionnary path, before prank update
arg3 = data dictionnary : pickle file containing data dictionnary path, after prank update
"""

print(f"\nUpdating the config file after alignment using {sys.argv[3]}\n")

with open(sys.argv[2], 'r') as config_in:
    config_dict = json.loads(config_in.read())

config_dict["data"]["aln"] = sys.argv[1]                #!#

config_dict_updated = json.dumps(config_dict)

with open(sys.argv[2],'w') as config_out:
    config_out.write(config_dict_updated)
