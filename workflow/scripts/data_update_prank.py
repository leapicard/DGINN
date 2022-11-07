import pickle, sys

"""
arg1 = output from prank aligner : fasta alignment file path
arg2 = data dictionnary : pickle file containing data dictionnary path, before prank update
arg3 = data dictionnary : pickle file containing data dictionnary path, after prank update
"""


with open(sys.argv[2], 'rb') as fichier:
    data = pickle.load(fichier)

data["aln"] = sys.argv[1]

with open(sys.argv[3],'wb') as fichier_data:
    pickle.dump(data,fichier_data,pickle.HIGHEST_PROTOCOL)