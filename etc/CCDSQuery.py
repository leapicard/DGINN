import CCDSQueryFunc as CQF
import argparse

if __name__ == "__main__":

	parser = argparse.ArgumentParser(prog=__file__, description='''This program get sequences' genes from HGNC Biomart.''')
	required = parser.add_argument_group('Mandatory input infos for running')
	required.add_argument('-in', '--inFile', metavar="<filename>", type = CQF.check, required=True, dest = 'inFile', help = 'Table of HGNC approved symbols (one per line) and corresponding CCDS accessions for the genes of interest, obtained from HGNC Biomart')
	#required.add_argument('-sp', '--species', metavar="<filename>", type = str, required=True, dest = 'species', help = 'Species for the request in the HGNC Biomart')
	#required.add_argument('-spn', '--spName', metavar="<filename>", type = str, required=True, dest = 'spName', help = "Species's Name of the genes")

	parameters = parser.parse_args()
	parameters.species = "human"
	parameters.spName = "HomSap"
	CQF.getCCDS(parameters.inFile, parameters.species, parameters.spName)
