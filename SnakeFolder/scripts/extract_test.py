from blast_test import parseBlast

"""This file pools the necessary functions to treat the input file of genes and their CCDS accessions."""

def makeAccnsFile(BlastRes, output_file):
	listAcc = parseBlast(BlastRes)
	geneAllAccns = "\n".join(set(listAcc))
		
	with open(output_file, "w") as out:
		out.write(geneAllAccns)
		out.close()	

if __name__ == "__main__" :
    queryFile = snakemake.input[0]
    o = snakemake.output[0]
    makeAccnsFile(queryFile,o)