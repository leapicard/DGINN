import tree_test
import os, sys, logging
from Bio import SeqIO


def filterData(sptree, filePath, o):
	"""
	Function which execute functions to delete genes which aren't in the species tree.

	@param1 sptree: Path of a newick file
	@param2 filePath: Path of a Fasta file
	@param3 o: Path of a directory
	@return path: Path of a file
	"""
	
	corsg = tree_test.assocFile(sptree, filePath, o)

	path = tree_test.supData(filePath, corsg, o)

	return path, corsg

