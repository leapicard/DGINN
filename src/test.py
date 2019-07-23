#import PSPFunc, PosSelFunc, GeneAnalysis, SiteAnalysis, BranchAnalysis, ExtractFunc, os, re
#import PSPFunc
#from importlib import reload
#reload(PSPFunc)
import TreeFunc

class basicData:
	def __init__(self, outDir, baseName, alnFormat):
		self.o = outDir
		self.baseName = baseName
		self.alnFormat = alnFormat

class logger:
	def __init__(self):
		self = self
	
	def info(self, var):
		print(var)
	
	def debug(self, var):
		print(var)

def add(a):
	s = 0
	for i in a:
		s += i
	return s


if __name__ == "__main__":
	"""kh = "/home/lpicard/22genes/22geneslist_complete_GBP5_CCDS_results_201904081404/GBP5_primates_filtered_longestORFs_prank_mincov_prank.gard"
	aln = "/home/lpicard/22genes/22geneslist_complete_GBP5_CCDS_results_201904081404/GBP5_primates_filtered_longestORFs_prank_mincov_prank.best.fas"
	pvalue = 0.05
	o = "/home/lpicard/22genes/22geneslist_complete_GBP5_CCDS_results_201904081404/"
	logger = logger()
	
	out = AnalysisFunc.procGARD(kh, aln)
	print(out)
	
	out2 = AnalysisFunc.parseGard(out, aln, pvalue, o, logger)
	print(out2)"""
	
	ORF = "/home/lpicard/22genes/201903_cov/22geneslist_complete_GADD45A_CCDS_results_201903241129/GADD45A_primates_filtered_longestORFs_prank_mincov.fasta"
	recTree = "/home/lpicard/22genes/201903_cov/22geneslist_complete_GADD45A_CCDS_results_201903241129/GADD45A_primates_filtered_longestORFs_prank_mincov_reconciliation.nhx"
	nbSp = 8
	o = "/home/lpicard/"
	logger = logger()
	
	TreeFunc.treeParsing(ORF, recTree, nbSp, o, logger)
	#TreeFunc.treeParsing2(ORF, recTree, nbSp, o, logger)
