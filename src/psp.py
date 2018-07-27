import geneAnalysis as geneA
import siteAnalysis, branchAnalysis, PSPFunc, AccessionFunc

def pspAnalysis(data, parms, dAlTree, debug, logger):
	"""
	procedure which execute functions for psp step

	@param1 data: basicData object
	@param2 logger: Logging object
	"""

	dCtrls, lModels = PSPFunc.getParams(parms["models"], parms["bppml"], parms["mixedlikelihood"], parms["opbFile"], parms["gnhFile"], parms["busted"], parms["meme"], parms["opb"], parms["gnh"])

	for aln in dAlTree:
		outDir = AccessionFunc.createGeneDir(data.o, aln.split("/")[-1].split(".")[0])
		cladoFile =  PSPFunc.supBoot(outDir, data.baseName, dAlTree[aln], logger)
						
		### Terminal output for user
		logger.info("Location of entry files: %s" %outDir)
		logger.info("Alignement: %s" %aln)
		logger.info("Alignement is in %s format." %data.alnFormat)
		logger.info("Tree: %s" %dAlTree[aln])

		nodes = PSPFunc.nbNode(dAlTree[aln], outDir, logger)
						
		### Run the different analysis as determined by control file
		logger.info("POSITIVE SELECTION ANALYSIS: ")
		logger.info("Analysis to be run:")

		dAnalysis = {"BUSTED":"Whole-Gene", "bppml":"Site (Optimization)", "bppmixedlikelihood":"Site (Results)", "OPB":"Branch", "GNH":"Branch-site on positively selected branches", "MEME":"Branch-site"}
		for key in dCtrls.keys():
			logger.info(dAnalysis[key])
						
		#try:		
		geneA.hyphyBusted(dCtrls, aln, cladoFile, outDir, data.baseName, debug, logger)
		#except Exception:
		#	logger.info("BUSTED uncountered an unexpected error, skipping.")
		#try:
		siteAnalysis.bppSite(dCtrls, aln, data.alnFormat, dAlTree[aln], nodes, lModels, outDir, data.baseName, debug, logger)
		#except Exception:
		#	logger.info("Bio++ Site uncountered an unexpected error, skipping.")
		#try:
		lPSNodes = branchAnalysis.bppBranch(dCtrls, outDir, data.baseName, aln, data.alnFormat, dAlTree[aln], logger)		
		branchAnalysis.bppBranchSite(dCtrls, lPSNodes, outDir, data.baseName, logger, aln, data.alnFormat, dAlTree[aln])
		#except Exception:
			#logger.info("Bio++ Branch Analysis uncountered an unexpected error, skipping.")
		#try:
		branchAnalysis.memeBranchSite(outDir, dCtrls, data.baseName, aln, cladoFile, logger)
		#except Exception:
			#logger.info("MEME uncountered an unexpected error, skipping.")
						
	logger.info("End analysis")