import PSPFunc, ExtractFunc, GeneAnalysis, SiteAnalysis, BranchAnalysis, os
from time import localtime, strftime
import logging

def pspAnalysis(data, parms, aln, tree):
	"""
	procedure which execute functions for psp step

	@param1 data: basicData object
	"""
	logger=logging.getLogger("main.positiveSelection")
	dCtrls, lModels = PSPFunc.getParams(parms["models"], parms["paml"], parms["bppml"], parms["mixedlikelihood"], parms["busted"], parms["meme"], parms["opb"], parms["gnh"])
	timeStamp = strftime("%Y%m%d%H%M", localtime())
	outDir = data.o+"positive_selection_results_"+timeStamp+"/"
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	cladoFile =  PSPFunc.supBoot(outDir, data.baseName, tree, logger)
					
	### Terminal output for user
	logger.info("Location of entry files: {:s}".format(outDir))
	logger.info("Alignement: {:s}".format(aln))
	logger.info("Alignement is in {:s} format.".format(data.alnFormat))
	logger.info("Tree: {:s}".format(tree))

	nodes = PSPFunc.nbNode(tree, logger)
					
	### Run the different analysis as determined by control file
	logger.info("POSITIVE SELECTION ANALYSIS: ")
	logger.info("Analysis to be run:")

	dAnalysis = {"paml": "Site (codeml)", "BUSTED":"Whole-Gene", "bppml":"Site (Optimization)", "bppmixedlikelihood":"Site (Results)", "OPB":"Branch", "GNH":"Branch-site on positively selected branches", "MEME":"Branch-site"}
	for key in dCtrls.keys():
		logger.info(dAnalysis[key])
	
	if "BUSTED" in dCtrls:
		GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data.baseName, logger)			
		"""try:		
			GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data.baseName, logger)
		except Exception:
			logger.info("BUSTED uncountered an unexpected error, skipping.")"""
	
	if "paml" in dCtrls and dCtrls["paml"] not in ["False", False] and len(lModels) > 1:
		SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data.baseName, logger)
		"""try:
			SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data.baseName, logger)
		except Exception:
			logger.info("PAML (codeml) Site uncountered an unexpected error, skipping.")"""
			
			
	if "bppml" and "bppmixedlikelihood" in dCtrls and len(lModels) > 1:
		try:
			SiteAnalysis.bppSite(dCtrls["bppml"], dCtrls["bppmixedlikelihood"], aln, data.alnFormat, tree, lModels, outDir, data.baseName, logger)
		except Exception:
			logger.error("Bio++ Site uncountered an unexpected error, skipping.")
	elif "bppml" or "bppmixedlikelihood" not in dCtrls:
		logger.error("Part of parameters for Bio++ site analysis are completed but not all.")
		logger.error("Analysis ignored (if unexpected, check paths to Bio++/bpp parameter files).")
	elif "bppml" and "bppmixedlikelihood" not in dCtrls:
		next
	
	lPSNodes = []
	if "OPB" in dCtrls:
		try:
			params = BranchAnalysis.bppBranch(dCtrls["OPB"], outDir, data.baseName, aln, data.alnFormat, tree, logger)	
		except Exception:
			logger.error("Bio++ Branch Analysis uncountered an unexpected error, skipping.")
		try:
			lPSNodes = BranchAnalysis.parseNodes(params, logger)	
		except Exception:
			logger.error("Could not find info about positively selected branches in Bio++ results. If unexpected, check Bio++ parameter files for OPB.")
	
	if "OPB" and "GNH" in dCtrls and len(lPSNodes) > 1:
		try:
			BranchAnalysis.bppBranchSite(dCtrls["GNH"], lPSNodes, outDir, data.baseName, aln, data.alnFormat, tree, logger)
		except Exception:
			logger.error("Bio++ Pseudo Branch-Site Analysis uncountered an unexpected error, skipping.")
	
	if "MEME" in dCtrls:
		try:
			BranchAnalysis.memeBranchSite(aln, cladoFile, outDir, data.baseName, logger)
		except Exception:
			logger.error("MEME uncountered an unexpected error, skipping.")
					
	logger.info("End analysis")
	return(outDir)
