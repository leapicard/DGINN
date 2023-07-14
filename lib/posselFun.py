import PSPFunc, GeneAnalysis, SiteAnalysis, branchAnalysis
from time import localtime, strftime
import os


def pspAnalysis(data, parms, aln, tree):
    """
    procedure which execute functions for psp step

    @param1 data: basicData object

        @return Output directory name
    """
    #logger=logging.getLogger("main.positiveSelection")
    dCtrls, lModels = PSPFunc.getParams(parms["models"], 
                        parms["paml"], 
                        parms["bppml"], 
                        parms["mixedlikelihood"], 
                        parms["busted"], 
                        parms["meme"], 
                        parms["opb"], 
                        parms["gnh"])
    timeStamp = strftime("%Y%m%d%H%M", localtime())
    
    outDir = data["o"]+"positive_selection_results_"+timeStamp+"/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    # cladoFile =  PSPFunc.supBoot(outDir, data["baseName"], tree, logger)<-- Old line, to uncomment if logger is back
    cladoFile =  PSPFunc.supBoot(outDir, data["baseName"], tree)
                    
    ### Terminal output for user
    """	
    logger.info("Output directory: {:s}".format(outDir))
    logger.info("Alignement: {:s}".format(aln))
    logger.info("Alignement is in {:s} format.".format(data["alnFormat"]))
    logger.info("Tree: {:s}".format(tree))"""

    ### Run the different analysis as determined by control file
    """
    logger.info("Starting positive selection analyses.")
    logger.info("POSITIVE SELECTION ANALYSIS: ")
    logger.info("Analysis to be run:")"""
    """
    dAnalysis = {"paml": "Site (codeml)", 
             "BUSTED":"Whole-Gene", 
             "bppml":"Site (Bio++ - Optimization)", 
             "bppmixedlikelihood":"Site (Bio++ - Results)", 
             "OPB":"Branch", 
             "GNH":"Branch-site on positively selected branches", 
             "MEME":"Branch-site"}
    for key in dCtrls.keys():
        logger.info(dAnalysis[key])"""
    
    if "BUSTED" in dCtrls:
        #GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data["baseName"], logger) <-- Old line, to uncomment if logger is back
        GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data["baseName"])			
        """try:		
            GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data.baseName, logger)
        except Exception:
            logger.info("BUSTED encountered an unexpected error, skipping.")"""		

    if "MEME" in dCtrls:
        try:
            # BranchAnalysis.memeBranchSite(aln, cladoFile, outDir, data["baseName"], logger)<-- Old line, to uncomment if logger is back
            branchAnalysis.memeBranchSite(aln, 
                              cladoFile, 
                              outDir, 
                              data["baseName"])
        except Exception:
            pass
            #logger.error("MEME encountered an unexpected error, skipping.")
            
    if "bppml" in dCtrls:
#	  try:
            if not dCtrls["bppmixedlikelihood"]:
              dCtrls["bppmixedlikelihood"]=dCtrls["bppml"]
            SiteAnalysis.bppSite(dCtrls["bppml"], 
                            dCtrls["bppmixedlikelihood"], 
                            aln, 
                            data["alnFormat"], 
                            tree, 
                            lModels, 
                            outDir, 
                            data["baseName"])
#	  except Exception:
#	    logger.error("Bio++ Site encountered an unexpected error, skipping.")
    
    lPSNodes = []
    if "OPB" in dCtrls:
#		try:
            params = branchAnalysis.bppBranch(dCtrls["OPB"], 
                              outDir, 
                              data["baseName"], 
                              aln, 
                              data["alnFormat"], 
                              tree)	
        # except Exception:
        # 	logger.error("Bio++ Branch Analysis encountered an unexpected error, skipping.")
    
    if "OPB" and "GNH" in dCtrls and len(lPSNodes) > 1:
#		try:
            branchAnalysis.bppBranchSite(dCtrls["GNH"], lPSNodes, outDir, data["baseName"], aln, data["alnFormat"], tree)
        # except Exception:
        # 	logger.error("Bio++ Pseudo Branch-Site Analysis encountered an unexpected error, skipping.")
    
    if "paml" in dCtrls:
        SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data["baseName"])
        """try:
            SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data.baseName, logger)
        except Exception:
            logger.info("PAML (codeml) Site encountered an unexpected error, skipping.")"""

    print("Finished positive selection analyses.")
    return(outDir)


def pspAnalysis2(data, parms, aln, tree):
    dCtrls, lModels = PSPFunc.getParams2(parms["models"], 
                        parms["paml"], parms["bppml"], parms["mixedlikelihood"], 
                        parms["busted"], parms["meme"], parms["opb"], parms["gnh"])
    timeStamp = strftime("%Y%m%d%H%M", localtime())
    
    outDir = data["o"]+"positive_selection_results_"+timeStamp+"/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    cladoFile =  PSPFunc.supBoot(outDir, data["baseName"], tree)
    
    if "busted" in dCtrls:
        GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, data["baseName"])

    if "meme" in dCtrls:
        try:
            branchAnalysis.memeBranchSite(aln, cladoFile, outDir, data["baseName"])
        except Exception:
            pass
            
    if "bppml" in dCtrls:
        if not dCtrls["mixedlikelihood"]:
            dCtrls["mixedlikelihood"]=dCtrls["bppml"]
        SiteAnalysis.bppSite(dCtrls["bppml"], dCtrls["mixedlikelihood"], aln, 
                        data["alnFormat"], tree, lModels, outDir, data["baseName"])

    lPSNodes = []
    if "opb" in dCtrls:
        params = branchAnalysis.bppBranch(dCtrls["opb"], outDir, data["baseName"], 
                            aln, data["alnFormat"], tree)	

    if "opb" and "gnh" in dCtrls and len(lPSNodes) > 1:
        branchAnalysis.bppBranchSite(dCtrls["gnh"], lPSNodes, outDir, data["baseName"], aln, data["alnFormat"], tree)
    
    if "paml" in dCtrls:
        SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data["baseName"])

    print("Finished positive selection analyses.")
    return(outDir)

