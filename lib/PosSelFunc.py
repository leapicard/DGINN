import PSPFunc, GeneAnalysis, SiteAnalysis, BranchAnalysis, os
from time import localtime, strftime
import logging
import subprocess, re

def pspAnalysis(params):
    """
    procedure which execute functions for psp step

    @param1 data: basicData object

        @return Output directory name
    """
    logger=logging.getLogger("main.positiveSelection")

    tree = os.path.join(params["outdir"],params["queryName"]+"_tree.dnd")
    aln = os.path.join(params["outdir"],params["queryName"]+"_align.fasta")
    
    timeStamp = strftime("%Y%m%d%H%M", localtime())

    ### set up new directory
    
    outDir = os.path.join(params["outdir"],params["queryName"] + "_positive_selection")
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    params["outDir"]=outDir
    
    cladoFile =  PSPFunc.supBoot(params)
                    
    ### Terminal output for user
    logger.info("Output directory: {:s}".format(outDir))
    logger.info("Alignement: {:s}".format(aln))
    logger.info("Tree: {:s}".format(tree))

    ### Run the different analysis as determined by control file
    logger.info("Starting positive selection analyses.")
    logger.info("POSITIVE SELECTION ANALYSIS: ")
    logger.info("Analysis to be run:")

    ###########################################################
    #### HYPHY

    if params["busted"]:
      try:		
        GeneAnalysis.hyphyBusted(aln, cladoFile, outDir, logger)
      except Exception:
        logger.info("BUSTED encountered an unexpected error, skipping.")

    if params["meme"]:
       try:
            BranchAnalysis.memeBranchSite(aln, cladoFile, outDir, logger)
       except Exception:
            logger.error("MEME encountered an unexpected error, skipping.")


    ###########################################################
    #### BPP

    lModels = list(map(str.strip,re.compile(r"[,;]").split(params["models"])))

    if params["bppml"] or params["opb"] or params["gnh"]:
      outBPP = outDir+"/bpp_site"
      if not os.path.exists(outBPP):
        subprocess.Popen("mkdir "+outBPP, shell =  True).wait()

    if params["bppml"]:
      params["bppml"]=outBPP+"/base.bpp"
      PSPFunc.pspFileCreation(params["bppml"],"bppml")
      
      if not params["mixedlikelihood"]:
        params["mixedlikelihood"]=params["bppml"]  # not used
      else:
        params["mixedlikelihood"]=outBPP+"/base_mll.bpp"
        PSPFunc.pspFileCreation(params["mixedlikelihood"],"bppmixedlikelihood")
                              
#      try:
      SiteAnalysis.bppSite(aln, 
                             tree, 
                             outDir, 
                             params["bppml"], 
                             params["mixedlikelihood"], 
                             lModels, 
                             logger)
      # except Exception:
      #   logger.error("Bio++ Site encountered an unexpected error, skipping.")
    
    lPSNodes = []
    if params["opb"]:
        params["opb"]=outBPP+"/base_opb.bpp"
        PSPFunc.pspFileCreation(params["opb"],"opb")
        try:
          params = BranchAnalysis.bppBranch(aln, 
                                            tree, 
                                            outDir, 
                                            params["opb"], 
                                            logger)	
        except Exception:
          logger.error("Bio++ Branch Analysis encountered an unexpected error, skipping.")
    
    if params["opb"] and params["gnh"] and len(lPSNodes) > 1:
      params["gnh"]=outBPP+"/base_gnh.bpp"
      PSPFunc.pspFileCreation(params["gnh"],"gnh")
      try:
            BranchAnalysis.bppBranchSite(aln, tree, outDir, params["gnh"], lPSNodes, logger)
      except Exception:
       	logger.error("Bio++ Pseudo Branch-Site Analysis encountered an unexpected error, skipping.")


    ###########################################################
    #### PAML
        
    if params["paml"]:
        SiteAnalysis.pamlSite(aln, 
                              tree, 
                              outDir, 
                              params["paml"], 
                              lModels,
                              logger)

        """try:
            SiteAnalysis.pamlSite(aln, tree, lModels, dCtrls["paml"], outDir, data.baseName, logger)
        except Exception:
            logger.info("PAML (codeml) Site encountered an unexpected error, skipping.")"""



    logger.info("Finished positive selection analyses.")
    return(outDir)
