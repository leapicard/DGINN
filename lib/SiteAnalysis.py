import sys, os, logging, subprocess
import PSPFunc
from ete3 import EvolTree

def bppSite(bppFile, bppMixed, alnFile, alnFormat, treeFile, lModels, outDir, baseName, logger):	
	# outDir=os.getcwd()+"/"  # used to debug
	logger.info(os.getcwd())
	### SITE ANALYSIS: BIO++
	dallMod={"M0":0,"M1":1,"M2":2,"M7":3,"M8a":4,"M8":5,"DFP07_0":1,"DFP07":2}
	lModels.sort(key = lambda x : dallMod[x])

	logger.info("Bio++ Site Analysis")
	logger.info("Models to be run: {:s}".format(", ".join(model for model in lModels)))
	logger.info("Bppml parameter file: {:s}".format(bppFile))
	
	## Bppml
	""" 
	Optimize tree and model using bppml
	Variables to include are
		INPUTFILE - alignement file
		FORMAT - format of the aln file (here, phyx)
		TREEFILE - tree file for the analyzed aln
		MODEL - choose which model you want run on the data YNGP_M0 through 8, same models as PAML, and DFP07 models
		IGNORE - parameters to ignore for optimization, for example if one is fixated (ex: omegas in M8a)
		OUTTREE - name of the optimized output tree
		OUTPARAMS - name of the output file summarizing parameters
		BACKUP - name of log file
	"""

	# Bppml output file names - dictionaries that associate model number with output file name for the model
	outSite = outDir+"bpp_site/"
	if not os.path.exists(outSite):
		subprocess.Popen("mkdir "+outSite, shell =  True).wait()

	outFileName = outSite+baseName
	dModelTrees = {model:outFileName+"_"+model+".dnd" for model in lModels}
	dModelParams = {model:outFileName+"_"+model+".params" for model in lModels}
	dModelLog = {model:outFileName+"_optimization_"+model for model in lModels}
	dModelSyntax = {model:["YNGP_"+model,"frequencies=F3X4(initFreqs=observed)"] for model in lModels if model[0]=="M"}		# dictionary model number - [MODEL name, MODEL arguments for bppml]
	dModelSyntax.update({model:[model[:5],"protmodel=JTT92", "frequencies=F3X4(initFreqs=observed)"] for model in lModels if model[:5]=="DFP07"})
	# take into account the specificities of each model (number of classes n for example)
	for model in lModels:
          if model in ["M7","M8","M8a"]:
            dModelSyntax[model].append("n=4")
            dModelSyntax[model].append("q=1")
          if model[0]=="M" and len(model) > 2:
            dModelSyntax[model][0] = dModelSyntax[model][0][:-1]
            dModelSyntax[model].append("omegas=1")
          if model[:5]=="DFP07":
            dModelSyntax[model].append(["p0=1","p0=0.1"][model=="DFP07"])
	dLogLlh = {}		# dictionary(model:logllh)

	for model in lModels:
	  prevmodel, dnewpar = getNewParfromOptim(model, lModels, dModelLog, logger)
	  if prevmodel != "":
            fnew=open(dModelLog[model],"w")
            for k,v in dnewpar.items():
              fnew.write(k+"="+v.strip()+"\n")
            fnew.close()
            
            lignore = setIgnoreParams(model, prevmodel, lModels, logger)
            ignore = ",".join(lignore)
	  else:
	    ignore = ""
            
	  if model!="M0" and "M0" in lModels:
	    treeFile = dModelTrees["M0"]+"_1"

          # create dictionary with all elements of the two argument lists to build commands
	  modelDesc=dModelSyntax[model][0]+"("+",".join(dModelSyntax[model][1:])+")"
	  dBppCmd = {"INPUTFILE":alnFile, 
                   "FORMAT":alnFormat, 
                   "TREEFILE":treeFile, 
                   "MODEL":modelDesc, 
                   "IGNORE":ignore, 
                   "OUTTREE":dModelTrees[model], 
                   "OUTPARAMS":dModelParams[model], 
                   "BACKUP":dModelLog[model], 
                   "param":bppFile}

          # running bppml
	  logger.info("Running {:s} optimization".format(model))

          # join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
	  argsMx = "\""+"\" \"".join([k+"="+v for k, v in dBppCmd.items()])+"\""
	  logger.debug("bppml "+argsMx)
	  runMx = subprocess.Popen("bppml "+argsMx, shell=True, stdout=subprocess.PIPE).wait()
	  logger.debug(subprocess.PIPE)

	  # fill dictionary with loglikelihoods of each model
	  if os.path.exists(dModelParams[model]):
            with open(dModelParams[model], "r") as params:
              dLogLlh[model] = float(params.readline().strip().split("= ")[-1])
              params.close()
              logger.info("Log Likelihood = {}".format(dLogLlh[model]))
	  else:
            logger.info("Possible failed optimization, likelihood has not been calculated.")


	# perform LRT
        # M1 vs M2
	if "M1"  in lModels and "M2" in lModels:
              if "M1"  in dLogLlh and "M2" in dLogLlh:
                      LR12, p12 = PSPFunc.LRT(dLogLlh["M1"], dLogLlh["M2"], 2)
                      logger.info("LRT of M1 vs M2: {}".format(p12))
              else:
                      logger.info("Possible failed optimization, likelihoods of M1 and M2 have not been computed.")
	if "M7"  in lModels and "M8" in lModels:
              if "M7"  in dLogLlh and "M8" in dLogLlh:
                      LR78, p78 = PSPFunc.LRT(dLogLlh["M7"], dLogLlh["M8"], 2)
                      logger.info("LRT of M7 vs M8: {}".format(p78))
              else:
                      logger.info("Possible failed optimization, likelihoods of M7 and M8 have not been computed.")
	if "M8" in lModels and "M8a" in lModels:
              if "M8"  in dLogLlh and "M8a" in dLogLlh:
                      LR88a, p88a = PSPFunc.LRT(dLogLlh["M8a"], dLogLlh["M8"], 1)
                      ts88a = 0.5*p88a + 0.5
                      logger.info("LRT of M8 vs M8a: {} (Treshold: {})".format(p88a, ts88a))
              else:
                      logger.info("Possible failed optimization, likelihoods have not been computed.")
	if "DFP07"  in lModels and "DFP07_0" in lModels:
              if "DFP07"  in dLogLlh and "DFP07_0" in dLogLlh:
                      LRDFP, pDFP = PSPFunc.LRT(dLogLlh["DFP07_0"], dLogLlh["DFP07"], 1)
                      tsDFP = 0.5*pDFP + 0.5
                      logger.info("LRT of DFP07 vs DFP07_07: {} (Treshold: {})".format(pDFP, tsDFP))
              else:
                      logger.info("Possible failed optimization, likelihoods have not been computed.")

      # Bppmixedlikelihoods
	""" 
        Optimize tree and model using bppml
        Variables to include are
              INPUTFILE - alignement file
              FORMAT - format of the aln file (here, phyx)
              TREEFILE - tree file for the analyzed aln
              PARAMS - .params file from model optimization (bppml)
              OUTINFO - name of the results file (info about sites etc.)
	"""

	for model in lModels:
              # use tree optimized in M0 for each model
              if "M0" in lModels:
                      treeFile = dModelTrees["M0"]+"_1"
              else:
                      treeFile = dModelTrees[model]+"_1"

              if model in ["M0","DFP07_0"]:
                      continue

              # dictionary(model:results file name)
              dModelResults = {model:outSite+baseName+"_results_"+model+".log" for model in lModels}

              dMixCmd = {"INPUTFILE":alnFile, 
                         "FORMAT":alnFormat, 
                         "TREEFILE":treeFile, 
                         "params":dModelParams[model], 
                         "OUTINFO":dModelResults[model], 
                         "param":bppMixed}

              logger.info("Running mixed likelihoods with model {:s}".format(model))
              argsMx = "\""+"\" \"".join([k+"="+v for k, v in dMixCmd.items()])+"\""
              logger.debug("bppmixedlikelihoods "+argsMx)
              runMx = subprocess.Popen("bppmixedlikelihoods "+argsMx, shell=True, stdout=subprocess.PIPE).wait()
              logger.debug(subprocess.PIPE)


def getNewParfromOptim(model, lModels, dModelLog, logger):          ### new values for parameters
  # Use previous backup file (in order M0->M1->M2->M7->(M8a)->M8) to accelerate optimization
  # dictionary of equivalences of specific parameter names between models
  dequiv={}
  ## omega from M0->M1->M2->M7(->M8a)->M8 &M0->DFP07
  dequiv["omega"] = {"M1":{"YNGP_M1.omega":"omega"},
		     "M2":{"YNGP_M2.omega0":"omega"},
		     "M0":{"YN98.omega":"omega"},
		     "M7":{"YNGP_M7.p":"[omega/(1-omega),1][omega==1]"},
		     "M8a":{"YNGP_M8.p":"[omega/(1-omega),1][omega==1]"},
		     "M8":{"YNGP_M8.p":"[omega/(1-omega),1][omega==1]"},
                     "DFP07_0":{"DFP07.omega":"omega"},
                     "DFP07":{"DFP07.omega":"omega"}}
  dequiv["p0"] = {"DFP07_0":{"DFP07.p0":"1"},
                  "DFP07":{"DFP07.p0":"0.1"}}

  dnewpar={}
  
  prevmodel = ""
  if not os.path.exists(dModelLog[model]):
    if model[0]=="M":
      for prevmodel in ["M8a","M7","M2","M1","M0"]:
        if not prevmodel in lModels or not os.path.exists(dModelLog[prevmodel]+".def"):
          prevmodel=""
        else:
          break
    elif model[:5]=="DFP07":
      for prevmodel in ["DFP07_0","M0"]:
        if not prevmodel in lModels or not os.path.exists(dModelLog[prevmodel]+".def"):
          prevmodel=""
        else:
          break


    if prevmodel!="":
      logger.info("Optimization for model " + model + " uses optimized parameters from model " + prevmodel)  
      fprev=open(dModelLog[prevmodel]+".def","r")
      lprev=list(fprev.readlines())
      fprev.close()

      dprevpar={l[:l.find("=")]:l[l.find("=")+1:] for l in lprev}
      if prevmodel=="M8a":  # Trick to get right parameters prefix
        prevmodel="M8"
      if model=="M8a":  # Trick to get right parameters prefix
        model="M8"
        

      logger.info("Optimization for model " + model + " uses optimized parameters from model " + prevmodel)  
      # first copy all parameters
      for st,val in dprevpar.items(): 
        if prevmodel=="M0":
          if model[0]=="M":
            nst=st.replace("YN98","YNGP_"+model)
          else:
            nst=st.replace("YN98","DFP07")
        else:
            nst=st.replace(prevmodel,model)

        if not nst in dnewpar.keys():
          dnewpar[nst]=val

      # And then for specific parameters
      for key, par in dequiv.items():
        if model in par.keys() and prevmodel in par.keys():
          parav=par[prevmodel]
          parap=par[model]
          for oname,oval in dprevpar.items():
            ## look which oname is in equivalence list
            for kparav in parav.keys():
              if oname.startswith(kparav+"_"):
                for npar, nexp in parap.items():
                  nname=oname.replace(kparav,npar)
                  nval=str(eval(nexp.replace(key,oval).strip()))
                  if True:#not nname in dnewpar.keys():
                    dnewpar[nname]=nval

              #            break
              # write in backup file
      logger.debug(str(dnewpar))

  return prevmodel, dnewpar

def setIgnoreParams(model, prevmodel, lModels, logger):
  # if M0 optimization in models, use tree optimized in M0 for subsequent model optimizations
  lignore=[]
  if model!="M0" and "M0" in lModels:
    lignore.append("BrLen")
  
  if model == "M8a":
    lignore.append("YNGP_M8.omegas*")

  if model=="DFP07_0":
    lignore.append("DFP07.p0_1")
    
  # do not re-optimize root & equilibrium if done before
  if prevmodel!="":
    logger.info("Optimization for model " + model + " does not re-optimize root frequencies" )  
    lignore.append("Ancient")

    logger.info("Optimization for model " + model + " does not re-optimize equilibrium frequencies" )  
    lignore.append("*_Full.theta*")
            
  return lignore


def pamlSite(alnFile, treeFile, lModels, pamlParams, outDir, baseName, logger):

      tree = EvolTree(treeFile)
      os.mkdir(outDir+"paml_site/")
      tree.workdir = outDir+"paml_site/"
      tree.link_to_alignment(alnFile, "Fasta")
      logger.info("PAML codeml")

      dModelRun = {}
      for model in lModels:
        if model in ["M0","M1","M2","M7","M8"]:
          logger.info("Running {:s}".format(model))
          dModelRun[model] = tree.run_model(model)

      if "M1" in dModelRun and "M2" in dModelRun:
              p12 = tree.get_most_likely("M2", "M1")
              logger.info("LRT of M1 vs M2 = {}".format(p12))
      if "M7"  in dModelRun and "M8" in dModelRun:
              p78 = tree.get_most_likely("M8", "M7")
              logger.info("LRT of M7 vs M8 = {}".format(p78))
      if "M8a"  in dModelRun and "M8" in dModelRun:
              p88a = tree.get_most_likely("M8a", "M8")
              logger.info("LRT of M8 vs M8a = {}".format(p88a))


