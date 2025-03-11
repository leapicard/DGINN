import os, logging, subprocess
import PSPFunc
import AnalysisFunc

def bppSite(alnFile, treeFile, outDir, bppFile, bppMixed, lModels, logger):	
  # outDir=os.getcwd()+"/"  # used to debug
  ### SITE ANALYSIS: BIO++
  lModels = [m if (m[-2:]=="_C" or m.find("_G")!=-1) else m+"_G" for m in lModels] #gamma(n=4) distrib is default
  dlModels = {"C":[m[:-2] for m in lModels if m[-2:]=="_C"]}  #constant distrib
  dlModels["G"] = [m[:m.rfind("_")] for m in lModels if m[-2:]!="_C"]  #gamma distrib
  dallMod={"M0":0,"M1":1,"M2":2,"M7":3,"M8a":4,"M8":5,"M10":4,"DFP07_0":1,"DFP07":2}
  for k in ["C","G"]:
    dlModels[k].sort(key = lambda x : dallMod[x])

  logger.info("Bio++ Site Analysis")
  if len(dlModels["C"])!=0:
          logger.info("Models to be run with Constant rate: {:s}".format(", ".join(model for model in dlModels["C"])))
  if len(dlModels["G"])!=0:
          logger.info("Models to be run with Gamma rate: {:s}".format(", ".join(model for model in dlModels["G"])))
  logger.info("Bppml parameter file: {:s}".format(bppFile))

  ## Bppml
  """ 
  Optimize tree and model using bppml
  Variables to include are
    INPUTFILE - alignement file
    FORMAT - format of the aln file (here, phyx)
    TREEFILE - tree file for the analyzed aln
    MODEL - choose which model you want run on the data YNGP_M0 through 8, same models as PAML, and DFP07 models
    IGNORE - parameters to ignore for optimization, for example if one is fixed.
    OUTTREE - name of the optimized output tree
    OUTPARAMS - name of the output file summarizing parameters
    BACKUP - name of log file
  """

  # Bppml output file names - dictionaries that associate model number with output file name for the model
  outSite = outDir+"/bpp_site/"
  if not os.path.exists(outSite):
    subprocess.Popen("mkdir "+outSite, shell =  True).wait()

  dModelTrees = {model:outSite+model+".dnd" for model in lModels}
  dModelParams = {model:outSite+model+".params" for model in lModels}
  dModelLog = {model:outSite+"optimization_"+model for model in lModels}
  dModelSyntax = {}
  for k,lmod in dlModels.items():
    dModelSyntax[k]={model:["YNGP_"+model,"frequencies=F3X4","initFreqs=observed", "data=1"] for model in lmod if model[0]=="M"}		# dictionary model number - [MODEL name, MODEL arguments for bppml]
    dModelSyntax[k].update({model:[model[:5],"protmodel=JTT92", "frequencies=F3X4","initFreqs=observed", "data=1"] for model in lmod if model[:5]=="DFP07"})
# take into account the specificities of each model (number of classes n for example)
  for k,lmod in dlModels.items():
   for model in lmod:
      if model in ["M7","M8","M8a"]:
        dModelSyntax[k][model].append("n=4")
      if model in ["M7","M8"]:
        dModelSyntax[k][model].append("q=1")
      if model=="M8":
        dModelSyntax[k][model].append("omegas=2")
      if model in ["M10"]:
        dModelSyntax[k][model].append("nbeta=4")
        dModelSyntax[k][model].append("ngamma=4")
        
      if model[:5]=="DFP07":
        dModelSyntax[k][model].append(["p0=1","p0=0.1"][model=="DFP07"])
        
  dLogLlh = {}		# dictionary(model:logllh)

  for k,lmod in dlModels.items():
    dLogLlh[k]={}                
    for model in lmod:
      if model!="M0" and "M0" in lmod:
          treeFile = dModelTrees["M0_"+k]+"_1"

      ## already done
      #if os.path.exists(dBppCmd["OUTTREE"]+"_1"):
      #  logger.info("{:s} optimization already done because {:s} exists".format(model+"_"+k,treefile))
      #  break
      
      prevmodel, dnewpar = getNewParfromOptim(model+"_"+k, lModels, dModelLog, logger)
      if prevmodel != "":
          fnew=open(dModelLog[model+"_"+k],"w")
          for k2,v in dnewpar.items():
            fnew.write(k2+"="+v.strip()+"\n")
          fnew.close()
        
          lignore = setIgnoreParams(model+"_"+k, prevmodel, lModels, logger)
          ignore = ",".join(lignore)
      else:
          ignore = ""
            

      # create dictionary with all elements of the two argument lists to build commands
      modelDesc=dModelSyntax[k][model][0]+"("+",".join(dModelSyntax[k][model][1:])+")"
      distribDesc=["Constant()","Gamma(n=4)"][k=="G"]
      dBppCmd = {"INPUTFILE":alnFile, 
                      "FORMAT":"Fasta", 
                      "TREEFILE":treeFile, 
                      "MODEL":modelDesc, 
                      "DISTRIB":distribDesc, 
                      "IGNORE":ignore, 
                      "OUTTREE":dModelTrees[model+"_"+k], 
                      "OUTPARAMS":dModelParams[model+"_"+k], 
                      "BACKUP":dModelLog[model+"_"+k], 
                      "param":bppFile}

      ## invalidate scenario if not mixed
      if model=="M0":
        dBppCmd["scenario1"]=""
        
      # running bppml
      logger.info("Running {:s} optimization".format(model+"_"+k))
      
      # join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
      argsMx = "\""+"\" \"".join([k+"="+str(v) for k, v in dBppCmd.items()])+"\""
      logger.debug("bppml "+argsMx)
      runMx = subprocess.Popen("bppml "+argsMx, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
      logger.debug(runMx.communicate()[0].decode("utf-8"))

      # fill dictionary with loglikelihoods of each model
      if os.path.exists(dModelParams[model+"_"+k]):
        with open(dModelParams[model+"_"+k], "r") as params:
          dLogLlh[k][model] = float(params.readline().strip().split("=")[-1])
          params.close()
          logger.info("Log Likelihood = {}".format(dLogLlh[k][model]))
      else:
          logger.info("Possible failed optimization, likelihood has not been calculated.")


  # perform LRT
  for k,lModels in dlModels.items():
    if not k in dLogLlh:
                  continue
    # M1 vs M2
    if "M1"  in lModels and "M2" in lModels:
                if "M1"  in dLogLlh[k] and "M2" in dLogLlh[k]:
                        LR12, p12 = PSPFunc.LRT(dLogLlh[k]["M1"], dLogLlh[k]["M2"], 2)
                        logger.info("LRT of M1_%s vs M2_%s: %f"%(k,k,p12))
                else:
                        logger.info("Possible failed optimization, likelihoods of M1_%s and M2_%s have not been computed."%(k,k))
    if "M7"  in lModels and "M8" in lModels:
                if "M7"  in dLogLlh[k] and "M8" in dLogLlh[k]:
                        LR78, p78 = PSPFunc.LRT(dLogLlh[k]["M7"], dLogLlh[k]["M8"], 2)
                        logger.info("LRT of M7_%s vs M8_%s: %f"%(k,k,p78))
                else:
                        logger.info("Possible failed optimization, likelihoods of M7_%s and M8_%s have not been computed."%(k,k))
    if "M7"  in lModels and "M10" in lModels:
                if "M7"  in dLogLlh[k] and "M10" in dLogLlh[k]:
                        LR710, p710 = PSPFunc.LRT(dLogLlh[k]["M7"], dLogLlh[k]["M10"], 3)
                        logger.info("LRT of M7_%s vs M10_%s: %f"%(k,k,p710))
                else:
                        logger.info("Possible failed optimization, likelihoods of M7_%s and M10_%s have not been computed."%(k,k))
    if "M8" in lModels and "M8a" in lModels:
                if "M8"  in dLogLlh[k] and "M8a" in dLogLlh[k]:
                        LR88a, p88a = PSPFunc.LRT(dLogLlh[k]["M8a"], dLogLlh[k]["M8"], 1)
                        ts88a = 0.5*p88a + 0.5
                        logger.info("LRT of M8_%s vs M8a_%s: %f (Treshold: %f)"%(k, k, p88a, ts88a))
                else:
                        logger.info("Possible failed optimization, likelihoods of M8_%s and M8a_%s have not been computed."%(k,k))
    if "DFP07"  in lModels and "DFP07_0" in lModels:
                if "DFP07"  in dLogLlh[k] and "DFP07_0" in dLogLlh[k]:
                        LRDFP, pDFP = PSPFunc.LRT(dLogLlh[k]["DFP07_0"], dLogLlh[k]["DFP07"], 1)
                        tsDFP = 0.5*pDFP + 0.5
                        logger.info("LRT of DFP07_%s vs DFP07_0_%s:%f (Treshold: %f)"%(k,k,pDFP, tsDFP))
                else:
                        logger.info("Possible failed optimization, likelihoods of DFP07_%s and DFP07_0_%s have not been computed."%(k,k))

                        
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
  
  for k,lModels in dlModels.items():
    for model in lModels:
              # use tree optimized in M0 for each model
              if "M0" in lModels:
                      treeFile = dModelTrees["M0_"+k]+"_1"
              else:
                      treeFile = dModelTrees[model+"_"+k]+"_1"

              if model in ["M0","DFP07_0"]:
                      continue

              # dictionary(model:results file name)
              dModelResults = {model+"_"+k:outSite+"results_"+model+"_"+k+".log" for model in lModels}

              dMixCmd = {"INPUTFILE":alnFile, 
                         "FORMAT":"Fasta", 
                         "TREEFILE":treeFile, 
                         "param":",".join([dModelParams[model+"_"+k],bppMixed]),
                         "OUTINFO":dModelResults[model+"_"+k], 
                         }

              logger.info("Running mixed likelihoods with model {:s}".format(model+"_"+k))
              argsMx = "\""+"\" \"".join([k2+"="+v for k2, v in dMixCmd.items()])+"\""
              logger.info("bppmixedlikelihoods "+argsMx)
              runMx = subprocess.Popen("bppmixedlikelihoods "+argsMx, shell=True, stdout=subprocess.PIPE).wait()
              logger.debug(subprocess.PIPE)


def getNewParfromOptim(model, lModels, dModelLog, logger):          ### new values for parameters
  # Use previous backup file (in order M0->M1->M2->M7->(M8a)->M8) to accelerate optimization
  # dictionary of equivalences of specific parameter names between models
  dequiv={}
  distr=model[-2:]  # dist _C or _G
  ## omega from M0->M1->M2->M7(->M8a)->M8 &M0->DFP07
  ## omega from M0->M1->M2->M7->M10
  dequiv["omega"] = {"M1":{"YNGP_M1.omega":"omega"},
         "M2":{"YNGP_M2.omega0":"omega"},
         "M0":{"YN98.omega":"omega"},
         "M7":{"YNGP_M7.p":"[omega/(1-omega),1][omega==1]"},
         "M8a":{"YNGP_M8.p":"[omega/(1-omega),1][omega==1]"},
         "M8":{"YNGP_M8.p":"[omega/(1-omega),1][omega==1]"},
                      "DFP07_0":{"DFP07.omega":"omega"},
                     "DFP07":{"DFP07.omega":"omega"}}
  dequiv["p0"] = {"DFP07_0":{"DFP07.p0":"1"},
                  "DFP07":{"DFP07.p0":"0.1"}} #0.1 to avoid optim gap for p0=1

  dnewpar={}
  
  prevmodel = ""
  if not os.path.exists(dModelLog[model]):
    if model[0]=="M":
      for prevmodel in [x+distr for x in ["M10", "M8a","M7","M2","M1","M0"]]:
        if not prevmodel in lModels or not os.path.exists(dModelLog[prevmodel]+".def"):
          prevmodel=""
        else:
          break
    elif model[:5]=="DFP07":
      for prevmodel in [x+distr for x in ["DFP07_0","M0"]]:
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
      modeldeb=model[:-2]
      prevmodeldeb=prevmodel[:-2]

      # first copy all parameters
      for st,val in dprevpar.items(): 
        if prevmodeldeb=="M0":
          if model[0]=="M":
            nst=st.replace("YN98","YNGP_"+modeldeb)
          else:
            nst=st.replace("YN98","DFP07")
        else:
            nst=st.replace(prevmodeldeb,modeldeb)

        if not nst in dnewpar:
          dnewpar[nst]=val

      # And then for specific parameters
      for key, par in dequiv.items():
        if modeldeb in par and prevmodeldeb in par:
          parav=par[prevmodeldeb]
          parap=par[modeldeb]
          for oname,oval in dprevpar.items():
            ## look which oname is in equivalence list
            for kparav in parav:
              if oname.startswith(kparav+"_"):
                for npar, nexp in parap.items():
                  nname=oname.replace(kparav,npar)
                  nval=str(eval(nexp.replace(key,oval).strip()))
                  if True:#not nname in dnewpar:
                    dnewpar[nname]=nval

              #            break
              # write in backup file
      logger.debug(str(dnewpar))

  return prevmodel, dnewpar

def setIgnoreParams(model, prevmodel, lModels, logger):
  # if M0 optimization in models, use tree optimized in M0 for subsequent model optimizations
  lignore=[]
  distr=model[-2:]  # dist _C or _G
  
  if model=="DFP07_0"+distr:
    lignore.append("DFP07.p0_1")
    
  # do not re-optimize root & equilibrium if done before
  if prevmodel!="":
    logger.info("Optimization for model " + model + " does not re-optimize root frequencies" )  
    lignore.append("Ancient")

    logger.info("Optimization for model " + model + " does not re-optimize equilibrium frequencies" )  
    lignore.append("*_Full.theta*")

    logger.info("Optimization for model " + model + " does not re-optimize branch lengths" )  
    lignore.append("BrLen")

  return lignore

def assignCodemlCtl(alnFile, treeFile, outFile, model):
       codemlparam={"seqfile":alnFile,  "treefile":treeFile,  "outfile":outFile,  "aaDist":0,  "fix_blength":0,  "cleandata":0,  "clock":0,  "CodonFreq":2,  "getSE":0,  "icode":0,  "fix_kappa":0,  "kappa":2,  "Mgene":0,  "model":0,  "ncatG":8,  "noisy":0,  "NSsites":0,  "fix_omega":0,  "omega":0.7,  "RateAncestor":0,  "runmode":0,  "seqtype":1,  "Small_Diff":1e-6,  "verbose":2}

       if model=="M0":
         codemlparam["NSsites"]=0
       elif model=="M1":
         codemlparam["NSsites"]=1
       elif model=="M2":
         codemlparam["NSsites"]=2
       elif model=="M7":
         codemlparam["NSsites"]=7
       elif model=="M8":
         codemlparam["NSsites"]=8 
       elif model=="M8a":
         codemlparam["NSsites"]=8
         codemlparam["fix_omega"]=1
         codemlparam["omega"]=1
         
       return(codemlparam)

def fasta2paml(aln, outPhy):
  fin = open(aln, "r")
  dseq={}
  l = fin.readline()
  while l:
    l2 = l.replace("!", "N")
    if l2.startswith(">"):
      name=l2[1:].strip()
      dseq[name]=""
    else:
      dseq[name]+=l2.strip().upper()
    l = fin.readline()
      
  fin.close()

  lseq=len(dseq[list(dseq.keys())[0]])
  fout = open(outPhy, "w")
  fout.write('%d\t%d\n\n'%(len(dseq),lseq))

  for k,v in dseq.items():
    fout.write(k+"\n")
    fout.write(v+"\n")
  fout.close()
  
def pamlSite(alnFile, treeFile, outDir, pamlParams, lModels, logger):

      lModels = [m if (m[-2:]=="_C" or m.find("_G")!=-1) else m+"_C" for m in lModels] #constant distrib only available
      dlModels = {"C":[m[:-2] for m in lModels if m[-2:]=="_C"]}  #constant distrib
#      dlModels["G"] = [m[:m.rfind("_")] for m in lModels if m[-2:]!="_C"]  #gamma distrib

      dLogLlh = {}		# dictionary(model:logllh)
      outP=outDir+"/paml_site/"
      if not os.path.exists(outP):
        os.mkdir(outP)

      here=os.getcwd()
      outPhy=alnFile[:alnFile.rfind(".")]+".phy"
      fasta2paml(alnFile, outPhy)
      logger.info("Codeml Site Analysis")
      
      for k,lModels in dlModels.items():
       dLogLlh[k]={}
       outpaml=outP+k+"/"
       if not os.path.exists(outpaml):
        os.mkdir(outpaml)

       for model in lModels:
        if model in ["M0","M1","M2","M7","M8","M8a"]:
          logger.info("Running {:s}".format(model+"_"+k))
          outpamlmod=os.path.join(outpaml,model)
          if not os.path.exists(outpamlmod):
            os.mkdir(outpamlmod)

          codemlparam = assignCodemlCtl(outPhy, treeFile, os.path.join(outpamlmod,"result.txt"), model)
          parname = os.path.join(outpamlmod,"codeml.ctl")
          f=open(parname,"w")
          f.write("\n".join([k+" = "+str(v) for k,v in codemlparam.items()]))
          f.close()
          os.chdir(outpamlmod) #Necessary because codeml does not like long paths...
          runCodeml = subprocess.Popen("codeml codeml.ctl", shell =  True, stdout=subprocess.PIPE).wait()
          lgl=pamlGetLogL(os.path.join(outpamlmod,"result.txt"))
          if lgl!=None:
                  dLogLlh[k][model]=lgl
          logger.info("Log Likelihood = {}".format(dLogLlh[k].get(model,None)))
          os.chdir(here)
          
       for k,lModels in dlModels.items():
          if not k in dLogLlh:
            continue
          if "M1"  in lModels and "M2" in lModels:
                        if "M1"  in dLogLlh[k] and "M2" in dLogLlh[k]:
                          LR12, p12 = PSPFunc.LRT(dLogLlh[k]["M1"], dLogLlh[k]["M2"], 2)
                        else:
                          logger.info("Possible failed optimization, likelihoods of M1_%s and M2_%s have not been computed."%(k,k))
          if "M7" in lModels and "M8" in lModels:
                    if "M7"  in dLogLlh[k] and "M8" in dLogLlh[k]:
                      LR78, p78 = PSPFunc.LRT(dLogLlh[k]["M7"], dLogLlh[k]["M8"], 2)
                      logger.info("LRT of M7_%s vs M8_%s: %f"%(k,k,p78))
                    else:
                      logger.info("Possible failed optimization, likelihoods of M7_%s and M8_%s have not been computed."%(k,k))
          if "M7"  in lModels and "M10" in lModels:
                      if "M7"  in dLogLlh[k] and "M10" in dLogLlh[k]:
                        LR710, p710 = PSPFunc.LRT(dLogLlh[k]["M7"], dLogLlh[k]["M10"], 3)
                        logger.info("LRT of M7_%s vs M10_%s: %f"%(k,k,p710))
                      else:
                        logger.info("Possible failed optimization, likelihoods of M7_%s and M10_%s have not been computed."%(k,k))
          if "M8" in lModels and "M8a" in lModels:
                      if "M8"  in dLogLlh[k] and "M8a" in dLogLlh[k]:
                        LR88a, p88a = PSPFunc.LRT(dLogLlh[k]["M8a"], dLogLlh[k]["M8"], 1)
                        ts88a = 0.5*p88a + 0.5
                        logger.info("LRT of M8_%s vs M8a_%s: %f (Treshold: %f)"%(k, k, p88a, ts88a))
                      else:
                        logger.info("Possible failed optimization, likelihoods of M8_%s and M8a_%s have not been computed."%(k,k))


def pamlGetLogL(file):
#        """ Get LogL from rst file."""
#        try:
                f=open(file,"r")
                for l in f.readlines():
                        if l[:3]=="lnL":
                                x=float(l.split()[-2])
                                f.close()
                                return x
                f.close()
#        except:
#                return
