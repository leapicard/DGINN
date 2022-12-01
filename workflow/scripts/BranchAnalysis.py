import os, shlex, logging, subprocess, re
from AnalysisFunc import cmd
from subprocess import PIPE
from shutil import copyfile
from SiteAnalysis import getNewParfromOptim, setIgnoreParams
import PSPFunc

def bppBranch(OPBFile, outDir, baseName, alnFile, alnFormat, treeFile):
	### BRANCH ANALYSIS: BIO++ ONE PER BRANCH
	
	#logger.info("One Per Branch (BIO++)")
	#logger.info("OPB parameter file: {:s}".format(OPBFile))
	
	outOPB = outDir+"bpp_branch/"
	if not os.path.exists(outOPB):
		os.makedirs(outOPB)
	
	model = "M2"	

	outFileName = outOPB+baseName
	outTree = outFileName+"_"+model+".dnd"
	outParams = outFileName+"_"+model+".params"
	outBackup = outFileName+"_optimization_"+model
	
	# create dictionary with all elements of the two argument lists to build commands
	dBppCmd = {"INPUTFILE": alnFile, 
		   "FORMAT": alnFormat, 
		   "TREEFILE": treeFile, 
		   "OUTTREE": outTree, 
		   "OUTPARAMS": outParams, 
		   "BACKUP": outBackup, 
	           "model1":"YNGP_"+model+"(frequencies=F3X4, initFreqs=observed, data=1)",
		   "param": OPBFile,
                   "process1":"OnePerBranch(model=1, tree=1, rate=1, root_freq=1, shared_parameters=(*kappa, *Full.theta*))"}
	# running bppml
	#logger.info("Running Branch optimization")

        ### look for previous M0 optim
	outSite = outDir+"bpp_site/"
	outSiteFileName = outSite+baseName

	lModels=["M0", "M2"]
	dPrevModelLog = {model:outSiteFileName+"_optimization_"+model for model in lModels}
	prevmodel, dnewpar = getNewParfromOptim(model, lModels, dPrevModelLog) # prevmodel, dnewpar = getNewParfromOptim(model, lModels, dPrevModelLog, logger)
	lignore=[]
	if prevmodel!="":        
	  fnew=open(outBackup,"w")
	  for k,v in dnewpar.items():
	     fnew.write(k+"="+v.strip()+"\n")
	  fnew.close()
	lignore= setIgnoreParams(model, prevmodel, lModels) #lignore= setIgnoreParams(model, prevmodel, lModels, logger)
        
	dBppCmd["IGNORE"]=",".join(lignore)
         
	# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
	argsOPB = "bppml \""+"\" \"".join([k+"="+v for k, v in dBppCmd.items()])+"\""
	#logger.debug(argsOPB)
	runOPB = cmd(argsOPB, False)

        # test each branch
        # Scan all parameter names
	fback = open(dBppCmd["BACKUP"]+".def","r")
	dparam={param.split("=")[0]:float(param.split("=")[1]) for param in fback.readlines()}
	fback.close()

	valM2=float(dparam["f(x)"])
        
        # cp outParams for each branch with M2 replaced with M1
        # only for the where theta1 < 0.999 & theta2 < 0.999

	fparam = open(outParams, "r")
	lcmd=[l for l in fparam.readlines()]
	fparam.close()

        ## Look for correspondance  model_nb <-> node_id        
	lprocess=[l for l in lcmd if l[:7]=="process"][0]
	lid=lprocess.split(".nodes_id=(")
	cormodid={}
	for i in range(0,len(lid),2):
	  mod=int(lid[i][lid[i].rfind("l")+1:])
	  idi=int(lid[i+1][:lid[i+1].find(")"):])
	  cormodid[idi]=mod

        ## Compute lk for each node with theta1_mod * theta2_mod < 0.999:
	fresbranch=open(outFileName+"_branch.txt","w")
	fresbranch.write("Id\tomega2\tprop\tM2\tM1\tLR\tp\n")
	del(dBppCmd["process1"])
	del(dBppCmd["model1"])
	for idi,mod in cormodid.items():
	  if dparam["YNGP_M2.theta2_%d"%mod] * dparam["YNGP_M2.theta1_%d"%mod] >= 0.999:
	     continue
	  fback = open(outBackup+"_%d"%mod,"w")
	  [fback.write(key+"="+str(val)+"\n") for key,val in dparam.items() if key!="YNGP_M2.theta2_%d"%mod]
	  fback.write("YNGP_M2.theta2_%d=1\n"%mod)
	  fback.write("YNGP_M2.omega2_%d=1\n"%mod)
	  fback.close()

	  lignore2=lignore[:]+[key for key in dparam if key not in ["YNGP_M2.theta1_%d"%mod, "YNGP_M2.omega0_%d"%mod]] 
	  dBppCmd["IGNORE"]=",".join(lignore2)
	  dBppCmd["params"] = outParams
	  dBppCmd["OUTPARAMS"] = outParams+"_%d"%idi
	  dBppCmd["BACKUP"] = outBackup+"_%d"%mod 
          
	  argsOPB = "bppml \""+"\" \"".join([k+"="+v for k, v in dBppCmd.items()])+"\""
#	  logger.debug(argsOPB)
	  runOPB = cmd(argsOPB, False)

	  fback = open(dBppCmd["BACKUP"]+".def","r")
	  dparam2={param.split("=")[0]:float(param.split("=")[1]) for param in fback.readlines()}
	  fback.close()

	  valM1=float(dparam2["f(x)"])
	  LR, p = PSPFunc.LRT(valM2,valM1,2)
	  fresbranch.write("%d\t%f\t%f\t%f\t%f\t%f\t%f\n"%(idi,dparam["YNGP_M2.omega2_%d"%mod],(1-dparam["YNGP_M2.theta2_%d"%mod]) * (1-dparam["YNGP_M2.theta1_%d"%mod]),valM2,valM1,LR,p))
	  if p<0.05:
	    #logger.info("Node {:d} is interesting (w = {:f})".format(idi, dparam["YNGP_M2.omega2_%d"%mod]))
	fresbranch.close()
	return(outParams)

def parseNodes(outParams):
	# parse branches under positive selection
	lPSNodes = []
	lNSNodes = []

	with open(outParams, "r") as op:
		for line in op.readlines():
			if line.startswith("model") and line.endswith(")\n"):
				line = line.split("=")
				w = float(line[-1][:-2])
				node = int(line[0][5:])-1
				
				if w > 1:
					#logger.info("Node {:d} is interesting (w = {:f})".format(node, w))
					lPSNodes.append(node)
				else:
					lNSNodes.append(node)
		op.close()
	#logger.info("Nodes under positive selection {}".format(lPSNodes))
	#logger.info("Nodes under neutral or negative selection {}".format(lNSNodes))

	return lPSNodes

def bppBranchSite(GNHFile, lPSNodes, outDir, baseName, alnFile, alnFormat, treeFile):

	### PSEUDO BRANCH-SITE ANALYSIS: BIO++ GENERAL NON HOMOLOGOUS
	
	#logger.info("General Non Homogenous on branches with w > 1 (BIO++)")
	#logger.info("GNH parameter file: {:s}".format(GNHFile))
	
	for node in lPSNodes:
		outGNH = outDir+"bpp_gnh/"
		if not os.path.exists(outGNH):
			os.makedirs(outGNH)
		
		outFileName = outGNH+baseName
		outTree = outFileName+"_pseudo_branchsite.dnd"
		outParams = outFileName+"_pseudo_branchsite.params"
		outBackup = outFileName+"_pseudo_branchsite_optimization"
		
		# create dictionary with all elements of the two argument lists to build commands
		dBppCmd = {"INPUTFILE":alnFile, 
				   "FORMAT":alnFormat, 
				   "TREEFILE":treeFile, 
				   "NODE":str(node), 
				   "OUTTREE":outTree, 
				   "OUTPARAMS":outParams, 
				   "BACKUP":outBackup, 
				   "param":GNHFile}
			
		# running bppml
		#logger.info("Running Pseudo Branch-Site optimization")
		
		# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
		argsGNH = 'bppml "'+'" "'.join([k+"="+v for k, v in dBppCmd.items()])+'"'
		#logger.debug(argsGNH)
		runGNH = cmd(argsGNH, False)


def memeBranchSite(aln, cladoFile, outDir, baseName):

	### BRANCH-SITE ANALYSIS: HYPHY MEME
	
	#logger.info("Episodic selection (MEME, HYPHY)")
	
	outBSA = outDir+"meme/"
	if not os.path.exists(outBSA):
		subprocess.Popen("mkdir "+outBSA, shell =  True).wait()

	dopt = {}        
	dopt["--output"] = outBSA+baseName+".MEME.json"
	dopt["--alignment"] = aln
	dopt["--tree"] = cladoFile
	
	lopt = " ".join([k + " " + v for k,v in dopt.items()])
        
	# run MEME
	#logger.info("hyphy meme "+ lopt)

	fout = open(outBSA+baseName+"_meme.out","w")
	runMeme = subprocess.Popen("hyphy meme "+ lopt, shell = True, stdout = fout, bufsize=0).wait()
	fout.close()
        
	os.rename(dopt["--output"], outBSA+baseName+".MEME.json")
	
	
