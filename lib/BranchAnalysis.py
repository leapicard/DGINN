import os, shlex, logging, subprocess, re
from AnalysisFunc import cmd
from subprocess import PIPE

def bppBranch(OPBFile, outDir, baseName, alnFile, alnFormat, treeFile, logger):
	### BRANCH ANALYSIS: BIO++ ONE PER BRANCH
	
	logger.info("One Per Branch (BIO++)")
	logger.info("OPB parameter file: {:s}".format(OPBFile))
	
	outOPB = outDir+"bpp_branch/"
	if not os.path.exists(outOPB):
		os.makedirs(outOPB)
	
	outFileName = outOPB+baseName
	outTree = outFileName+"_branch.dnd"
	outParams = outFileName+"_branch.params"
	outBackup = outFileName+"_branch_optimization"
	
	# create dictionary with all elements of the two argument lists to build commands
	dBppCmd = {"INPUTFILE": alnFile, 
		   "FORMAT": alnFormat, 
		   "TREEFILE": treeFile, 
		   "OUTTREE": outTree, 
		   "OUTPARAMS": outParams, 
		   "BACKUP": outBackup, 
		   "param": OPBFile}
		
	# running bppml
	logger.info("Running Branch optimization")
	# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
	argsOPB = "bppml \""+"\" \"".join([k+"="+v for k, v in dBppCmd.items()])+"\""
	logger.debug(argsOPB)
	runOPB = cmd(argsOPB, False)

	return(outParams)

def parseNodes(outParams, logger):
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
					logger.info("Node {:d} is interesting (w = {:f})".format(node, w))
					lPSNodes.append(node)
				else:
					lNSNodes.append(node)
		op.close()
	logger.info("Nodes under positive selection {}".format(lPSNodes))
	logger.info("Nodes under neutral or negative selection {}".format(lNSNodes))

	return lPSNodes

def bppBranchSite(GNHFile, lPSNodes, outDir, baseName, alnFile, alnFormat, treeFile, logger):

	### PSEUDO BRANCH-SITE ANALYSIS: BIO++ GENERAL NON HOMOLOGOUS
	
	logger.info("General Non Homogenous on branches with w > 1 (BIO++)")
	logger.info("GNH parameter file: {:s}".format(GNHFile))
	
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
		logger.info("Running Pseudo Branch-Site optimization")
		
		# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
		argsGNH = 'bppml "'+'" "'.join([k+"="+v for k, v in dBppCmd.items()])+'"'
		logger.debug(argsGNH)
		runGNH = cmd(argsGNH, False)


def memeBranchSite(aln, cladoFile, outDir, baseName, logger):

	### BRANCH-SITE ANALYSIS: HYPHY MEME
	
	logger.info("Episodic selection (MEME, HYPHY)")
	
	outBSA = outDir+"meme/"
	if not os.path.exists(outBSA):
		subprocess.Popen("mkdir "+outBSA, shell =  True).wait()

	dopt = {}        
	dopt["--output"] = outBSA+baseName+".MEME.json"
	dopt["--alignment"] = aln
	dopt["--tree"] = cladoFile
	
	lopt = " ".join([k + " " + v for k,v in dopt.items()])
        
	# run MEME
	logger.info("hyphy meme "+ lopt)

	fout = open(outBSA+baseName+"_meme.out","w")
	runMeme = subprocess.Popen("hyphy meme "+ lopt, shell = True, stdout = fout, bufsize=0).wait()
	fout.close()
        
	os.rename(dopt["--output"], outBSA+baseName+".MEME.json")
	
	
