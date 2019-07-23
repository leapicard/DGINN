import os, shlex, logging, subprocess, re
from AnalysisFunc import cmd

def bppBranch(OPBFile, outDir, baseName, alnFile, alnFormat, treeFile, logger):
	### BRANCH ANALYSIS: BIO++ ONE PER BRANCH
	lPSNodes = []
	
	logger.info("One Per Branch (BIO++)")
	logger.info("OPB parameter file: {:s}".format(OPBFile))
	
	outOPB = outDir+"branch_analysis/"
	if not os.path.exists(outOPB):
		os.makedirs(outOPB)
	
	outFileName = outOPB+baseName
	outTree = outFileName+"_branch.dnd"
	outParams = outFileName+"_branch.params"
	outBackup = outFileName+"_branch_optimization"
	
	# create dictionary with all elements of the two argument lists to build commands
	dBppCmd = {"INPUTFILE": alnFile, "FORMAT": alnFormat, "TREEFILE": treeFile, "OUTTREE": outTree, "OUTPARAMS": outParams, "BACKUP": outBackup, "param": OPBFile}
		
	# running bppml
	logger.info("Running Branch optimization")
	# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
	argsOPB = "bppml \""+"\" \"".join([k+"="+v for k, v in dBppCmd.items()])+"\""
	logger.debug(argsOPB)
	runOPB = cmd(argsOPB, False)
	
	# parse branches under positive selection
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

	logger.info("Nodes under positive selection {}".format(lPSNodes))
	logger.info("Nodes under neutral or negative selection {}".format(lNSNodes))

	return lPSNodes

def bppBranchSite(GNHFile, lPSNodes, outDir, baseName, alnFile, alnFormat, treeFile, logger):

	### PSEUDO BRANCH-SITE ANALYSIS: BIO++ GENERAL NON HOMOLOGOUS
	
	logger.info("General Non Homogenous on branches with w > 1 (BIO++)")
	logger.info("GNH parameter file: {:s}".format(GNHFile))
	
	for node in lPSNodes:
		outGNH = outDir+"pseudo_branchsite_analysis/"
		if not os.path.exists(outGNH):
			os.makedirs(outGNH)
		
		outFileName = outGNH+baseName
		outTree = outFileName+"_pseudo_branchsite.dnd"
		outParams = outFileName+"_pseudo_branchsite.params"
		outBackup = outFileName+"_pseudo_branchsite_optimization"
		
		# create dictionary with all elements of the two argument lists to build commands
		
		dBppCmd = {"INPUTFILE":alnFile, "FORMAT":alnFormat, "TREEFILE":treeFile, "NODE":str(node), "OUTTREE":outTree, "OUTPARAMS":outParams, "BACKUP":outBackup, "param":GNHFile}
			
		# running bppml
		logger.info("Running Pseudo Branch-Site optimization")
		# join each couple of the cmd dictionary so that it reads "k1 = v1" "k2 = v2" etc...
		argsGNH = 'bppml "'+'" "'.join([k+"="+v for k, v in dBppCmd.items()])+'"'
		logger.debug(argsGNH)
		runGNH = cmd(argsGNH, False)


def memeBranchSite(aln, cladoFile, outDir, baseName, logger):

	### BRANCH-SITE ANALYSIS: HYPHY MEME
	
	logger.info("Branch-site (MEME, HYPHY)")
	
	outBSA = outDir+"branch_site/"
	if not os.path.exists(outBSA):
		subprocess.Popen("mkdir "+outBSA, shell =  True).wait()

	outMeme = open(outBSA+baseName+"_meme.out", "w")
	errMeme = open(outBSA+baseName+"_meme.err", "w")
	memeFile = outBSA+baseName+"_meme.bf"
	memeFile = memeFile.replace(":", "\:")

	with open(memeFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"Universal\";\n")
		bf.write("inputRedirect[\"02\"] = \"{:s}\";\n".format(aln.replace(":", "\:")))
		bf.write("inputRedirect[\"03\"] = \"{:s}\";\n".format(cladoFile.replace(":", "\:")))
		bf.write("inputRedirect[\"04\"] = \"All\";\n")
		bf.write("inputRedirect[\"05\"] = \"0.1\";\n")
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"SelectionAnalyses\" + DIRECTORY_SEPARATOR + \"MEME.bf\", inputRedirect);")

	logger.info("MEME batch file: {:s}".format(memeFile))
	
	# run MEME
	logger.debug("HYPHYMP "+memeFile)

	runBusted = subprocess.Popen("HYPHYMP "+memeFile, shell=True, stdout=outMeme, stderr=errMeme).wait()
	
	patterns = re.compile(r'\b\.bf\b|\b\.json\b')
	memeFileList = [ outDir+fn for fn in os.listdir(outDir) if re.search(patterns, fn) ]

	for fileM in memeFileList:
		fileM2 = fileM.replace(".fasta.", ".").replace(outDir, outBSA)
		os.rename(fileM, fileM2)
		if os.path.exists(fileM2):
			print("Output: "+fileM2)
	outMeme.close()
	errMeme.close()
