import logging, os, re, subprocess

def hyphyBusted(alnFile, cladoFile, outDir, baseName, logger):
### WHOLE-GENE ANALYSIS: HYPHY BUSTED
	logger.info("Whole Gene (BUSTED, HYPHY)")
		
	# create BUSTED batch file for the gene
	bustedFile = outDir+baseName+"_busted.bf"
	with open(bustedFile, "w") as bf:
		bf.write("inputRedirect = {};\n")
		bf.write("inputRedirect[\"01\"] = \"Universal\";\n")
		bf.write("inputRedirect[\"02\"] = \"{:s}\";\n".format(alnFile.replace(":", "\:")))
		bf.write("inputRedirect[\"03\"] = \"{:s}\";\n".format(cladoFile.replace(":", "\:")))
		bf.write("inputRedirect[\"04\"] = \"All\";\n")
		bf.write("inputRedirect[\"05\"] = \"\";\n")
		bf.write("ExecuteAFile(HYPHY_LIB_DIRECTORY + \"TemplateBatchFiles\" + DIRECTORY_SEPARATOR + \"SelectionAnalyses\" + DIRECTORY_SEPARATOR + \"BUSTED.bf\", inputRedirect);")
	logger.info("Batch file: {:s}".format(bustedFile))
	
	# run BUSTED
	### find way to suppress stderr
	outWG = outDir+"busted/"
	if not os.path.exists(outWG):
		subprocess.Popen("mkdir "+outWG, shell =  True).wait()
	outBusted = open(outWG+baseName+"_busted.out", "w")
	errBusted = open(outWG+baseName+"_busted.err", "w")
	
	cmd = "HYPHYMP {:s}".format(bustedFile.replace(":", "\:"))
	logger.debug(cmd)
	
	runBusted = subprocess.Popen(cmd, shell=True, stdout=outBusted, stderr=errBusted).wait()
	
	# get BUSTED result files and move them to appropriate directory
	
	patterns = re.compile(r'\b\.bf\b|\b\.json\b')
	bustedFileList = [outDir+fn for fn in os.listdir(outDir) if re.search(patterns, fn)]
	
	for fileB in bustedFileList:
		fileB2 = fileB.replace(".fasta.", ".").replace(outDir, outWG)
		os.rename(fileB, fileB2)
		if os.path.exists(fileB2):
			logger.info("Output: "+fileB2)
	outBusted.close()
	errBusted.close()
		
