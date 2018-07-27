import logging, os, re, subprocess

def hyphyBusted(dCtrls, alnFile, cladoFile, outDir, baseName, debug, logger):
### WHOLE-GENE ANALYSIS: HYPHY BUSTED
	if "BUSTED" in dCtrls:
		logger.info("Whole Gene (BUSTED, HYPHY)")
		logger.info("BUSTED parameter file: %s" %dCtrls["BUSTED"])
		
		# create BUSTED batch file for the gene
		bustedFile = outDir+baseName+"_busted.bf"
		with open(bustedFile, "w") as bf:
			bf.write("inputRedirect = {};\n")
			bf.write("inputRedirect[\"01\"] = \"Universal\";\n")
			bf.write("inputRedirect[\"02\"] = \"%s\";\n" %alnFile)
			bf.write("inputRedirect[\"03\"] = \"%s\";\n" %cladoFile)
			bf.write("inputRedirect[\"04\"] = \"All\";\n")
			bf.write("ExecuteAFile(\"%s\", inputRedirect);\n" %dCtrls["BUSTED"])
		logger.info("Batch file: %s" %bustedFile)
		
		# run BUSTED
		### find way to suppress stderr
		outWG = outDir+"whole_gene/"
		if not os.path.exists(outWG):
			subprocess.Popen("mkdir "+outWG, shell =  True).wait()
		outBusted = open(outWG+baseName+"_busted.out", "w")
		errBusted = open(outWG+baseName+"_busted.err", "w")
		
		if debug:
			logger.debug("HYPHYMP "+bustedFile)

		runBusted = subprocess.Popen("HYPHYMP "+bustedFile, shell=True, stdout=outBusted, stderr=errBusted).wait()
		
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
		