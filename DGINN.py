#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package DGINN.py
# @author Lea Picard & Quentin Ganivet

"""
Main script called to run the pipeline.
"""

#################### MODULES ####################
import argparse, os, subprocess, logging, sys, shlex, time
from Bio import SeqIO
from pathlib import Path
from statistics import median
sys.path.insert(0, str((Path(__file__).parents[0] / "lib").resolve()))
import Init, BlastFunc, ExtractFunc, FastaResFunc, AnalysisFunc, \
	LoadFileFunc, PosSelFunc, TreeFunc, DataObject

#################### GLOBAL #####################
version="1.0"
VERSION_DATE='2019/11/06'

#################### MAIN CODE ##################
if __name__ == "__main__":

	#List of steps
	lSteps = ["blast", 
			  "accessions", 
			  "fasta", 
			  "orf", 
			  "alignment", 
			  "tree", 
			  "duplication", 
			  "recombination", 
			  "positiveSelection"]

	#Initializing necessary variables for pipeline execution
	parse = Init.parameterCommandLine(version, "DGINN")
	parameters = parse.parse_args()
	debug = parameters.debug
	hostfile = parameters.hostfile
	"""
	threads = int(parameters.threads)
	if threads < 2:
		threads = 2
	"""
	
	parameters = Init.paramDef(parameters.params, 
							   parameters.infile, 
							   parameters.queryName)
	Data, logger = Init.initLogger(parameters, 
								   debug, 
								   version)

	funcNeeded = Init.initPipeline(parameters["step"], lSteps)
	Data.sptree, parameters["duplication"] = TreeFunc.treeCheck(Data.sptree, 
															 parameters["duplication"], 
															 logger)
	dAlTree = {}
	Data.setGenAttr(parameters["step"])
	
	firstStep = ""
	if parameters["step"] in ["blast", "accessions", "fasta"]:
		firstStep = "orf"
	elif parameters["step"] in ["recombination", "positiveSelection"]:
		firstStep = ""
	else:
		firstStep = parameters["step"]
	
	for i in range(len(lSteps)):
		if funcNeeded[i] == True:
				
			if lSteps[i] == "blast":
				Data.baseName = LoadFileFunc.baseNameInit(Data.baseName, 
														  Data.queryFile, 
														  Data.aln, 
														  logger)		
				
				Data = BlastFunc.treatBlast(Data, 
											parameters["evalue"], 
											parameters["percID"], 
											parameters["mincov"], 
											parameters["APIKey"], 
											parameters["remote"], 
											parameters["entryQuery"])
			
			elif lSteps[i] == "accessions":
				if parameters["step"] == "accessions":
					Data = LoadFileFunc.accnEntry(Data)	

				ExtractFunc.treatAccns(Data)
				
			elif lSteps[i] == "fasta":
				if parameters["step"] == "fasta":
					Data = LoadFileFunc.getSeqEntry(Data, parameters["duplication"])

				FastaResFunc.fastaCreation(Data, 
										   logger, 
										   parameters["remote"], 
										   parameters["APIKey"], 
										   parameters["step"], 
										   parameters["duplication"])

			elif lSteps[i] == "orf":
				if parameters["step"] == "orf":
					Data = LoadFileFunc.orfEntry(Data, parameters["duplication"])

					if lSteps[i] == firstStep:
						LoadFileFunc.spTreeCheck(Data, 
												 firstStep, 
											 	 parameters["duplication"])	

				AnalysisFunc.orfFinder(Data)
				
			elif lSteps[i] == "alignment":
				if parameters["step"] == "alignment":
					Data = LoadFileFunc.prankEntry(Data, parameters["duplication"])

					if lSteps[i] == firstStep:
						LoadFileFunc.spTreeCheck(Data, 
												 firstStep, 
											 	 parameters["duplication"])	

				AnalysisFunc.alnPrank(Data, logger)
				fasCov = AnalysisFunc.covAln(Data.aln, 
											 parameters["mincov"], 
											 Data.queryName, 
											 Data.o)
				newAln = AnalysisFunc.runPrank(fasCov, 
											   Data.geneName, 
											   Data.o)
				Data.aln = newAln

			elif lSteps[i] == "tree":
				if parameters["step"] == "tree":
					Data = LoadFileFunc.phymlRecEntry(Data, logger)

					if lSteps[i] == firstStep:
						LoadFileFunc.spTreeCheck(Data, 
												 firstStep, 
											 	 parameters["duplication"])	

				dAlTree = AnalysisFunc.phyMLTree(Data, logger)
				dAlTree = AnalysisFunc.checkPhyMLTree(Data, 
													  dAlTree, 
													  logger)

			elif lSteps[i] == "duplication" and parameters["duplication"]:
				if parameters["step"] == "duplication":
					Data, dAlTree = LoadFileFunc.duplPSEntry(Data, logger)

					if lSteps[i] == firstStep:
						LoadFileFunc.spTreeCheck(Data, 
												 firstStep, 
											 	 parameters["duplication"])	

				dAlTree = TreeFunc.treeTreatment(Data, 
												 dAlTree, 
												 parameters["nbspecies"], 
												 logger)

			elif lSteps[i] == "recombination" and parameters["duplication"]:
				if parameters["step"] == "recombination":
					Data = LoadFileFunc.phymlRecEntry(Data, logger)
					dAlTree[Data.aln] = ""

				dAlTree = AnalysisFunc.gardRecomb(Data, 
												  parameters["hyphySeuil"], 
												  dAlTree, 
												  hostfile, 
												  logger)

			elif lSteps[i] == "positiveSelection" and parameters["positiveSelection"]:
				logger.info("Starting positive selection analyses.")
				
				if parameters["step"] == "positiveSelection":
					Data, dAlTree = LoadFileFunc.pspEntry(Data, 
														  parameters, 
														  logger)
			
				listArgsPosSel =  []
				for aln in dAlTree:
					listArgs = [Data, 
								parameters, 
								aln, 
								dAlTree[aln], 
								logger]
					listArgsPosSel.append(listArgs)
					
					with open(Data.o+"files_list.txt", "w") as fAT:
						fAT.write(aln+"\t"+dAlTree[aln])
						
					PosSelFunc.pspAnalysis(Data, 
										   parameters, 
										   aln, 
										   dAlTree[aln], 
										   logger)
				
				logger.info("Finished positive selection analyses.")	
				
