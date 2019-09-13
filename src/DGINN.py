#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package DGINN.py
# @author Lea Picard & Quentin Ganivet

"""
Main script called to run the pipeline.
"""

#################### MODULES ####################
import Init, BlastFunc, ExtractFunc, FastaResFunc, AnalysisFunc, LoadFileFunc, PosSelFunc, TreeFunc, DataObject
import argparse, os, subprocess, logging, sys, shlex, time
from Bio import SeqIO
from statistics import median

#################################################
## Variables Globales
version="1.0"
VERSION_DATE='2018/06/25'

################### MAIN CODE ###################
if __name__ == "__main__":

	#List of steps
	lSteps = ["blast", "accessions", "fasta", "orf", "prank", "phyml", "duplication", "recombination", "positiveSelection"]

	#Initializing necessary variables for pipeline execution
	parse = Init.parameterCommandLine(version, __file__)
	parameters = parse.parse_args()
	debug = parameters.debug
	hostfile = parameters.hostfile
	threads = int(parameters.threads)
	if threads < 2:
		threads = 2
		
	parameters = Init.paramDef(parameters.params, parameters.infile, parameters.queryName)
	Data, logger = Init.initLogger(parameters, debug, version)
	funcNeeded = Init.initPipeline(parameters["step"], lSteps)
	Data.sptree, parameters["treerecs"] = TreeFunc.treeCheck(Data.sptree, parameters["treerecs"], logger)
	Data.setGenAttr(parameters["step"])
	dAlTree = {}
	
	firstStep = ""
	if parameters["step"] in ["blast", "accessions", "fasta"]:
		firstStep = "orf"
	elif parameters["step"] in ["recombination", "positiveSelection"]:
		firstStep = ""
	else:
		firstStep = parameters["step"]
	
	for i in range(len(lSteps)):
		if funcNeeded[i] == True:
			if lSteps[i] == firstStep:
				print("once")
				LoadFileFunc.spTreeCheck(Data, firstStep, parameters["treerecs"])
				
			if lSteps[i] == "blast":
				Data.baseName = LoadFileFunc.baseNameInit(Data.baseName, Data.queryFile, Data.aln, logger)			
				
				Data = BlastFunc.treatBlast(Data, parameters["evalue"], parameters["percID"], parameters["mincov"], parameters["APIKey"], parameters["remote"], parameters["entryQuery"])
			
			elif lSteps[i] == "accessions":
				if parameters["step"] == "accessions":
					Data = LoadFileFunc.accnEntry(Data)	

				ExtractFunc.treatAccns(Data)
				
			elif lSteps[i] == "fasta":
				if parameters["step"] == "fasta":
					Data = LoadFileFunc.getSeqEntry(Data, parameters["treerecs"])

				FastaResFunc.fastaCreation(Data, logger, parameters["remote"], parameters["APIKey"], parameters["step"], parameters["treerecs"])

			elif lSteps[i] == "orf":
				if parameters["step"] == "orf":
					Data = LoadFileFunc.orfEntry(Data, parameters["treerecs"])

				AnalysisFunc.orfFinder(Data)
				
			elif lSteps[i] == "prank":
				if parameters["step"] == "prank":
					Data = LoadFileFunc.prankEntry(Data, parameters["treerecs"])

				AnalysisFunc.alnPrank(Data, logger)
				fasCov = AnalysisFunc.covAln(Data.aln, parameters["mincov"], Data.queryName, Data.o)
				newAln = AnalysisFunc.runPrank(fasCov, Data.geneName, Data.o)
				Data.aln = newAln

			elif lSteps[i] == "phyml":
				if parameters["step"] == "phyml":
					Data = LoadFileFunc.phymlEntry(Data, parameters["treerecs"], logger)

				dAlTree = AnalysisFunc.phyMLTree(Data, logger)
				dAlTree = AnalysisFunc.checkPhyMLTree(Data, dAlTree, logger)

			elif lSteps[i] == "duplication" and parameters["treerecs"]:
				if parameters["step"] == "duplication":
					Data, dAlTree = LoadFileFunc.treeEntry(Data, logger)

				dAlTree = TreeFunc.treeTreatment(Data, dAlTree, parameters["nbspecies"], logger)

			elif lSteps[i] == "recombination" and parameters["gard"]:
				if parameters["step"] == "recombination":
					Data, dAlTree = LoadFileFunc.gardEntry(Data, parameters, logger)

				dAlTree = AnalysisFunc.gardRecomb(Data, parameters["hyphySeuil"], dAlTree, hostfile, logger)

			elif lSteps[i] == "positiveSelection" and parameters["positiveSelection"]:
				logger.info("Starting positive selection analyses.")
				if parameters["step"] == "positiveSelection":
					Data, dAlTree = LoadFileFunc.pspEntry(Data, parameters, logger)
			
				listArgsPosSel =  []
				for aln in dAlTree:
					listArgsPosSel.append([Data, parameters, aln, dAlTree[aln], logger])
					with open(Data.o+"files_list.txt", "w") as fAT:
						fAT.write(aln+"\t"+dAlTree[aln])
						
					PosSelFunc.pspAnalysis(Data, parameters, aln, dAlTree[aln], logger)
				
				logger.info("Finished positive selection analyses.")	
				
