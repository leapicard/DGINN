#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package DGINN.py
# @author Lea Picard & Quentin Ganivet

"""
Main script called to run the pipeline.
"""

#################### MODULES ####################
import Init, BlastFunc, ExtractFunc, FastaResFunc, AnalysisFunc, LoadFileFunc, PosSelFunc, TreeFunc
import G2P_Object
import argparse, os, subprocess, logging, sys, shlex, time
from Bio import SeqIO
from statistics import median
from pathos.pools import ProcessPool

#################################################
## Variables Globales
version="1.0"
VERSION_DATE='2018/06/25'

################### MAIN CODE ###################
if __name__ == "__main__":

	#List of steps
	lSteps = ["blast", "extract", "getSequences", "orf", "prank", "phyml", "duplication", "recombination", "positiveSelection"]

	#Initializing necessary variables for pipeline execution
	parse = Init.parameterCommandLine(version, __file__)
	parameters = parse.parse_args()
	debug = parameters.debug
	hostfile = parameters.hostfile
	threads = int(parameters.threads)
	if threads < 2:
		threads = 2
	parameters = Init.paramDef(parameters.params, parameters.infile)
	Data, logger = Init.initLogger(parameters, debug, version)
	funcNeeded = Init.initPipeline(parameters["step_id"], lSteps)
	Data.sptree, parameters["treerecs"] = TreeFunc.treeCheck(Data.sptree, parameters["treerecs"], logger)
	dAlTree = {}
	
	for i in range(len(lSteps)):

		if funcNeeded[i] == True:

			if lSteps[i] == "blast":
				Data.setGenAttr()
				Data.baseName = LoadFileFunc.baseNameInit(Data.baseName, Data.CCDSFile, Data.aln, logger)			
				
				Data = BlastFunc.treatBlast(Data, parameters["evalue"], parameters["perc_id"], parameters["mincov"], parameters["API_Key"], parameters["remote"], parameters["entry_query"])

			elif lSteps[i] == "extract":
				if parameters["step_id"] == "extract":
					Data = LoadFileFunc.extractEntry(Data)	

				ExtractFunc.treatAccns(Data, logger)
				
			elif lSteps[i] == "getSequences":
				if parameters["step_id"] == "getSequences":
					Data = LoadFileFunc.getSeqEntry(Data, parameters["step_id"])

				FastaResFunc.fastaCreation(Data, logger, parameters["remote"], parameters["API_Key"], parameters["treerecs"])

			elif lSteps[i] == "orf":
				if parameters["step_id"] == "orf":
					Data = LoadFileFunc.orfEntry(Data, parameters["treerecs"])

				AnalysisFunc.orfFinder(Data, logger)
				
			elif lSteps[i] == "prank":
				if parameters["step_id"] == "prank":
					Data = LoadFileFunc.prankEntry(Data, parameters["treerecs"])

				AnalysisFunc.alnPrank(Data, logger)
				fasCov = AnalysisFunc.covAln(Data.aln, parameters["mincov"], Data.geneName, Data.o)
				newAln = AnalysisFunc.runPrank(fasCov, Data.geneName, Data.o)
				Data.aln = newAln

			elif lSteps[i] == "phyml":
				if parameters["step_id"] == "phyml":
					Data = LoadFileFunc.phymlEntry(Data, parameters["treerecs"], logger)

				dAlTree = AnalysisFunc.phyMLTree(Data, logger)

			elif lSteps[i] == "duplication" and parameters["treerecs"]:
				if parameters["step_id"] == "duplication":
					Data, dAlTree = LoadFileFunc.treeEntry(Data, logger)

				dAlTree = TreeFunc.treeTreatment(Data, dAlTree, parameters["nbspecies"], logger)

			elif lSteps[i] == "recombination" and parameters["gard"]:
				if parameters["step_id"] == "recombination":
					Data, dAlTree = LoadFileFunc.gardEntry(Data, parameters, logger)

				dAlTree = AnalysisFunc.gardRecomb(Data, parameters["hyphySeuil"], dAlTree, hostfile, logger)

			elif lSteps[i] == "positiveSelection" and parameters["positiveSelection"]:
				logger.info("Starting positive selection analyses.")
				if parameters["step_id"] == "positiveSelection":
					Data, dAlTree = LoadFileFunc.pspEntry(Data, parameters, logger)
			
				listArgsPosSel =  []
				for aln in dAlTree:
					listArgsPosSel.append([Data, parameters, aln, dAlTree[aln], logger])
					with open(Data.o+"files_list.txt", "w") as fAT:
						fAT.write(aln+"\t"+dAlTree[aln])
						
				if len(listArgsPosSel) > 1:
					pool = ProcessPool(nodes=threads)
					logger.info("Pool started on {:d} threads".format(threads))				
					results = pool.amap(PosSelFunc.pspAnalysis, listArgsPosSel)
				
					while not results.ready():
						time.sleep(5)
						
				else:
					PosSelFunc.pspAnalysis(Data, parameters, aln, dAlTree[aln], logger)
				
				logger.info("Finished positive selection analyses.")
				
				
