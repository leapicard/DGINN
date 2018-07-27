#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package G2P.py
# @author Lea Picard & Quentin Ganivet

"""
Main script called to run the pipeline.
"""

#################### MODULES ####################
import Init, Blast, AccessionFunc, FastaResFunc, AnalysisFunc, LoadFileFunc, psp, treeFunc
import G2P_Object
import argparse, os, subprocess, logging, sys, shlex
from time import localtime, strftime
from statistics import median

#################################################
## Variables Globales
version="2.1"
VERSION_DATE='2018/06/25'

################### MAIN CODE ###################
if __name__ == "__main__":

	#List of steps
	lSteps = ["blast", "accession", "fasta", "orf", "prank", "phyml", "tree", "gard", "psp"]

	#Initialisation des variables néccéssaires au bon fonctionnement du pipline.
	parse = Init.parameterCommandLine(version, __file__)
	parameters = parse.parse_args()
	debug = parameters.debug
	parameters = Init.paramDef(parameters.inFile)
	Data, logger = Init.initLogger(parameters, debug, version)
	funcNeeded = Init.initPipeline(parameters["step_id"], lSteps)
	Data.sptree, parameters["treerecs"] = treeFunc.treeCheck(Data.sptree, parameters["treerecs"], logger)
	dAlTree = {}

	for i in range(len(lSteps)):

		if funcNeeded[i] == True:

			if lSteps[i] == "blast":

				try:
					SeqIO.parse(open(Data.CCDSFile), "fasta")
				except:
					logger.info("The InFile isn't a Fasta file.")
					exit()
				Data.setGenAttr()
				Data.baseName = baseNameInit(Data.baseName, Data.CCDSFile, Data.aln, logger)
				Data = Blast.treatBlast(Data, parameters["evalue"], parameters["perc_id"], parameters["mincov"], parameters["API_Key"], parameters["remote"], parameters["entry_query"])

			elif lSteps[i] == "accession":

				if parameters["step_id"] == "accession":

					Data = LoadFileFunc.accessionEntry(Data)	

				AccessionFunc.treatAccns(Data, logger)

			elif lSteps[i] == "fasta":

				if parameters["step_id"] == "fasta":

					Data = LoadFileFunc.fastaEntry(Data, parameters["step_id"])

				FastaResFunc.fastaCreation(Data, logger, parameters["remote"], parameters["API_Key"], parameters["treerecs"])


			elif lSteps[i] == "orf":

				if parameters["step_id"] == "orf":

					Data = LoadFileFunc.orfEntry(Data, parameters["treerecs"])

				AnalysisFunc.orfFinder(Data, logger)

			elif lSteps[i] == "prank":

				if parameters["step_id"] == "prank":

					Data = LoadFileFunc.prankEntry(Data, parameters["treerecs"])

				AnalysisFunc.alnPrank(Data, logger)

			elif lSteps[i] == "phyml":

				if parameters["step_id"] == "phyml":

					Data = LoadFileFunc.phymlEntry(Data, parameters["treerecs"], logger)

				dAlTree = AnalysisFunc.phyMLTree(Data, logger)

			elif lSteps[i] == "tree" and parameters["treerecs"] == "True":

				if parameters["step_id"] == "tree":

					Data, dAlTree = LoadFileFunc.treeEntry(Data, logger)

				dAlTree = treeFunc.treeTreatment(Data, dAlTree, logger)

			elif lSteps[i] == "gard" and parameters["gard"] == "True":

				if parameters["step_id"] == "gard":

					Data, dAlTree = LoadFileFunc.gardEntry(Data, parameters, logger)

				dAlTree = AnalysisFunc.gardRecomb(Data, parameters["hyphySeuil"], dAlTree, logger)

			elif lSteps[i] == "psp" and parameters["psp"] == "True":

				if parameters["step_id"] == "psp":

					Data, dAlTree = LoadFileFunc.pspEntry(Data, parameters, logger)

				psp.pspAnalysis(Data, parameters, dAlTree, debug, logger)
