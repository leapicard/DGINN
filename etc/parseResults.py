#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package Parse Results
# @author Lea Picard

#################### MODULES ####################
## Python modules
import parseResFunc as PRS
import argparse, os
from time import localtime, strftime
from collections import OrderedDict
from Bio import SeqIO, AlignIO

#################################################
## Variables Globales
version="2.0"
VERSION_DATE='2020/09/16'


################### MAIN CODE ###################
if __name__ == "__main__":
	
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This program outputs a summary of the results obtained through running DGINN on a list of genes of interest.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug', dest='debug', action='store_true', help='enter verbose/debug mode')

	filesReq = parser.add_argument_group('Mandatory input infos for running')
	filesReq.add_argument('-in', '--inFile', metavar="<filename>", type=PRS.check, required=True, dest = 'inFile', help =\
						'List of all the directories containing the results from DGINN analyses on different genes, and their corresponding alignments.')
	
	files = parser.add_argument_group('Optional input infos (default values)')
	files.add_argument('-o', '--outdir', metavar="<path/to/directory>", type=str, default="", required=False, dest = 'outDir', help =\
						'folder for analysis results (path - by default output file will be saved in the incoming directory)')
	files.add_argument('-pr', '--postrate', metavar="<value>", type=float, default=0.95, required=False, dest = 'pr', help =\
						'folder for analysis results (path - by default output file will be saved in the incoming directory)')	

	
	# Check parameters and get arguments
	args = parser.parse_args()
	inFile = args.inFile
	inDir = "/".join(inFile.split("/")[0:-1])
	outDir = args.outDir
	pr = args.pr
	if outDir == "":
		outDir = inDir
	#debug = args.debug
	
	lDirs = [line.rstrip() for line in open(inFile)]
	dSub2Cut = OrderedDict({sub.split("\t")[0]:sub.split("\t")[1] for sub in lDirs})
	#dSub2Cut = OrderedDict({sub:"/".join(sub.split("/")[0:-1])[0:-21]+".best.fas" for sub in lDirs if os.path.exists("/".join(sub.split("/")[0:-1])[0:-21]+".best.fas")})

	timeStamp = strftime("%Y%m%d%H%M", localtime())
	resFile = inDir+"/DGINN_{}summary.tab".format(timeStamp)
	res = open(resFile, "w")
	head = ""
	allRes = []
	cov = inDir+"/DGINN_{}coverage.tab".format(timeStamp)
	dAlnCov = {}
	llhFile = inDir+"/DGINN_{}likelihoods.tab".format(timeStamp)
	llh = open(llhFile, "w")
	llh.write("File\tMethod\tM1\tM2\tM7\tM8\n")
	
	for posDir, aln in dSub2Cut.items():
	#posDir = "/home/lea/Documents/genes/2020_shScreen/TRIM69_CCDS_results_202004012008/TRIM69_sequences_filtered_longestORFs_mafft_mincov_prank_results_202004201724/positive_selection_results_202004201724/"
	#aln = "/home/lea/Documents/genes/2020_shScreen/TRIM69_CCDS_results_202004012008/TRIM69_sequences_filtered_longestORFs_mafft_mincov_prank.best.fas"

	#if os.path.exists(posDir):
		baseName = aln.split("/")[-1].split(".")[0]
		M0fileBpp = posDir+"/bpp_site/"+baseName+"_M0.params"
		M0filePaml = posDir+"/paml_site/M0/rst1"
		dGene = OrderedDict()

		try:
			with open(M0fileBpp, "r") as M0:
				wM0Bpp = str(M0.read().split("omega=")[1].split(")")[0])
		except:
			wM0Bpp = "na"

		try:
			with open(M0filePaml, "r") as M0:
				wM0Paml = str(M0.read().split("\t")[-4])

		except:
			wM0Paml = "na"

		baseName2 = baseName.split("_")[0]
		bust = PRS.ResBusted(baseName2, posDir)
		dGene.update(bust)
		meme = PRS.ResMeme(baseName2, posDir)
		dGene.update(meme)
		bpp1v2, bpp7v8, bppLlh = PRS.ResBpp(baseName2, posDir, pr)
		dGene.update(bpp1v2)
		dGene.update(bpp7v8)
		paml1v2, paml7v8, pamlLlh = PRS.ResPaml(posDir, pr)
		dGene.update(paml1v2)
		dGene.update(paml7v8)

		f = SeqIO.parse(open(aln),'fasta')
		ids = [fasta.id.split("_")[1].upper() for fasta in f]
		f = SeqIO.parse(open(aln),'fasta')
		sp = set([fasta.id.split("_")[0] for fasta in f])
		nbSp = len(sp)
		f = SeqIO.parse(open(aln),'fasta')
		lLen = [len(str(fasta.seq)) for fasta in f]
		alnLen = max(lLen)/3

		
		splitName = baseName.split("_")
		shortName = baseName.split("_")[0]

		for elt in splitName:
			if "gp" in elt:
				shortName += "_"+elt.replace("D", "").replace("gp", "-")
			if "part" in elt:
				shortName += "_"+elt
			if "mincov" in elt:
				shortName += "_all"
			if "remaining" in elt:
				shortName += "_rem"
			if ":" in elt:
				shortName += "_"+elt.replace("frag", "")

		if "frag" in baseName:
			add = "["+"-".join(baseName.split("frag")[1].split("to"))+"]"
			mainName = max(set(ids), key=ids.count)+add
			shortName += add
		else:
			mainName = max(set(ids), key=ids.count)
		
		lBase = ["File", "Name", "Gene", "GeneSize", "NbSpecies", "omegaM0Bpp", "omegaM0codeml"]
		head = "\t".join(lBase + list(dGene.keys()))
		lBaseRes = [baseName, shortName, mainName, str(alnLen), str(nbSp), wM0Bpp, wM0Paml]
		resLine = "\t".join(lBaseRes + list(dGene.values()))
		allRes.append(resLine)

		if type(bppLlh) is dict:
			bppLlhRes = [str(value) for value in bppLlh.values()]
			llh.write(baseName+"\tBPP\t"+"\t".join(bppLlhRes)+"\n")
		else:
			llh.write(baseName+"\tBPP\tna\n")
		if type(bppLlh) is dict:
			pamlLlhRes = [str(value) for value in pamlLlh.values()]
			llh.write(baseName+"\tPAML\t"+"\t".join(pamlLlhRes)+"\n")
		else:
			llh.write(baseName+"\tPAML\tna\n")

		dAlnCov[baseName+"\t"+shortName] = PRS.getCov(aln)

	dCovSort = {x:dAlnCov[x] for x in sorted(dAlnCov, key=lambda k: len(dAlnCov[k]), reverse = True)}
	with open(cov, "w") as outf:
		outf.write("\n".join([k+"\t"+"\t".join(map(str, v)) for k, v in dCovSort.items()]))

	res.write(head+"\n")
	res.write("\n".join(allRes))
	res.close()
	llh.close()	
	print("Parsed results found in {:s}".format(resFile))
