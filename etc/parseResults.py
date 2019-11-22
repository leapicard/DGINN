#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package Parse Results
# @author Lea Picard



#################### MODULES ####################
## Python modules
import argparse, sys, os, numpy, pandas
from Bio import SeqIO
from collections import defaultdict, OrderedDict
from time import localtime, strftime
from multiprocessing import Pool
from statistics import median
from scipy import stats

#################################################
## Variables Globales
version="1.0"
VERSION_DATE='2019/03/10'

#################################################
## Functions

def dict2fasta(dico):
	txtoutput = ""
	for key, value in dico.items():
		txtoutput += ">{:s}\n{:s}\n".format(str(key),str(value))
	return(txtoutput)

def check(arg):	
	if not os.path.exists(arg):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: arg does not exist, please enter valid argument.
		raise argparse.ArgumentTypeError("\"{0}\" does not exist, please enter valid argument.".format(arg))
	
	return os.path.abspath(arg)

def LRT(ll1, ll2, df):
	"""
	Calculates likelihood ratio test between two models.
	:params ll1, ll2: likelihood of the two models studied
	:param df: degrees of freedom of difference between the two models
	"""
	LR = abs(2*(ll1-ll2))
	stats.chisqprob = lambda chisq, df: stats.chi2.sf(LR, df)
	p = stats.chisqprob(LR, df)
	return(LR, p)
	

################### MAIN CODE ###################
if __name__ == "__main__":
	
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This program creates a summary of DGINN's results.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	parser.add_argument('-dd', '--debug', dest='debug', action='store_true', help='enter verbose/debug mode')

	filesReq = parser.add_argument_group('Mandatory input infos for running')
	filesReq.add_argument('-in', '--inDir', metavar="<filename>", type=check, required=True, dest = 'inDir', help =\
						'Path to directory of results')
	
	files = parser.add_argument_group('Optional input infos (default values)')
	files.add_argument('-o', '--outdir', metavar="<path/to/directory>", type=str, default="", required=False, dest = 'outDir', help =\
						'folder for analysis results (path - by default output file will be saved in the incoming directory)')	

	
	# Check parameters and get arguments
	args = parser.parse_args()
	inDir = args.inDir
	outDir = args.outDir
	if outDir == "":
		outDir = inDir
	debug = args.debug
	
	resFile = inDir+"_summary.res"
	res = open(resFile, "w")
	lFiles = [os.path.join(inDir, file) for file in os.listdir(inDir) if file.endswith(".phylip_phyml_tree.txt")]
	
	# Duplication
	res.write("DUPLICATION\n")
	lDupl = [f.split(".")[0]+".best.fas" for f in lFiles if "_D" in f.split("/")[-1] and not ":" in f.split("/")[-1]]
	
	if len(lDupl) > 0:
		res.write("{:d} duplications found.\n".format(len(lDupl)))
	else:
		res.write("No duplication found.\n")
			
	# Recombination
	res.write("\nRECOMBINATION\n")
	lRecomb = [f.split(".")[0]+".best.fas" for f in lFiles if ":" in f.split("/")[-1]]
	if len(lRecomb) > 0:
		res.write("Recombination breakpoints identified.\n")
	else:
		res.write("No trace of recombination detected.\n")
		
	# Positive selection
	res.write("\nPOSITIVE SELECTION\n\n")
	subDirs = [os.path.join(inDir, dI) for dI in os.listdir(inDir) if os.path.isdir(os.path.join(inDir,dI))]
	dDirAln = {}
	print(subDirs)
	exit()
	
	for sub in subDirs:
		dDirAln = {sub:aln.split(".")[0]+".best.fas" for aln in lFiles if aln.split(".")[0] == sub[0:-21]}
		print(dDirAln)
		exit()
		aln = dDirAln[sub]
		posDir = [os.path.join(sub, dir) for dir in os.listdir(sub) if os.path.isdir(os.path.join(sub, dir))][0]
		baseName = aln.split("/")[-1].split(".")[0]
		geneName = baseName.split("_")[0]
		WG = posDir+"/whole_gene/"+baseName+"_busted.out"
		BS = posDir+"/branch_site/"+baseName+"_meme.out"
		dSA = {}
		dSAres = {}
		dLogLlh = {}
		
		res.write(baseName+"\n")
		
		if aln in lDupl:
			dupl = SeqIO.parse(open(aln),'fasta')
			ids = [ortho.id for ortho in dupl]
			duplNb = aln.split("_")[-2].split(".")[0]
			res.write("\n{:s}: {:d} orthologs\n".format(duplNb,len(ids)))
			res.write(", ".join(ids)+"\n")
			
		if aln in lRecomb:
			recomb = aln.split("_")[-1].split(".")[0][4:]
			res.write("\nRecombination fragment: {:s}\n".format(recomb))
		
		BppM2Pss = posDir+"/"+baseName+"_bppM2_pss.fasta"
		BppM8Pss = posDir+"/"+baseName+"_bppM8_pss.fasta"
		MemePss = posDir+"/"+baseName+"_meme_pss.fasta"
		
		res.write("\nMethod\tp-value LRT\tsites\n")
		
		if os.path.exists(WG):
			with open(WG, "r") as wg:
				try:
					LRTwg = wg.read().split("**")[-2].replace(" ", "")
					res.write("BUSTED\t{:s}\n".format(LRTwg))
				except:
					res.write("BUSTED\tna\n")
				
		for model in ["M1","M2","M7","M8"]:
			dSA[model] = posDir+"/site_analysis/"+baseName+"_"+model+".params"
			dSAres[model] = posDir+"/"+baseName+"_results"+model+".log"
			if os.path.exists(dSA[model]):
				with open(dSA[model], "r") as params:
					dLogLlh[model] = float(params.readline().strip().split("= ")[-1])
		
		if "M1" and "M2" in dLogLlh.keys():
			try:
				if "M1" and "M2" in dLogLlh:
					LR, p = LRT(dLogLlh["M1"], dLogLlh["M2"], 2)
					
					if p < 0.05 and os.path.exists(dSAres["M2"]):
						df = pandas.read_csv(dSAres["M2"], sep='\t')
						w = float(df.columns[-2].split("=")[-1])*0.8
						lRes = df[df.ix[:,-1]>w].ix[:,0].tolist()
						res.write("Bppml (M1 vs M2)\tp={:f}\t{:d} sites {}\n".format(p, len(lRes), lRes))
						
						dId2pss = OrderedDict()
						for f in SeqIO.parse(open(aln),'fasta'):
							fID, fSeq = f.id, str(f.seq)
							fSeq2 = ""
							for i in lRes:
								fSeq2 += fSeq[i-1]
							dId2pss[fID] = fSeq2
						with open(BppM2Pss, "w") as out:
							out.write(dict2fasta(dId2pss))
						
					else:
						res.write("Bppml (M1 vs M2)\tp={:f}\n".format(p))
			except:
				res.write("na")
				
		if "M7" and "M8" in dLogLlh.keys():
			try:
				if "M7" and "M8" in dLogLlh:
					LR, p = LRT(dLogLlh["M7"], dLogLlh["M8"], 2)
					
					if p < 0.05 and os.path.exists(dSAres["M8"]):
						df = pandas.read_csv(dSAres["M8"], sep='\t')
						w = float(df.columns[-2].split("=")[-1])*0.8
						lRes = df[df.ix[:,-1]>w].ix[:,0].tolist()
						res.write("Bppml (M7 vs M8)\tp={:f}\t{:d} sites {}\n".format(p, len(lRes), lRes))
						
						dId2pss = OrderedDict()
						for f in SeqIO.parse(open(aln),'fasta'):
							fID, fSeq = f.id, str(f.seq)
							fSeq2 = ""
							for i in lRes:
								fSeq2 += fSeq[i-1]
							dId2pss[fID] = fSeq2
						with open(BppM8Pss, "w") as out:
							out.write(dict2fasta(dId2pss))
						
					else:
						res.write("Bppml (M7 vs M8)\tp={:f}\n".format(p))
			except:
				res.write("na")
									
		if os.path.exists(BS):
			try:
				with open(BS, "r") as bs:
					BSlines = bs.readlines()
				
				lRes = []
				nbSites = 0
				for line in BSlines:
					if line.startswith("### ** Found"):
						nbSites = int(line.split("_")[1])
					
				if nbSites > 0:
					lRes = [int(line.split("|")[1].replace(" ", "")) for line in BSlines if line.startswith("|") and line.split("|")[1].replace(" ", "").isdigit()]
						
				if len(lRes) > 0:
					res.write("MEME\t\t{:d} sites (p<0.05) {}\n".format(len(lRes), lRes))
					dId2pss = OrderedDict()
					for f in SeqIO.parse(open(aln),'fasta'):
						fID, fSeq = f.id, str(f.seq)
						fSeq2 = ""
						for i in lRes:
							fSeq2 += fSeq[i-1]
						dId2pss[fID] = fSeq2
					with open(MemePss, "w") as out:
						out.write(dict2fasta(dId2pss))
				else:
					res.write("MEME\t\tNo positively selected sites detected\n")
			except:
				res.write("na")
		res.write("\n")
		
	res.close()
	print("Parsed results found in {:s}".format(resFile))
