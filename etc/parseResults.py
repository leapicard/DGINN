#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package Parse Results
# @author Lea Picard



#################### MODULES ####################
## Python modules
import argparse, sys, os, numpy, pandas
from Bio import SeqIO
from collections import defaultdict, OrderedDict, Counter
from time import localtime, strftime
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
	stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
	LR = abs(2*(ll1-ll2))
	p = stats.chisqprob(LR, df)
	return(LR, p)

def countDupl(lFiles):
	lDupl = [f.split(".")[0]+".best.fas" for f in lFiles if "_D" in f.split("/")[-1] and not ":" in f.split("/")[-1]]
	
	if len(lDupl) == 1:
		return("{:d} duplication found.\n".format(len(lDupl)))
	elif len(lDupl) > 1:
		return("{:d} duplications found.\n".format(len(lDupl)))
	else:
		return("No duplication found.\n")

def resDupl(lFiles):
	lDupl = [f.split(".")[0]+".best.fas" for f in lFiles if "_D" in f.split("/")[-1] and not ":" in f.split("/")[-1]]
	for f in lFiles:
		if "remaining" in f.split("/")[-1] and ":" not in f.split("/")[-1]:
			lDupl.append(f.split(".")[0]+".best.fas")
	lRes = []
	dName2DuplNb = {}
	
	for aln in lDupl:
		dupl = SeqIO.parse(open(aln),'fasta')
		ids = [ortho.id.split("_")[1] for ortho in dupl]
		
		if "_part" in aln.split("/")[-1]:
			duplNb = "part"+aln.split("part")[1][0:1]
		elif "_D" in aln.split("/")[-1]:
			duplNb = aln.split("_")[-2].split(".")[0]
		elif "remaining" in aln.split("/")[-1]:
			duplNb = "remaining"
		
		dName2DuplNb[aln] = duplNb	
		mainName = max(set(ids), key=ids.count)
		lRes.append("{:s}: {:d} orthologs (mainly: {:s})\n".format(duplNb, len(ids) ,mainName))
	
	return(dName2DuplNb[aln])

def resRecomb(lFiles):
	lRecomb = [f.split(".")[0].split("/")[-1] for f in lFiles if ":" in f.split("/")[-1]]
	res = ""
	
	if len(lRecomb) > 0:
		lRecomb = sorted(lRecomb)	
		lNames = []
		
		for aln in lRecomb:
			if "mincov" in aln:
				lNames.append(aln.split("_")[5])
			elif "remaining" in aln:
				lNames.append(aln.split("_")[5])
			elif "_D" in aln:
				lNames.append(aln.split("_")[4])
		
		dBkpts = Counter(lNames)
		for event in dBkpts:
			res += "Recombination events identified in {:s} ({:d} breakpoint(s))\n".format(event, dBkpts[event] - 1)
	
	else:
		res = "No recombination events identified in any of the alignments.\n"
		
	return(res)
	
def ResBusted(baseName, posDir):
	WG = posDir+"/busted/"+baseName+"_busted.out"
	LRT = ""
	
	if os.path.exists(WG):
		with open(WG, "r") as wg:
			try:
				p = wg.read().split("**")[-2].replace(" ", "").split("=")[1]
				
				if float(p) < 0.05:
					LRT = "Y\t{}".format(p)
				else:
					LRT = "N\t{}".format(p)
				
				return("BUSTED\t{:s}\n".format(LRT))
			except:
				return("BUSTED\tna\n")
	else:
		return("BUSTED\tna\n")

def ResMeme(baseName, posDir):
	BS = posDir+"/meme/"+baseName+"_meme.out"
	MemePss = posDir+"/"+baseName+"_meme_pss.fasta"
	
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
				completeRes = "MEME\t\t\t{:d}\t{}\n".format(len(lRes), lRes)
				completeRes = completeRes.replace("[", "").replace("]", "")
				return(completeRes)
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
				return("MEME\t\t\t0\n")
		except:
			return("MEME\t\t\tna\n")
	else:
		return("MEME\t\t\tna\n")

def ResBppExtract(models, dLogLlh, dSAres, posDir, baseName):
	model1, model2 = models.split(" ")[0], models.split(" ")[1]
	PssFile = posDir+"/"+baseName+"_bpp"+model2+"_pss.fasta"
	
	if model1 and model2 in dLogLlh.keys() and isinstance(dLogLlh[model1], float) and isinstance(dLogLlh[model2], float):
		try:
			if model1 and model2 in dLogLlh:
				LR, p = LRT(dLogLlh[model1], dLogLlh[model2], 2)

				if p < 0.05 and os.path.exists(dSAres[model2]):
					df = pandas.read_csv(dSAres[model2], sep='\t')
					w = float(df.columns[-2].split("=")[-1])*0.8
					lRes1 = df[df.ix[:,-1]>w].ix[:,0].tolist()
					lRes2 = df[df.ix[:,-2]>0.95].ix[:,0].tolist()
					lRes = list(set(lRes1).intersection(lRes2))
					pR = "Y\t{}".format(p)
					
					completeRes = "Bpp{:s}{:s}\t{:s}\t{:d}\t{}\n".format(model1, model2, pR, len(lRes), lRes)
					completeRes = completeRes.replace("[", "").replace("]", "")
					return(completeRes)
					
					dId2pss = OrderedDict()
					for f in SeqIO.parse(open(aln),'fasta'):
						fID, fSeq = f.id, str(f.seq)
						fSeq2 = ""
						for i in lRes:
							fSeq2 += fSeq[i-1]
						dId2pss[fID] = fSeq2
					with open(PssFile, "w") as out:
						out.write(dict2fasta(dId2pss))
					
				else:
					pR = "N\t{}".format(p)
					return("Bpp{:s}{:s}\t{:s}\n".format(model1, model2, pR))
		except:
			return("Bpp{:s}{:s}\tna\n".format(model1, model2))
	else:
		return("Bpp{:s}{:s}\tna\n".format(model1, model2))

def ResBpp(baseName, posDir):
	dSA = {}
	dSAres = {}
	dLogLlh = OrderedDict()
	lModels = ["M1","M2","M7","M8"]
	
	for model in lModels:
		dSA[model] = posDir+"/bpp_site/"+baseName+"_"+model+".params"
		dSAres[model] = posDir+"/"+baseName+"_results"+model+".log"
		
		if os.path.exists(dSA[model]):
			with open(dSA[model], "r") as params:
				dLogLlh[model] = float(params.readline().strip().split("= ")[-1])
		
		else:
			dLogLlh[model] = "na"
				
	m1m2 = ResBppExtract("M1 M2", dLogLlh, dSAres, posDir, baseName)
	m7m8 = ResBppExtract("M7 M8", dLogLlh, dSAres, posDir, baseName)
	
	return(m1m2, m7m8, dLogLlh)

def ResPamlExtract(models, dModelLlh, dModelFile):
	model1, model2 = models.split(" ")[0], models.split(" ")[1]
	
	if model1 and model2 in dModelLlh and all(isinstance(value, float) for value in dModelLlh.values()):
		LR, p = LRT(dModelLlh[model1], dModelLlh[model2], 2)

		if p < 0.05:
			with open(dModelFile[model2], "r") as modFile:
				content = modFile.read()
				PSS = [int(line.strip().split(" ")[0]) for line in content.split("BEB")[1].split("The grid")[0].split(")")[-1].split("\n") if "*" in line]
			
			pR = "Y\t{}".format(p)
			completeRes = "Paml{:s}{:s}\t{:s}\t{:d}\t{}\n".format(model1, model2, pR, len(PSS), PSS)
			completeRes = completeRes.replace("[", "").replace("]", "")
			return(completeRes)
		
		else:
			pR = "N\t{}".format(p)
			return("Paml{:s}{:s}\t{:s}\n".format(model1, model2, pR))
	
	else:
		return("Paml{:s}{:s}\tna\n".format(model1, model2))

def ResPaml(posDir):
	PAML = posDir+"/paml_site/"
	lModels = ["M1","M2","M7","M8"]
	
	if os.path.exists(PAML):
		pamlDirs = sorted([os.path.join(PAML, dI) for dI in os.listdir(PAML) if os.path.isdir(os.path.join(PAML, dI))])
		dModelLlh = OrderedDict()
		dModelFile = {}
		lnL = 0

		for pDir in pamlDirs:
			for model in lModels:
				if model in pDir:
					dModelFile[model] = pDir+"/out"
					
					with open(pDir+"/rst1", "r") as modelFile:
						line = modelFile.read().strip()
						try:
							dModelLlh[model] = float(line.split("\t")[-1])
						except ValueError:
							dModelLlh[model] = "na"
								
		for model in lModels:
			if model not in dModelLlh.keys():
				dModelLlh[model] = "na"

		try:
			m1m2 = ResPamlExtract("M1 M2", dModelLlh, dModelFile)
		except:
			m1m2 = "PamlM1M2\tna\n"
		try:
			m7m8 = ResPamlExtract("M7 M8", dModelLlh, dModelFile)
		except:
			m7m8 = "PamlM7M8\tna\n"
			
		return(m1m2, m7m8, dModelLlh)
	
	else:
		return("PamlM1M2\tna\n", "PamlM7M8\tna\n", "na")

################### MAIN CODE ###################
if __name__ == "__main__":
	
	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This program outputs a summary of the results obtained through running DGINN on a list of genes of interest.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug', dest='debug', action='store_true', help='enter verbose/debug mode')

	filesReq = parser.add_argument_group('Mandatory input infos for running')
	filesReq.add_argument('-in', '--inFile', metavar="<filename>", type=check, required=True, dest = 'inFile', help =\
						'List of all the directories containing the results from DGINN analyses on different genes.')
	
	files = parser.add_argument_group('Optional input infos (default values)')
	files.add_argument('-o', '--outdir', metavar="<path/to/directory>", type=str, default="", required=False, dest = 'outDir', help =\
						'folder for analysis results (path - by default output file will be saved in the incoming directory)')	

	
	# Check parameters and get arguments
	args = parser.parse_args()
	inFile = args.inFile
	inDir = "/".join(inFile.split("/")[0:-1])
	outDir = args.outDir
	if outDir == "":
		outDir = inDir
	#debug = args.debug
	
	lDirs = [line.rstrip() for line in open(inFile)]
	dSub2Cut = OrderedDict({sub:"/".join(sub.split("/")[0:-1])[0:-21]+".best.fas" for sub in lDirs if os.path.exists("/".join(sub.split("/")[0:-1])[0:-21]+".best.fas")})
	
	resFile = inDir+"/summary.res"
	res = open(resFile, "w")
	res.write("FullName\tGene\tGeneSize\tOmega\tMethod\tPosSel\tPValue\tNbSites\tPSS\n")
	
	llhFile = inDir+"/llh_comp.tab"
	llh = open(llhFile, "w")
	llh.write("File\tMethod\tM1\tM2\tM7\tM8\n")
	
	# Positive selection
	#res.write("\nPOSITIVE SELECTION\n\n")
	
	for posDir, aln in dSub2Cut.items():
		#print(aln)
		baseName = aln.split("/")[-1].split(".")[0]
		M0file = posDir+"/bpp_site/"+baseName+"_M0.params"
		
		try:
			with open(M0file, "r") as M0:
				omega = M0.read().split("omega=")[1].split(")")[0]
		except:
			M0 = "na"
		
		bust = ResBusted(baseName, posDir)
		meme = ResMeme(baseName, posDir)
		bpp1v2, bpp7v8, bppLlh = ResBpp(baseName, posDir)
		paml1v2, paml7v8, pamlLlh = ResPaml(posDir)
		
		f = SeqIO.parse(open(aln),'fasta')
		ids = [fasta.id.split("_")[1].upper() for fasta in f]
		f = SeqIO.parse(open(aln),'fasta')
		lLen = [len(str(fasta.seq)) for fasta in f]
		
		if "frag" in baseName:
			add = "["+"-".join(baseName.split("frag")[1].split("to"))+"]"
			mainName = max(set(ids), key=ids.count)+add
		else:
			mainName = max(set(ids), key=ids.count)
		
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
		
		base = shortName+"\t"+mainName+"\t"+str(max(lLen)/3)+"\t"+omega+"\t"
		for x in [bust, bpp1v2, bpp7v8, paml1v2, paml7v8, meme]:
			res.write(base+x)
		
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
		
	
	res.close()
	llh.close()	
	print("Parsed results found in {:s}".format(resFile))
