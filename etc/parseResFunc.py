#################### MODULES ####################
## Python modules
import argparse, sys, os, numpy, pandas
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict, Counter
from time import localtime, strftime
from statistics import median
from scipy import stats
import glob

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

def ResBusted(baseName, posDir):
	WG = glob.glob(posDir+"/busted/*_busted.out")
	posSel = ""
	d = OrderedDict()
	if len(WG)>= 1 and os.path.exists(WG[0]):
		with open(WG[0], "r") as wg:
			try:
				p = wg.read().split("**")[-2].replace(" ", "").split("=")[1]

				if float(p) < 0.05:
					posSel = "Y"
				else:
					posSel = "N"
				d["BUSTED"] = posSel
				d["BUSTED_p-value"] = str(p)
			except:
				d = {}

	if len(d) == 0:
		d["BUSTED"] = "na"
		d["BUSTED_p-value"] = "na"

	return(d)

def ResMeme(baseName, posDir):
	BS = glob.glob(posDir+"/meme/*_meme.out")
	d = OrderedDict({"MEME_NbSites":"0", "MEME_PSS":"na"})

	if len(BS)>=1 and os.path.exists(BS[0]):
		try:
			with open(BS[0], "r") as bs:
				BSlines = bs.readlines()
			
			lRes = []
			nbSites = 0
			for line in BSlines:
				if line.startswith("### ** Found"):
					nbSites = int(line.split("_")[1])
				
			if nbSites > 0:
				lRes = [int(line.split("|")[1].replace(" ", "")) for line in BSlines if line.startswith("|") and line.split("|")[1].replace(" ", "").isdigit()]

			if len(lRes) > 0:
				d["MEME_NbSites"] = str(len(lRes))
				d["MEME_PSS"] = "{}".format(lRes).replace("[", "").replace("]", "")

		except:
			next

	"""
	if len(d) == 0:
		d["MEME_PSS"] = "na"
		d["MEME_NbSites"] = 0
	"""
	return(d)

def ResBppExtract(models, dLogLlh, dSAres, posDir, baseName, pr):
	model1, model2 = models.split(" ")[0], models.split(" ")[1]
	method = "Bpp{:s}{:s}".format(model1, model2)
	d = OrderedDict({method:"na", method+"_p-value":"na", method+"_NbSites":"0", method+"_PSS":"na"})
	if isinstance(dLogLlh[model1], float) and isinstance(dLogLlh[model2], float):
		try:
			if model1 and model2 in dLogLlh:
				LR, p = LRT(dLogLlh[model1], dLogLlh[model2], 2)
				if p < 0.05 and os.path.exists(dSAres[model2][0]):
					df = pandas.read_csv(dSAres[model2][0], sep='\t')
					w = float(df.columns[-2].split("=")[-1])*0.8

					lRes1 = df[df.iloc[:,-1]>w].iloc[:,0].tolist()
					lRes2 = df[df.iloc[:,-2]>pr].iloc[:,0].tolist()
					lRes = list(map(lambda x:x+1,set(lRes1).intersection(lRes2)))
					if len(lRes)!=0:
					        lResFinal = str("{}".format(lRes).replace("[", "").replace("]", ""))
					else:
                                                lResFinal = "na"
					posSel = "Y"
					nbPSS = str(len(lRes))
				else:
					posSel = "N"
					lResFinal = "na"
					nbPSS = "0"

				
				d[method] = posSel
				d[method+"_p-value"] = str(p)
				d[method+"_NbSites"] = nbPSS
				d[method+"_PSS"] = lResFinal
					
		except:
			next
	return(d)

def ResBpp(baseName, posDir, pr):
	dSA = {}
	dSAres = {}
	dLogLlh = {}
	lModels = ["M1","M2","M7", "M8a", "M8", "DFP07_0", "DFP07"]

	for model in lModels:
		dSA[model] = glob.glob(posDir+"/bpp_site/*"+model+".params")
		dSAres[model] = glob.glob(posDir+"/bpp_site/*_results_"+model+".log")
		if len(dSA[model])>=1 and os.path.exists(dSA[model][0]):
			with open(dSA[model][0], "r") as params:
				dLogLlh[model] = float(params.readline().strip().split("= ")[-1])

		else:
			dLogLlh[model] = "na"

	m1m2 = ResBppExtract("M1 M2", dLogLlh, dSAres, posDir, baseName, pr)
	m7m8 = ResBppExtract("M7 M8", dLogLlh, dSAres, posDir, baseName, pr)
	m8am8 = ResBppExtract("M8a M8", dLogLlh, dSAres, posDir, baseName, pr)
	dfp = ResBppExtract("DFP07_0 DFP07", dLogLlh, dSAres, posDir, baseName, pr)
	
	return(m1m2, m7m8, m8am8, dfp, dLogLlh)

def ResPamlExtract(models, dModelLlh, dModelFile, pr):
	model1, model2 = models.split(" ")[0], models.split(" ")[1]
	method = "codeml{:s}{:s}".format(model1, model2)
	d = OrderedDict({method:"na", method+"_p-value":"na", method+"_NbSites":"0", method+"_PSS":"na"})

	if isinstance(dModelLlh[model1], float) and isinstance(dModelLlh[model2], float):
		LR, p = LRT(dModelLlh[model1], dModelLlh[model2], 2)

		if p < 0.05:
			with open(dModelFile[model2], "r") as modFile:
				content = modFile.read()
				res = content.split("BEB")[1].split("The grid")[0].split("SE for w")[-1].split("\n")
				res = list(filter(None, res))
				PSS = [int(line.strip().split(" ")[0]) for line in res if float(list(filter(None, line.split(" ")))[2].replace("*", "")) > pr]
			
			posSel = "Y"
			nbPSS = str(len(PSS))
			lResFinal = str("{}".format(PSS).replace("[", "").replace("]", ""))
		
		else:
			posSel = "N"
			lResFinal = "na"
			nbPSS = "0"

		d[method] = posSel
		d[method+"_p-value"] = str(p)
		d[method+"_NbSites"] = nbPSS
		d[method+"_PSS"] = lResFinal
	
	return(d)

def ResPaml(posDir, pr):
	PAML = posDir+"/paml_site/"
	lModels = ["M1","M2","M7","M8","M8a"]
	dModelLlh = OrderedDict({model:"na" for model in lModels})
	dModelFile = {model:"na" for model in lModels}
	#dModelOmega = {model:"na" for model in lModels}

	if os.path.exists(PAML):
		pamlDirs = OrderedDict({dI:os.path.join(PAML, dI) for dI in sorted(os.listdir(PAML)) if os.path.isdir(os.path.join(PAML, dI))})

		for pDir in pamlDirs:
			for model in lModels:
				if model in pamlDirs[pDir].split("/"):
					dModelFile[model] = pamlDirs[pDir]+"/out"
					#print(dModelFile)
					if os.path.exists(pamlDirs[pDir]+"/rst1"):
						with open(pamlDirs[pDir]+"/rst1", "r") as modelFile:
							line = modelFile.read().strip()
							try:
								dModelLlh[model] = float(line.split("\t")[-1])
								#dModelOmega[model] = float(line.split("\t")[-4])
							except ValueError:
								dModelLlh[model] = "na"
		
		for model in lModels:
			if model not in dModelLlh.keys():
				dModelLlh[model] = "na"

	m1m2 = ResPamlExtract("M1 M2", dModelLlh, dModelFile, pr)
	m7m8 = ResPamlExtract("M7 M8", dModelLlh, dModelFile, pr)
	m8am8 = ResPamlExtract("M8a M8", dModelLlh, dModelFile, pr)
			
	return(m1m2, m7m8, m8am8, dModelLlh)
	
	#else:
	#	return("PamlM1M2\tna\n", "PamlM7M8\tna\n", "na")

def getCov(fAln):
	lCov = []
	try:
	        aln = AlignIO.read(open(fAln, "r"), "fasta")
	except ValueError:
                print("Alignment does not fit "  + fAln)
                return
	nbSeq = len(aln)
	alnLen = aln.get_alignment_length()
	for i in range(0, alnLen):
		nb = [record.seq[i] for record in aln].count("-")
		lCov.append((nbSeq-nb)/nbSeq)

	return(lCov[1::3])
	#return(lCov)
