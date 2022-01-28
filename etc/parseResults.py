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
import sys, glob

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
						'Threshold rate of omega>1 to admit positive selected sites.')	

	
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
	lDirs = [line for line in lDirs if line!='']
	dSub2Cut = OrderedDict({sub.split("\t")[0]:sub.split("\t")[1:] for sub in lDirs})
	#dSub2Cut = OrderedDict({sub:"/".join(sub.split("/")[0:-1])[0:-21]+".best.fas" for sub in lDirs if os.path.exists("/".join(sub.split("/")[0:-1])[0:-21]+".best.fas")})

	timeStamp = strftime("%Y%m%d%H%M", localtime())
	resFile = inDir+"/DGINN_{}_summary.tab".format(timeStamp)
	res = open(resFile, "w")
	head = ""
	allRes = []
	cov = inDir+"/DGINN_{}_coverage.tab".format(timeStamp)
	dAlnCov = {}
	# llhFile = inDir+"/DGINN_{}_likelihoods.tab".format(timeStamp)
	# llh = open(llhFile, "w")
	# lmeth=["M1","M2","M7","M8a","M8","M10","DFP07_0","DFP07"]
	# llh.write("File\tMethod\t"+"\t".join(lmeth)+"\n")
	
	for posDir, [aln] in dSub2Cut.items():
		posDir=posDir.rstrip("/")
		baseName = aln.split("/")[-1].split(".")[0]
		repDir = "/".join(posDir.split("/")[:-1])
		print(aln,repDir)
		if (not aln.startswith(repDir)):
		        allF = [repDir+"/"+f for f in os.listdir(repDir) if f.endswith("fas") or f.endswith("fasta")]
		        if len(allF)!=0:
		                aln=max(allF, key=os.path.getctime)
		# M0fileBpp = glob.glob(posDir+"/bpp_site/*_optimization_M0_G.def")
		# M0filePaml = posDir+"/paml_site/M0/rst1"
		dGene = OrderedDict()
		# try:
		# 	with open(M0fileBpp[0], "r") as M0:
                #           l=M0.readline()
                #           while l:
                #             if l.find("omega")!=-1:
                #               wM0Bpp = l.split("=")[1].strip()
                #               break
                #             l=M0.readline()
		# except:
		# 	wM0Bpp = "na"

		# try:
		# 	with open(M0filePaml, "r") as M0:
		# 	  wM0Paml = str(M0.read().split("\t")[-4])

		# except:
		# 	wM0Paml = "na"

                ### Get the method specific results
		baseName2 = baseName.split("_")[0]

		bust = PRS.ResBusted(baseName, posDir)
		dGene.update(bust)

		meme = PRS.ResMeme(baseName, posDir)
		dGene.update(meme)

		resBpp = PRS.ResBpp(baseName, posDir, pr)
		for suff,dicres in resBpp.items():
                  for (m1,m2),val in dicres.items():
                    pref="Bpp"+suff+"_"+m1+":"+m2
                    for k,v in val.items():
                      dGene[pref+"_"+k]=v
                      
		resPAML = PRS.ResPaml(posDir, pr)
		for val in resPAML.values():
		    dGene.update(val)

		lsplitid=[fasta.id.split("_") for fasta in SeqIO.parse(open(aln),'fasta')]
		minlid=min(map(len, lsplitid))
		sp=set([l[0] for l in lsplitid])
		ids=[l[1].upper() for l in lsplitid]
		if (len(set(ids))!=len(lsplitid)):
		  ids=["_".join(l[1:]).upper() for l in lsplitid]
		  if len(set(ids))!=len(lsplitid):
                    print("Unable to find good ids in alignment file " + aln)
                    sys.exit()

		lLen = [len(str(fasta.seq)) for fasta in SeqIO.parse(open(aln),'fasta')]
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
		
		lBaseRes = [baseName, shortName, mainName, str(int(alnLen)), str(len(sp))]#, wM0Bpp, wM0Paml]
#		resLine = "\t".join(lBaseRes + list(dGene.values()))
		allRes.append((lBaseRes,dGene))

		# bppLlhRes = [str(bppLlh[k]) for k in lmeth]
		# llh.write(baseName+"\tBPP\t"+"\t".join(bppLlhRes)+"\n")

		# if type(bppLlh) is dict:
		# 	pamlLlhRes = [str(value) for value in pamlLlh.values()]
		# 	llh.write(baseName+"\tPAML\t"+"\t".join(pamlLlhRes)+2*"\tna"+"\n")
		# else:
		# 	llh.write(baseName+"\tPAML\tna\n")

		dAlnCov[baseName+"\t"+shortName] = PRS.getCov(aln)

	dCovSort = {x:dAlnCov[x] for x in sorted(dAlnCov, key=lambda k: len(dAlnCov[k]), reverse = True)}
	with open(cov, "w") as outf:
		outf.write("\n".join([k+"\t"+"\t".join(map(str, v)) for k, v in dCovSort.items()]))

        # filter out unused headers
	if len(allRes)==0:
	  sys.exit(0)
	
	lBase = ["File", "Name", "Gene", "GeneSize", "NbSpecies"]#, "omegaM0Bpp", "omegaM0codeml"]
	allK = allRes[0][1].keys()
	allKok=[]
	for k in allK:
	  ok=False
	  for (a,meth) in allRes:
	    if meth[k]!="na":
	      ok=True
	      break
	  if ok:
	    allKok.append(k)
	head="\t".join(lBase+allKok)
	res.write(head+"\n")
	for (a,meth) in allRes:
	  line = "\t".join(a +  [meth[k] for k in allKok])
	  res.write(line+"\n")
	  
	res.close()
	# llh.close()	
	print("Parsed results found in {:s}".format(resFile))
