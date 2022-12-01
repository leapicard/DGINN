from scipy import stats
from collections import OrderedDict
import sys, re, os, logging, ete3

def getParams(models, paml, bppml, mixed, Busted, Meme, opb, gnh):
	# Check analyses to be run and where the parameters file are

	dCtrls = {}
	lModels = []
	if models != "":
		if bppml not in ["", "False", False]:# and mixed not in ["", "False", False]:
			dCtrls["bppml"] = bppml
			dCtrls["bppmixedlikelihood"] = mixed
		if paml not in ["", "False", False]:
			dCtrls["paml"] = paml
		lModels = re.compile("\s*,\s*").split(models)
	if opb != "" and opb != False:
		dCtrls["OPB"] = opb
	if gnh != "" and gnh != False:
		dCtrls["GNH"] = gnh
	if Busted:
		dCtrls["BUSTED"] = ""
	if Meme:
		dCtrls["MEME"] = ""

	return dCtrls, lModels

def supBoot(outDir, baseName, treeFile):
	# Suppress bootstrap numbers from treeFile (necessary for HYPHY)
	cladoFile = outDir+baseName+"_clado.tree"
	t = ete3.Tree(treeFile)
	t.write(format=9, outfile=cladoFile)
	return cladoFile

def nbNode(treeFile):
	# count number of nodes in tree file
	with open(treeFile, "r") as tree:
		data = tree.read()
		nodes = str(data.count("(")+data.count(")"))
		#logger.info("There are {:s} nodes in the provided tree.".format(nodes))
		tree.close()
	return nodes

def LRT(ll1, ll2, df):
	"""
	Calculates likelihood ratio test between two models.
	:params ll1, ll2: likelihood of the two models studied
	:param df: degrees of freedom of difference between the two models
	"""
	stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
	LR = max(2*(ll2-ll1),0)
	p = stats.chisqprob(LR, df)
	return(LR, p)

def NHXTree(tree):
	"""
	Take a newick tree and outputs it in NHX format with branch numbering.
	"""
	pattern = re.compile(r"[,\)]")
	par=pattern.findall(tree)
	lm=pattern.split(tree)
	
	nhxpar=["[&&NHX:ND=%d]"%i+par[i] for i in range(len(par))]
	
	lc=["".join(x) for x in zip(lm,nhxpar)]
	sout="".join(lc)+";"
	
	return sout
	
def pspFileCreation(path, option):
  if True:
    dparams={}
    dparams["alphabet"] = "Codon(letter=DNA)"
    dparams["input.data1"] = "alignment(file=$(INPUTFILE), format=$(FORMAT), sites_to_use=all, max_gap_allowed=50%, max_unresolved_allowed=100%)"
    dparams["input.tree1"] = "user(file=$(TREEFILE), format=Newick)"
    dparams["phylo1"] = "Single(process=1, data=1)"
    if option == "mixedlikelihood":
      #dparams["params"] = "$(PARAMS)"
      dparams["output.likelihoods.file"] = "$(OUTINFO)"
    else:
      dparams["root_freq1"] = "F3X4(init=observed, data=1)"
      dparams["rate_distribution1"] = "$(DISTRIB)"
      dparams["optimization"] = "FullD(derivatives=Newton)"
      dparams["optimization.ignore_parameters"] = "$(IGNORE)"
      dparams["optimization.max_number_f_eval"] =  "1000"
      dparams["optimization.tolerance"] = "0.00001"
      dparams["output.tree.file"] = "$(OUTTREE)"
      dparams["output.tree.format"] = "Newick"
      dparams["output.estimates"] = "$(OUTPARAMS)"
      dparams["optimization.backup.file"] = "$(BACKUP)"
      if option == "bppml":
        dparams["model1"] = "$(MODEL)"
        dparams["process1"] = "Homogeneous(model=1, tree=1, rate=1, root_freq=1)"
        dparams["scenario1"] = "split(model=1)"
      elif option == "gnh":
        dparams["nonhomogeneous"] = "general"
        dparams["nonhomogeneous.number_of_models"] = "2"
        dparams["model1"] = "YNGP_M1(frequencies=F3X4,initFreqs=observed, data=1)"
        dparams["model2"] = "YNGP_M2(frequencies=F3X4,initFreqs=observed,kappa=YNGP_M1.kappa_1,omega=YNGP_M1.omega_1, data=1)"
        dparams["process1"] = "NonHomogeneous(model1=1, model1.nodes_id=$(NODES1), model2=2, model2.nodes_id$(NODES2), tree=1, rate=1, root=1)"
        dparams["optimization.ignore_parameters"] = "BrLen,*kappa*,*theta*,Ancient"

    with open(path, "w") as bppFile:
      bppFile.write("\n".join([par + " = " + val for par,val in dparams.items()]))
      bppFile.close()
