
import sys, os, logging
# get the path of the directory containing the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# construct the path of the module's parent directory
module_dir = os.path.join(script_dir)
# add the module's parent directory to the system path
sys.path.append(module_dir)

from FormatFunc import isCCDSFasta, isAln, isTree
from PSPFunc import pspFileCreation
from time import localtime, strftime
from Bio import SeqIO

version="1.0"

"""This file pools the necessary functions to initialize the pipeline."""

def paramDef(params):
	"""
	Check the parameters in the file.

	@param inf: path's file
	@return defaultParam: dico of parameters
	"""
	#If infile(s) given through command line, takes priority
	params["infile"] = list(map(str.strip,params["infile"].split(",")))

	#If list of file given, split and check what each file is
	if len(params["infile"]) > 1:
		for entryfile in params["infile"]:
			if isCCDSFasta(entryfile):
				params["queryFile"] = os.path.abspath(entryfile)
			if isAln(entryfile):
				params["alnfile"] = os.path.abspath(entryfile)
			if isTree(entryfile):
				params["treefile"] = os.path.abspath(entryfile)
	else:
		params["queryFile"] = os.path.abspath(params["infile"][0])
	
	if params["queryFile"] != "":
		params["infile"] = params["queryFile"]	
	elif params["alnfile"] != "":
		params["infile"] = params["alnfile"]	
	
	answers = ["Y", "YES", "T", "TRUE"]
	negAnswers = ["N", "NO", "F", "FALSE"]
	
	for param in params.keys():
		if type(params[param]) is not list:
			if str(params[param]).upper() in answers:
				params[param] = True
			elif str(params[param]).upper() in negAnswers:
				params[param] = False
		
	#Check if parameters are correct ["duplication","recombination","positiveSelection"] => ["option"]
	lSteps = ["blast", "accessions", "fasta", "orf", "alignment", "tree", "option"]
			  
	if params["step"] not in lSteps:
		print("Step \""+params["step"]+"\" not available, set to blast by default.")
		params["step"] = "blast"
	
	if params["remote"] == "":
		print("Remote option needs to be a boolean, set to True by default.")
		params["remote"] = True
	
	if "positiveSelection" not in params:
		print("Positive selection analyses will not be executed, set to False by default.")
		params["positiveSelection"] = False
		
	elif params["positiveSelection"]:
		if params["step"] == "option":
			if params["treefile"] == "":
				print("The pipeline requires a phylogenetic tree. Please provide one.")
				sys.exit()
			elif params["alnfile"] == "":
				print("The pipeline requires a codon alignment. Please provide one.")
				sys.exit()
		for opt in ["meme", "busted", "models", "paml", "bppml", "mixedlikelihood", "opb", "gnh"]:
			if opt not in ["meme", "busted", "models", "paml"]:
				if (type(params[opt]) is not bool) and os.path.exists(params[opt].strip("\n")):
					params[opt] = params[opt].strip("\n")
				elif type(params[opt]) is bool and params[opt]:
					path = params["outdir"]+opt+"_params.bpp"
					pspFileCreation(path, opt)
					params[opt] = path
					
			elif opt == "models":
				ltemp = []
				lmodelok = ["M0", "M1", "M2", "M7", "M8", "M8a", "M10", "DFP07_0", "DFP07"]
				for M in map(str.strip,params[opt].split(",")):
					if M == "":
						next
					elif M not in lmodelok and M[:-2] not in lmodelok:
						print(M + " isn't a valid model.")
					else:
						ltemp.append(M)
				params["models"] = ",".join(ltemp)

	# elif params["step"] == "positiveSelection":
	# 	print("Error: positiveSelection option set to false and step set to positiveSelection.")
	# 	sys.exit()
		
	if params["step"] in ["blast","accessions","fasta"]:
		if params["infile"] == "" or params["blastdb"] == "":
			print("Infile and Blastdb are necessary.")
			sys.exit()

	return params


def initLogger(data, args, debug, version):
	"""
	Function initializing pipeline logger for optimal monitoring.

	@param1 args: Object containing pipeline parameters
	@param2 debug: if the option debug is set
	@param3 version: Pipeline version
	@return1 mainData: Filled basicData object
	@return2 logger: Logging object
	"""
	
	## Log
	### Set up the log directory
	filename = args["infile"].split("data/")[1].split(".")[0]
	timeStamp = strftime("%Y%m%d%H%M", localtime())
	# if args["infile"] != "":
	# 	if args["logfile"] == "":
	# 		args["logfile"] = args["outdir"]+filename+"_DGINN_"+timeStamp+".log"
	# 	else:
	# 		args["logfile"] = args["outdir"]+filename+args["logfile"]
	# else:
	# 	args["logfile"] = args["outdir"]+"DGINN_error.log"
			
	data["queryFile"] = args["infile"]
	data["o"] = args["outdir"]
	data["db"] = args["blastdb"]
	data["lBlastRes"] = []
	data["sptree"] = args["sptree"]
	data["aln"] = args["alnfile"]
	data["tree"] = args["treefile"]
	data["queryName"] = args["queryName"]
	
	# if args["step"] == "duplication" and args["duplication"] == False:
	# 	args["duplication"] = True
	# elif args["step"] == "recombination" and args["recombination"] == False:
	# 	args["recombination"] = True

	firstStep = ""
	if args["step"] in ["blast", "accessions", "fasta"]:
		firstStep = "orf"
	elif args["step"] in ["option"]:
		firstStep = ""
	# else:
	# 	firstStep = args["step"]

	data["firstStep"] = firstStep

	return data

def setGenAttr(data,params):	
	#Set attributes from the queryFile.
	step = params["step"]
	queryFile = data["queryFile"]

	if step == "blast":
		accns = list(SeqIO.parse(open(queryFile),'fasta'))
		accn = accns[0]
		data["queryName"] = accn.id
		params["queryName"] = accn.id
	else:
		if data["queryName"] == "" or params["queryName"] == "":
			print("queryName is not provided, coverage check will not run.")




	
	


