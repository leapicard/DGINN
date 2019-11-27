import argparse, logging, os, sys
import FormatFunc, PSPFunc, DataObject
from time import localtime, strftime

"""This file pools the necessary functions to initialize the pipeline."""

def check(arg):	
	"""
	Function checking whether a path to a file exists.
	
	@param arg: Path to a file
	@return os.path.abspath(arg): Normalized/absolute version of the path
	"""

	if not os.path.exists(arg):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: arg does not exist, please enter valid argument.
		raise argparse.ArgumentTypeError("\"{0}\" does not exist, please enter valid argument.".format(arg))
	
	return os.path.abspath(arg)


def parameterCommandLine(version, __file__):
	"""
	Function describing the available parameters and options to run the pipeline, for the user.

	@param version: Pipeline version
	@return parser: Object containing all parameters selected by the user
	"""

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description="DGINN, a pipeline for the Detection of Genetic Innovations.")
	parser.add_argument('-dd', '--debug', dest='debug', action='store_true', help='Enter verbose/debug mode')


	args = parser.add_argument_group('Mandatory parameters')
	args.add_argument('-p', '--params', metavar="<filename>", required=True, type=check, dest = 'params', 
					  help = 'Mandatory file with all the parameters necessary to run the pipeline.')
	
	xargs = parser.add_argument_group('Optional parameters')
	xargs.add_argument('-i', '--infile', metavar="<filename>", required=False, dest = 'infile', default="", 
					  help = 'Path or list of paths (absolute or relative) to the file(s) needed to start the pipeline (if indicated, will take priority over the parameters file)')
	xargs.add_argument('-q', '--query', metavar="<string>", required=False, dest = 'queryName', default="", 
					  help = 'Full identifier of the query in the format SpeciesName_GeneName_GeneID (if indicated, will take priority over the parameters file)')
	xargs.add_argument('-host', '--hostfile', metavar="<filename>", required=False, dest = 'hostfile', default="", 
					  help = 'Path to cluster hostfile if needed for mpi process')
	#xargs.add_argument('-th', '--threads', metavar="<integer>", required=False, dest = 'threads', default=2, 
	#				  help = 'Number of threads to use for parallelization (cluster use only)')
	
	return parser


def paramDef(params, inf, queryName):
	"""
	Check the parameters in the file.

	@param infile: path's file
	@return defaultParam: dico of parameters
	"""
	
	if not os.path.exists(params):
		print("The provided parameter file does not exist, try again.")
		sys.exit()

	#Parsing
	lParams = ["infile", 
			   "queryName", 
			   "queryFile", 
			   "blastdb", 
			   "outdir", 
			   "logfile", 
			   "evalue", 
			   "mincov", 
			   "percID", 
			   "step", 
			   "remote", 
			   "entryQuery", 
			   "sptree", 
			   "APIKey", 
			   "recombination", 
			   "duplication", 
			   "nbspecies", 
			   "positiveSelection", 
			   "basename", 
			   "hyphySeuil", 
			   "busted", 
			   "meme", 
			   "models", 
			   "paml", 
			   "bppml", 
			   "mixedlikelihood", 
			   "opb", 
			   "gnh"]
			   
	with open(params, "r") as content:
		dParams = {}
		for line in content:
			if line.startswith("#"):
				pass
			else
				temp = line.strip("\n").split(":")
				if temp[0] not in lParams:
					print(temp[0]+" is not a valid parameter.\n")
				else:
					dParams[temp[0]] = temp[1].strip()
	
	#If infile(s) given through command line, takes priority
	if inf != "":
		dParams["infile"] = inf.split(",")
	else:
		dParams["infile"] = dParams["infile"].split(",")
		
	#Idem queryName
	if queryName != "":
		dParams["queryName"] = queryName
	else:
		dParams["queryName"] = dParams["queryName"]
	
	#If list of file given, split and check what each file is
	if len(dParams["infile"]) > 1:
		for entryfile in dParams["infile"]:
			if FormatFunc.isCCDSFasta(entryfile):
				dParams["queryFile"] = entryfile
			if FormatFunc.isAln(entryfile):
				dParams["alnfile"] = entryfile
			if FormatFunc.isTree(entryfile):
				dParams["treefile"] = entryfile
	else:
		dParams["queryFile"] = dParams["infile"][0]
	
	if "queryFile" in dParams.keys() and dParams["queryFile"] != "":
		dParams["infile"] = dParams["queryFile"]	
	elif "alnfile" in dParams.keys() and dParams["alnfile"] != "":
		dParams["infile"] = dParams["alnfile"]	
	
	answers = ["Y", "YES", "T", "TRUE"]
	negAnswers = ["N", "NO", "F", "FALSE"]
	
	for param in dParams.keys():
		if type(dParams[param]) is not list:
			if dParams[param].upper() in answers:
				dParams[param] = True
			elif dParams[param].upper() in negAnswers:
				dParams[param] = False
		
	#Check if parameters are correct
	lSteps = ["blast", 
			  "accessions", 
			  "fasta", 
			  "orf", 
			  "prank", 
			  "phyml", 
			  "duplication", 
			  "recombination", 
			  "positiveSelection", 
			  ""]
			  
	if "step" not in dParams or dParams["step"] not in lSteps:
		print("Step \""+dParams["step"]+"\" not available, set to blast by default.")
		dParams["step"] = "blast"
	if dParams["step"] == "":
		dParams["step"] = "blast"
	
	if "remote" not in dParams or dParams["remote"] == "":
		print("Remote option needs to be a boolean, set to True by default.")
		dParams["remote"] = True
	
	if "positiveSelection" not in dParams:
		print("Positive selection analyses will not be executed, set to False by default.")
		dParams["positiveSelection"] = False
		
	elif dParams["positiveSelection"]:
		if dParams["step"] == "positiveSelection":
			if "treefile" not in dParams or dParams["treefile"] == "":
				print("The pipeline requires a phylogenetic tree. Please provide one.")
				sys.exit()
			elif "alnfile" not in dParams or dParams["alnfile"] == "":
				print("The pipeline requires a nucleotide or, preferably, codon alignment. Please provide one.")
				sys.exit()
				
		for opt in ["meme", "busted", "models", "paml", "bppml", "mixedlikelihood", "opb", "gnh"]:
			if opt not in dParams:
				dParams[opt] = ""
				
			elif opt not in ["meme", "busted", "models", "paml"]:
				if type(dParams[opt]) is not bool and os.path.exists(dParams[opt].strip("\n")):
					dParams[opt] = dParams[opt].strip("\n")
				elif dParams[opt]:
					path = "/".join(dParams["infile"].split("/")[:-1])+"/"+opt+"_params.bpp"
					print("{:s} parameter file: {:s}".format(opt, path))
					PSPFunc.pspFileCreation(path, opt)
					dParams[opt] = path
					
			elif opt == "models":
				ltemp = []
				for M in dParams[opt].split(","):
					if M.strip(" ") == "":
						next
					elif M.strip(" ") not in ["M0", "M1", "M2", "M7", "M8", "M8a"]:
						print(M+" isn't a valid model.")
					else:
						ltemp.append(M)
				dParams[opt] = ",".join(ltemp)

	elif dParams["step"] == "positiveSelection":
		print("Error: positiveSelection option set to false and step set to positiveSelection.")
		sys.exit()
		
	if dParams["step"] in ["blast","accessions","fasta"]:
		if dParams["infile"] == "" or dParams["blastdb"] == "":
			print("Infile and Blastdb are necessary.")
			sys.exit()

	#Creation of a dictionnary with all the parameters
	defaultParam = {"infile":"", 
					"queryName":"", 
					"queryFile":"", 
					"blastdb":"", 
					"outdir":"", 
					"logfile":"", 
					"evalue":1e-3, 
					"mincov":50, 
					"percID":70, 
					"entryQuery":"", 
					"APIKey":"", 
					"sptree":"", 
					"duplication":False, 
					"nbspecies":8, 
					"recombination":False, 
					"remote":False, 
					"step":"blast",
					"positiveSelection":False, 
					"alnfile":"", 
					"treefile":"", 
					"alnformat":"Fasta", 
					"basename":"", 
					"hyphySeuil":0.05, 
					"busted":False, 
					"meme":False, 
					"models":"", 
					"paml":"", 
					"bppml":"", 
					"mixedlikelihood":"", 
					"opb":False, 
					"gnh":False}
	
	for i in defaultParam:
		if i in dParams.keys() and dParams[i] != "":
			defaultParam[i] = dParams[i]
	
	return defaultParam

def initLogger(args, debug, version):
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
	timeStamp = strftime("%Y%m%d%H%M", localtime())
	if args["infile"] != "":
		if args["logfile"] == "":
			args["logfile"] = args["infile"].split(".")[0]+"_DGINN_"+timeStamp+".log"
		else:
			args["logfile"] = args["infile"].split(".")[0]+args["logfile"]
	else:
		if args["logfile"] == "":
			args["logfile"] = args["alnfile"].split(".")[0]+"_DGINN_"+timeStamp+".log"
		else:
			args["logfile"] = args["alnfile"].split(".")[0]+args["logfile"]
			
	# create logger
	#logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger("main")
	logger.setLevel(logging.INFO)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(args["logfile"])
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	if debug:
		logger.setLevel(logging.DEBUG)
		ch.setLevel(logging.DEBUG)
		fh.setLevel(logging.DEBUG)
	else:
		ch.setLevel(logging.INFO)
		fh.setLevel(logging.INFO)
	#create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)
	# Welcome message
	logger.info("Starting {:s} (v{:s})".format(__file__, version))

	mainData = DataObject.basicData(args["infile"], 
									args["outdir"], 
									args["blastdb"], 
									timeStamp, 
									args["sptree"], 
									args["alnfile"], 
									args["treefile"], 
									args["queryName"])

	if args["step"] == "duplication" and args["duplication"] == False:
		args["duplication"] == True
	elif args["step"] == "recombination" and args["recombination"] == False:
		args["recombination"] == True

	logger.info("Reading input file {:s}".format(mainData.queryFile))
	logger.info("Output directory: {:s}".format(mainData.o))
	logger.info("Analysis will begin at the {:s} step".format(args["step"]))

	return mainData, logger


def initPipeline(step, lStep):
	"""
	Function initializing the list of functions to execute in the pipeline.

	@param1 step: Step at which to start the pipeline
	@param2 lStep: Full list of steps from the pipeline
	@return toDo: List of boolean associated with each step
	"""

	#Create a list of step to do
	toDo = []
	i = 0

	while step.lower() != lStep[i].lower():
			toDo.append(False)
			i+=1

	while len(toDo) != len(lStep):
			toDo.append(True)

	return toDo
