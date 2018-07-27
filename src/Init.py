import argparse, logging, os
import G2P_Object
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
	parser = argparse.ArgumentParser(prog=__file__, description="This Script simplify the commandline of G2P.")
	parser.add_argument('-dd', '--debug', dest='debug', action='store_true', help='enter verbose/debug mode')


	required = parser.add_argument_group('Mandatory input infos for running')
	required.add_argument('-in', '--inFile', metavar="<filename>", required=True, type=check, dest = 'inFile', help = 'Table of parameters mandatory to start the pipeline G2P.')

	return parser

def paramDef(infile):
	"""
	Check the parameters in the file.

	@param infile: path's file
	@return defaultParam: dico of parameters
	"""
	if not os.path.exists(infile):
		print("The file does not exist, try again:\n")

		while os.path.exist(infile):
			infile = raw_input("Path:\n")

	#Parsing
	lParams = ["infile", "blastdb", "outdir", "logfile", "evalue", "mincov", "perc_id", "step_id", "remote", "entry_query", "sptree", "API_Key", "gard", "treerecs", "psp", "alnfile", "treefile", "alnformat", "basename", "hyphySeuil", "busted", "meme", "models", "bppml", "mixedlikelihood", "opbFile", "gnhFile", "opb", "gnh"]
	with open(infile, "r") as content:
		dParams = {}
		for line in content:
			temp = line.strip("\n").split(":")
			if temp[0] not in lParams:
				if temp[0].startswith("#") == False:
					print(temp[0]+" is not a valide parameter.\n")
			else:
				dParams[temp[0]] = temp[1]

	#Verifying if parameters are correct
	if "step_id" not in dParams or dParams["step_id"] not in ["blast", "accession", "fasta", "orf", "prank", "phyml", "tree", "gard", "psp", ""]:
		print(dParams["step_id"]+" not available, it will be blast by default.\n")
		dParams["step_id"] = "blast"

	if "remote" not in dParams or dParams["remote"] not in ["True", "False", ""]:
		print("Remote option need to be a boolean, set to True by default\n")
		dParams["remote"] = "True"

	if "psp" not in dParams or dParams["psp"] not in ["True", "False", ""]:
		print("PSP will not be executed, set to False by default\n")
		dParams["psp"] = "False"
	elif dParams["psp"] == "True":
		if dParams["step_id"] == "psp":
			if "treefile" not in dParams or dParams["treefile"] == "":
				print("You need to precise the tree file because psp needed")
				exit()
			elif "treefile" not in dParams or dParams["alnfile"] == "":
				print("You need to precise the aln file because psp needed")
				exit()
		for opt in ["meme", "busted", "models", "bppml", "mixedlikelihood", "opbFile", "gnhFile"]:
			if opt not in dParams:
				dParams[opt] = ""
			elif opt not in ["meme", "busted", "models"]:
				if os.path.exists(dParams[opt].strip("\n")) == True:
					dParams[opt] = dParams[opt].strip("\n")
				else:
					path = dParams["infile"].split("/")[:-1]+opt+"_params.bpp"
					print(opt+" doesn't exist, new parameters file: "+path)
					PSPFunc.pspFileCreation(path, opt)
					dParams[opt] = path
			elif opt == "models":
				ltemp = []
				for M in dParams[opt].split(","):
					if M.strip(" ") not in ["M0", "M1", "M2", "M7", "M8", "M8a"]:
						print(M+" isn't a valid model.")
					else:
						ltemp.append(M)
				dParams[opt] = ",".join(ltemp)

			elif dParams[opt] not in ["True", "False"]:
				dParams[opt] = "False"

	elif dParams["step_id"] == "psp":
		print("Error: psp option set to false and step_id set to psp ...")
		exit()
	if dParams["step_id"] in ["blast","accession","fasta"]:
		if dParams["infile"] == "" or dParams["blastdb"] == "":
			print("Infile and Blastdb are necessary !\n")
			exit()

	#Creation of a dictionnary with all the parameters
	defaultParam = {"infile":"", "blastdb":"", "outdir":"", "logfile":"", "evalue":1e-3, "mincov":50, "perc_id":70, "entry_query":"", "API_Key":"", "sptree":"", "treerecs":"False", "gard":"False", "remote":"False", "step_id":"blast","psp":"False", "alnfile":"", "treefile":"", "alnformat":"Fasta", "basename":"", "hyphySeuil":0.05, "busted":"False", "meme":"False", "models":"", "bppml":"", "mixedlikelihood":"", "opbFile":"", "gnhFile":"", "opb":"False", "gnh":"False"}
	for i in defaultParam:
		if i in dParams:
			if dParams[i] != "":
				defaultParam[i] = dParams[i].strip("\n")


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
			args["logfile"] = args["infile"].split(".")[0]+"_G2P_"+timeStamp+".log"
		else:
			args["logfile"] = args["infile"].split(".")[0]+args["logfile"]
	else:
		if args["logfile"] == "":
			args["logfile"] = args["alnfile"].split(".")[0]+"_G2P_"+timeStamp+".log"
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

	mainData = G2P_Object.basicData(args["infile"], args["outdir"], args["blastdb"], timeStamp, args["sptree"], args["alnfile"], args["treefile"])

	if args["step_id"] == "tree" and args["treerecs"] == False:
		args["treerecs"] == True
	elif args["step_id"] == "gard" and args["gard"] == False:
		args["gard"] == True

	logger.info("Reading input file {:s}".format(mainData.CCDSFile))
	logger.info("Output directory: {:s}".format(mainData.o))
	logger.info("Analysis will begin at the {:s} step".format(args["step_id"].title()))

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
	i =0
	while step.lower() != lStep[i]:
			toDo.append(False)
			i+=1

	while len(toDo) != len(lStep):
			toDo.append(True)

	return toDo
