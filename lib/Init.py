import argparse, logging, os, sys
import PSPFunc
from time import localtime, strftime

"""This file pools the necessary functions to initialize the pipeline."""

# --- Path functions ---


def check(arg):
    """
    Function checking whether a path to a file exists.

    @param arg: Path to a file
    @return os.path.abspath(arg): Normalized/absolute version of the path
    """

    if not os.path.exists(arg):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: arg does not exist, please enter valid argument.
        raise argparse.ArgumentTypeError(
            '"{0}" does not exist, please enter valid argument.'.format(arg)
        )

    return os.path.abspath(arg)


def parameterCommandLine(version, __file__):
    """
    Function describing the available parameters and options to run the pipeline, for the user.

    @param version: Pipeline version
    @return parser: Object containing all parameters selected by the user
    """

    # Parameters recovery
    parser = argparse.ArgumentParser(
        prog=__file__,
        description="DGINN, a pipeline for the Detection of Genetic Innovations.",
    )
    parser.add_argument(
        "-dd",
        "--debug",
        dest="debug",
        action="store_true",
        help="Enter verbose/debug mode",
    )

    args = parser.add_argument_group("Mandatory parameters")
    args.add_argument(
        "-p",
        "--params",
        metavar="<filename>",
        required=True,
        type=check,
        dest="params",
        help="Mandatory file with all the parameters necessary to run the pipeline.",
    )

    xargs = parser.add_argument_group("Optional parameters")
    xargs.add_argument(
        "-i",
        "--infile",
        metavar="<filename>",
        required=False,
        dest="infile",
        default="",
        help="Path or list of paths (absolute or relative) to the file(s) needed to start the pipeline (if indicated, will take priority over the parameters file)",
    )
    xargs.add_argument(
        "-q",
        "--query",
        metavar="<string>",
        required=False,
        dest="queryName",
        default="",
        help="Full identifier of the query in the format SpeciesName_GeneName_GeneID (if indicated, will take priority over the parameters file)",
    )
    xargs.add_argument(
        "-o",
        "--outdir",
        metavar="<path>",
        required=False,
        dest="outdir",
        default="",
        help="Path to the output directory (if indicated, will take priority over the parameters file)",
    )
    xargs.add_argument(
        "-host",
        "--hostfile",
        metavar="<filename>",
        required=False,
        dest="hostfile",
        default="",
        help="Path to cluster hostfile if needed for mpi process",
    )
    # xargs.add_argument('-th', '--threads', metavar="<integer>", required=False, dest = 'threads', default=2,
    # 				  help = 'Number of threads to use for parallelization (cluster use only)')

    return parser


def paramDef(dParams):
    """
    Check the parameters in the file.

    @param dParams: dico of parameters
    """

    # Set default values
    
    defaultParam = {
        "input": "",
        "queryName": "",
        "queryFile": "",
        "blastdb": "",
        "outdir": "",
        "logfile": "",
        "evalue": 1e-3,
        "mincov": 50,
        "percID": 70,
        "maxLen": "cutoff",
        "entryQuery": "",
        "APIKey": "",
        "phymlOpt": "",
        "sptree": "",
        "duplication": False,
        "LBopt": "cutoff",
        "nbspecies": 8,
        "recombination": False,
        "remote": False,
        "step": "blast",
        "positiveSelection": False,
        "alnfile": "",
        "treefile": "",
        "alnformat": "Fasta",
        "basename": "",
        "hyphySeuil": 0.05,
        "busted": False,
        "meme": False,
        "models": "",
        "paml": "",
        "bppml": "",
        "mixedlikelihood": "",
        "opb": False,
        "gnh": False,
    }

    for i in defaultParam:
        if not i in dParams.keys() or  dParams[i] == "" or dParams[i]==None:
          dParams[i] = defaultParam[i]
        elif type(dParams[i])!=type(defaultParam[i]):
          dParams[i]=type(defaultParam[i])(dParams[i])

    # homogenize boolean answers
    answers = ["Y", "YES", "T", "TRUE"]
    negAnswers = ["N", "NO", "F", "FALSE"]

    for param in dParams.keys():
        if isinstance(dParams[param], str):
            if dParams[param].upper() in answers:
                dParams[param] = True
            elif dParams[param].upper() in negAnswers:
                dParams[param] = False


    # Check if parameters are correct
    lSteps = [
        "blast",
        "accessions",
        "fasta",
        "orf",
        "alignment",
        "tree",
        "duplication",
        "recombination",
        "positiveSelection",
        "",
    ]

    if "remote" not in dParams or dParams["remote"] == "":
        print("Remote option needs to be a boolean, set to True by default.")
        dParams["remote"] = True

        
    if "positiveSelection" not in dParams:
        print(
            "Positive selection analyses will not be executed, set to False by default."
        )
        dParams["positiveSelection"] = False

    elif dParams["positiveSelection"]:
        for opt in [
            "meme",
            "busted",
            "models",
            "paml",
            "bppml",
            "mixedlikelihood",
            "opb",
            "gnh",
        ]:
            if opt not in dParams:
                dParams[opt] = ""

#             elif opt not in ["meme", "busted", "models", "paml"]:
#                 if type(dParams[opt]) is not bool and os.path.exists(
#                     dParams[opt].strip("\n")
#                 ):
#                     dParams[opt] = dParams[opt].strip("\n")
#                 elif dParams[opt]:
#                     path = dParams["outdir"] + "/" + dParams["quer + "_positive_selection.txt","r")
# (
#                         "/".join(dParams["input"].split("/")[:-1])
#                         + "/"
#                         + opt
#                         + "_params.bpp"
#                     )
#                     PSPFunc.pspFileCreation(path, opt)
#                     dParams[opt] = path

#             elif opt == "models":
#                 ltemp = []
#                 lmodelok = [
#                     "M0",
#                     "M1",
#                     "M2",
#                     "M7",
#                     "M8",
#                     "M8a",
#                     "M10",
#                     "DFP07_0",
#                     "DFP07",
#                 ]
#                 for M in map(str.strip, dParams[opt].split(",")):
#                     if M == "":
#                         next
#                     elif M not in lmodelok and M[:-2] not in lmodelok:
#                         print(M + " isn't a valid model.")
#                     else:
#                         ltemp.append(M)
#                 dParams["models"] = ",".join(ltemp)

    elif dParams["step"] == "positiveSelection":
        print(
            "Error: positiveSelection option set to false and step set to positiveSelection."
        )
        sys.exit()

    if dParams["step"] in ["blast", "accessions", "fasta"]:
        if dParams["blastdb"] == "":
            print("Blastdb is necessary.")
            sys.exit()


    return dParams


def initLogger(data, logfile, debug, version):
    """
    Function initializing pipeline logger for optimal monitoring.

    @param1 data: Data object containing pipeline parameters
    @param2 logfile: log file
    @param3 debug: if the option debug is set
    @param4 version: Pipeline version
    @return1 mainData: Filled basicData object
    """

    ## Log

    # create logger
    # logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("main")
    logger.setLevel(logging.INFO)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfile)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    if debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
        fh.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    # Welcome message
    logger.info("Starting {:s} (v{:s})".format(__file__, version))

    logger.info("Reading input file {:s}".format(data.queryFile))
    logger.info("Output directory: {:s}".format(data.o))


