import BlastFunc, TreeFunc, FormatFunc
import os, sys, logging

"""
This file pools the necessary functions to enter in the pipeline at a specific step.
"""


def checkPath(path, attrNames):
    """
    Check if the path exist

    @param1 path: path's file
    @param2 attrNames: Name of the attribute
    """
    if not os.path.exists(path):
        print("The file with the {:s} information does not exist.".format(attrNames))
        sys.exit()


def baseNameInit(baseName, queryFile, aln, step=""):
    """
    Initialization of the attribut basename

    @param1 baseName: String
    @param2 queryFile: Path
    @param3 aln: Path
        @param4
    @return baseName: String
    """

    logger = logging.getLogger(".".join(["main", step]))
    if baseName == "":
        if queryFile != "":
            baseName = queryFile.split(".")[0].split("/")[-1]
        elif aln != "":
            baseName = aln.split(".")[0].split("/")[-1]
        else:
            logger.error("Basename can't be initialized.")
            sys.exit()
    return baseName


def filterData(sptree, filePath, o):
    """
    Function which execute functions to delete genes which aren't in the species tree.

    @param1 sptree: Path of a newick file
    @param2 filePath: Path of a Fasta file
    @param3 o: Path of a directory
    @return path: Path of a file
    """
    corsg = TreeFunc.assocFile(sptree, filePath, o)

    path = TreeFunc.supData(filePath, corsg, o)

    return path, corsg


def parseLFile(path):
    """
    Function which return a list of path from a file

    @param path: Path
    @return list of path
    """
    with open(path, "r") as listFile:
        listFile = listFile.readlines()
        listFile.close()

    return [i.strip("\n") for i in listFile]


##==============================================================================================================================


def accnEntry(queryFile):
    """
    Function handling start of the pipeline at the Extract step.

    @param Data: basicData object
    @return data: basicData object
    """
    if FormatFunc.isBlastRes(queryFile):
        blastRes = queryFile
        query, lBlastRes = BlastFunc.parseBlast(blastRes)
        baseName = baseNameInit(
            "", query, "", "accessions"
        )
    else:
        logger = logging.getLogger("main.accessions")
        logger.error(
            "The provided file is not a tabular output of Blast+, exiting DGINN."
        )
        sys.exit()

    return query, lBlastRes


def getSeqEntry(Data):
    """
    Function handling start of the pipeline at the Fasta step.

    @param Data: basicData object
    @return data: basicData object
    """
    if FormatFunc.isAccns(Data.queryFile):
        Data.accnFile = Data.queryFile
        Data.baseName = baseNameInit(
            Data.baseName, Data.queryFile, Data.accnFile, "fasta"
        )

        Data.lBlastRes = [i.strip("\n") for i in open(Data.accnFile, "r").readlines()]
    else:
        logger = logging.getLogger("main.fasta")
        logger.error(
            "Provided file is not a list of NCBI accessions, terminating DGINN."
        )
        sys.exit()

    return Data


def orfEntry(Data):
    """
    Function handling start of the pipeline at the orf step.

    @param1 Data: basicData object
    @return data: basicData object
    """

    if FormatFunc.isFasta(Data.queryFile):
        Data.seqFile = Data.queryFile
        Data.baseName = baseNameInit(Data.baseName, Data.queryFile, Data.aln, "orf")
    else:
        logger = logging.getLogger("main.orf")
        logger.error(
            "The provided file is not a fasta of nucleotide sequences, exiting DGINN."
        )
        sys.exit()

    return Data


def prankEntry(Data):
    """
    Function handling start of the pipeline at the Prank step.

    @param1 Data: basicData object
    @return data: basicData object
    """
    if FormatFunc.isFasta(Data.queryFile):
        Data.ORFs = Data.queryFile
        Data.baseName = baseNameInit(
            Data.baseName, Data.queryFile, Data.ORFs, "alignment"
        )
        try:
            with open(Data.ORFs) as orf:
                Data.geneName = orf.readline().split("_")[1]
                orf.close()
        except IndexError:
            with open(Data.ORFs) as orf:
                Data.geneName = orf.readline().strip()
                orf.close()

    else:
        logger = logging.getLogger("main.alignment")
        logger.error("Provided file is not a fasta of sequences, terminating DGINN.")
        sys.exit()

    return Data


def phymlRecEntry(Data, step="tree"):
    """
    Function handling start of the pipeline at the phyml step.

    @param1 Data: basicData object
    @param2 treeOption: Boolean
    @return Data: basicData object
    """

    if FormatFunc.isAln(Data.queryFile):
        Data.aln = Data.queryFile
        Data.ORFs = Data.queryFile
        Data.baseName = baseNameInit(Data.baseName, Data.queryFile, Data.aln, step)
        try:
            with open(Data.aln) as orf:
                Data.geneName = orf.readline().split("_")[1]
                orf.close()
        except IndexError:
            with open(Data.aln) as orf:
                Data.geneName = orf.readline().strip()
                orf.close()

    else:
        logger = logging.getLogger(".".join(["main", step]))
        logger.error(
            "Provided file is not a multiple sequence alignment, terminating DGINN."
        )
        sys.exit()

    return Data


def spTreeCheck(Data, firstStep, treeOption):
    if not hasattr(Data, "cor") and treeOption:
        if firstStep == "orf":
            aln = Data.seqFile
        elif firstStep == "alignment":
            aln = Data.ORFs
        elif firstStep == "tree" or firstStep == "duplication":
            aln = Data.aln

        if not os.path.exists(Data.sptree):
            Data.sptree, treeOption = TreeFunc.treeCheck(Data.sptree, aln, treeOption)

        if Data.sptree != "":
            aln2, corSG = filterData(Data.sptree, aln, Data.o)

            if firstStep == "orf":
                Data.seqFile = aln2
            elif firstStep == "alignment":
                Data.ORFs = aln2
            elif firstStep == "tree" or firstStep == "duplication":
                Data.aln = aln2

            setattr(Data, "cor", corSG)


def duplPSEntry(Data, step="duplication"):
    """
    Function handling start of the pipeline at the tree step.

    @param Data: basicData object
    @return Data: basicData object
    """
    logger = logging.getLogger(".".join(["main", step]))
    dico = {}
    if Data.aln != "" and Data.tree != "":
        logger.info("Alignement file: " + Data.aln)
        logger.info("Gene Tree file: " + Data.tree)
        dico[Data.aln] = Data.tree
        Data.ORFs = Data.aln
        Data.baseName = baseNameInit(
            Data.baseName, Data.queryFile, Data.aln, "duplication"
        )
    else:
        logger.error("Alignment and/or gene tree file have not been provided.")
        sys.exit()

    return Data, dico


def pspEntry(Data, parameters):
    Data, dico = duplPSEntry(Data, "positiveSelection")

    Data.alnFormat = parameters["alnformat"].title()

    return Data, dico


##==============================================================================================================================
