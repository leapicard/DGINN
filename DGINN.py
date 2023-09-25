#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @package DGINN.py
# @author Lea Picard & Quentin Ganivet

"""
Main script called to run the pipeline.
"""

#################### MODULES ####################
import argparse, os, subprocess, logging, sys, shlex, time
from Bio import SeqIO
from pathlib import Path
from statistics import median

sys.path.insert(0, str((Path(__file__).parents[0] / "lib").resolve()))
import Init, BlastFunc, ExtractFunc, FastaResFunc, AnalysisFunc, LoadFileFunc, PosSelFunc, TreeFunc, DataObject

#################### GLOBAL #####################
version = "1.0"
VERSION_DATE = "2019/11/06"

#################### MAIN CODE ##################
if __name__ == "__main__":
    # List of steps
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
    ]

    # Initializing necessary variables for pipeline execution
    parse = Init.parameterCommandLine(version, "DGINN")
    parameters = parse.parse_args()
    debug = parameters.debug
    hostfile = parameters.hostfile
    """
    threads = int(parameters.threads)
    if threads < 2:
        threads = 2
    """

    parameters = Init.paramDef(
        parameters.params, parameters.infile, parameters.queryName, parameters.outdir
    )
    Data = Init.initLogger(parameters, debug, version)

    funcNeeded = Init.initPipeline(parameters["step"], lSteps)
    dAlTree = {}
    Data.setGenAttr(parameters["step"])

    firstStep = ""
    if parameters["step"] in ["blast", "accessions", "fasta"]:
        firstStep = "orf"
    elif parameters["step"] in ["recombination", "positiveSelection"]:
        firstStep = ""
    else:
        firstStep = parameters["step"]

    for i in range(len(lSteps)):
        if funcNeeded[i] == True:
            if lSteps[i] == "blast":
                Data.baseName = LoadFileFunc.baseNameInit(
                    Data.baseName, Data.queryFile, Data.aln
                )

                Data = BlastFunc.treatBlast(
                    Data,
                    parameters["evalue"],
                    parameters["percID"],
                    parameters["mincov"],
                    parameters["APIKey"],
                    parameters["remote"],
                    parameters["entryQuery"],
                )

            elif lSteps[i] == "accessions":
                if parameters["step"] == "accessions":
                    Data = LoadFileFunc.accnEntry(Data)

                ExtractFunc.treatAccns(Data)

            elif lSteps[i] == "fasta":
                if parameters["step"] == "fasta":
                    Data = LoadFileFunc.getSeqEntry(Data)

                FastaResFunc.fastaCreation(
                    Data,
                    parameters["remote"],
                    parameters["APIKey"],
                    parameters["maxLen"],
                    parameters["step"],
                    parameters["duplication"],
                )

            elif lSteps[i] == "orf":
                if parameters["step"] == "orf":
                    Data = LoadFileFunc.orfEntry(Data)

                    if lSteps[i] == firstStep:
                        LoadFileFunc.spTreeCheck(
                            Data, firstStep, parameters["duplication"]
                        )

                AnalysisFunc.orfFinder(Data)

            elif lSteps[i] == "alignment":
                if parameters["step"] == "alignment":
                    Data = LoadFileFunc.prankEntry(Data)

                    if lSteps[i] == firstStep:
                        LoadFileFunc.spTreeCheck(
                            Data, firstStep, parameters["duplication"]
                        )

                AnalysisFunc.alnMafft(Data)
                fasCov, nbOut = AnalysisFunc.covAln(
                    Data.aln, parameters["mincov"], Data.queryName, Data.o
                )
                Data.aln = AnalysisFunc.runPrank(fasCov, Data.o)
                Data.aln = AnalysisFunc.isoformAln(Data.aln, Data.o)

            elif lSteps[i] == "tree":
                if parameters["step"] == "tree":
                    Data = LoadFileFunc.phymlRecEntry(Data)

                    if lSteps[i] == firstStep:
                        LoadFileFunc.spTreeCheck(
                            Data, firstStep, parameters["duplication"]
                        )

                dAlTree = AnalysisFunc.phyMLTree(Data, parameters["phymlOpt"])

            elif lSteps[i] == "duplication" and parameters["duplication"]:
                if parameters["step"] == "duplication":
                    Data, dAlTree = LoadFileFunc.duplPSEntry(Data)

                    if lSteps[i] == firstStep:
                        dTree = dAlTree.pop(Data.aln)
                        LoadFileFunc.spTreeCheck(
                            Data, firstStep, parameters["duplication"]
                        )
                        dAlTree[Data.aln] = dTree

                dAlTree = AnalysisFunc.checkPhyMLTree(
                    Data, dAlTree, parameters["nbspecies"], parameters["LBopt"]
                )

                dAlTree = TreeFunc.treeTreatment(
                    Data, dAlTree, parameters["nbspecies"], parameters["phymlOpt"]
                )

            elif lSteps[i] == "recombination" and parameters["recombination"]:
                if parameters["step"] == "recombination":
                    Data = LoadFileFunc.phymlRecEntry(Data, "main.recombination")
                    dAlTree[Data.aln] = ""

                dAlTree = AnalysisFunc.gardRecomb(
                    Data, parameters["hyphySeuil"], dAlTree, hostfile
                )

            elif lSteps[i] == "positiveSelection" and parameters["positiveSelection"]:
                if parameters["step"] == "positiveSelection":
                    Data, dAlTree = LoadFileFunc.pspEntry(Data, parameters)

                listArgsPosSel = []
                fAT = open(Data.o + "files_list.txt", "w")
                for aln in dAlTree:
                    listArgs = [Data, parameters, aln, dAlTree[aln]]
                    listArgsPosSel.append(listArgs)

                    if len(dAlTree[aln]):  # only if it exists
                        fAT.write(aln + "\t" + dAlTree[aln])
                        outDir = PosSelFunc.pspAnalysis(
                            Data, parameters, aln, dAlTree[aln]
                        )
                        fAT.write("\t" + outDir + "\n")

    logger = logging.getLogger("main")
    logger.info("Finished DGINN analyses, exiting.")
