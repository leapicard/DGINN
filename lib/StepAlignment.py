"""
Script running the alignment analysis step.
"""

import logging
import os
import Init
import AnalysisFunc
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.alignment")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
    config["step"] = snakemake.rule

    aligner = config.get("aligner","macse")

    config["input"] = os.path.join(config["outdir"],config["queryName"]+"_orf.fasta")

    parameters = Init.paramDef(config)

    # Run step
    if config["step"]=="tree": # alignment already done
        os.symlink(config["input"],config["output"])
        sys.exit(0)


    ## Fix isoforms within each species through species specific MSA


    ## Global Mafft

    ORFs=parameters["input"]

    ## cluster isoforms based on mafft alignment
    outIso = AnalysisFunc.isoformMafft(ORFs, parameters)

    outMafft = AnalysisFunc.runMafft(outIso, parameters)

    ## discard sequences according to coverage on query sequence
    fasCov, nbOut = AnalysisFunc.covAln(outMafft, parameters)

    if aligner == "prank":
      outFile = AnalysisFunc.runPrank(fasCov, parameters)
    elif aligner == "macse":
      outFile = AnalysisFunc.runMacse(fasCov, parameters)
    else:
      print("Unknown aligner: " + aligner)
    
      
    if os.path.exists(outFile):
      outAli = outFile
    else:
      print(aligner + " did not run on file " + fasCov)
      outAli = fasCov

    ## Rerun isoform clustering
    outIso2 = AnalysisFunc.isoformAln(outAli, parameters)

    os.rename(outIso2, config["output"])
