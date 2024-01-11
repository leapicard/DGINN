"""
Script running the recombination analysis step.
"""

import logging
import os
import subprocess

import AnalysisFunc
import Init
import yaml
from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.recombination")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)
    config["input"] = list(snakemake.input)

    cq = config["allquery"][config["queryName"]]
    if len(config["input"]) > 1 and cq != "void":
        config["input"] = cq
    else:
        config["input"] = str(snakemake.input)
    config["step"] = snakemake.rule

    parameters = Init.paramDef(config)

    # Run step (inactivated for now)

    gardRes = config["outdir"] + "/" + config["queryName"] + "_bp.txt"
    frec = open(gardRes, "w")
    frec.write("breakpoints\n")
#    frec.write("[500]\n")
    frec.write("[]\n")
    frec.close()

    ## lQuer is the list of new queryNames
    ## lAln is the list of new alignments
    [lQuer, lAln] = AnalysisFunc.parseGard(gardRes, parameters)

    ## rerun snakemake if needed
    if len(lQuer) > 1:  # several sub
        dpar = {
            k: v
            for k, v in config.items()
            if k not in ["output", "queryName", "infile"]
        }
        dpar["queryName"] = lQuer
        dpar["infile"] = lAln
        dpar["recombination"] = False
        dpar["step"] = "alignment"
        newconfig = "." + config["queryName"] + "_config_rec.yaml"
        with open(newconfig, "w") as lout:
            yaml.dump(dpar, lout)

        subprocess.run(
            [
                "snakemake",
                "--cores=%d" % (max(1, len(lQuer))),
                "--nolock",
                "--configfile=" + newconfig,
                "--until=alignment",
            ]
        )

        os.remove(newconfig)
        lAln = [config["outdir"] + "/" + quer + "_align.fasta" for quer in lQuer]

    ## register resulting [queryName, alignment]s in output
    f = open(config["output"], "w")
    for i in range(len(lQuer)):
        f.write(lQuer[i] + "\t" + lAln[i] + "\n")

    f.close()
