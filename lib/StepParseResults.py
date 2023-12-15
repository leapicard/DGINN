"""
Script running the positive selection analysis step.
"""

import Init
import yaml, os, subprocess

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["output"] = str(snakemake.output)
    config["input"] = str(snakemake.input)
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["step"] = snakemake.rule

    # Run step

    subprocess.run(['python',"etc/parseResults.py","-in=" + config["input"],"-o=" + config["output"]])



