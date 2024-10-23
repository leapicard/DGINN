"""
Script running the positive selection analysis step.
"""
import logging
import os, subprocess

from Logging import setup_logger

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    # Logging
    logger = logging.getLogger("main.parse_results")
    setup_logger(logger, snakemake.log[0])

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":", 1)[0]
    config["output"] = str(snakemake.output)

    config["input"] = os.path.join(config["outdir"], config["queryName"] + "_positive_selection.txt")

    # Run step

    subprocess.run(
        [
            "python3",
            "etc/parseResults.py",
            "-in=" + config["input"],
            "-o=" + config["output"],
        ]
    )
