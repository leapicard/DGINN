"""
Script running the ORF analysis step.
"""

import Steps
import LoadFileFunc
import AnalysisFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]
    steps = Steps.Steps(
        step=snakemake.rule,
        data_file=snakemake.input[0],
        query_file=snakemake.input[1],
        log_file=snakemake.log[0],
        config=snakemake.config,
    )

    # Run step
    steps.Data = LoadFileFunc.orfEntry(steps.Data)
    LoadFileFunc.spTreeCheck(steps.Data, "orf", snakemake.config["duplication"])
    AnalysisFunc.orfFinder(steps.Data)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
