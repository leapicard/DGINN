"""
Script running the blast analysis step.
"""

import Steps
import BlastFunc
import pickle

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]
    steps = Steps.Steps(
        step=snakemake.rule,
        data_file=None,
        query_file=None,
        log_file=snakemake.log[0],
        config=snakemake.config,
    )

    # Run treatBlast function
    blast_params = steps.get_params(
        ["evalue", "percID", "mincov", "APIKey", "remote", "entryQuery"]
    )
    steps.Data = BlastFunc.treatBlast(steps.Data, **blast_params)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
