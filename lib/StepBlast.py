"""
Script running the blast analysis step.
"""

import Steps
import BlastFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["outputs"] = snakemake.output
    steps = Steps.Steps(
        step=snakemake.rule,
        data_file=None,
        query_file=None,
        log_file=snakemake.log[0],
        config=config,
    )

    # Run treatBlast function
    blast_params = steps.get_params(
        ["evalue", "percID", "mincov", "APIKey", "remote", "entryQuery"]
    )
    steps.Data = BlastFunc.treatBlast(steps.Data, **blast_params)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
