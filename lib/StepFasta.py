"""
Script running the fasta analysis step.
"""

import Steps
import LoadFileFunc
import FastaResFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["outputs"] = snakemake.output
    steps = Steps.Steps(
        step=snakemake.rule,
        data_file=snakemake.input[0],
        query_file=snakemake.input[1],
        log_file=snakemake.log[0],
        config=snakemake.config,
    )

    # Run step
    steps.Data = LoadFileFunc.getSeqEntry(steps.Data)
    fasta_parameters = steps.get_params(["remote", "APIKey", "maxLen", "duplication"])
    FastaResFunc.fastaCreation(steps.Data, "fasta", **fasta_parameters)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
