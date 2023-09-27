"""
Script running the recombination analysis step.
"""

import Steps
import AnalysisFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["outputs"] = snakemake.output
    steps = Steps.Steps(
        step=snakemake.rule,
        data_file=snakemake.input[0],
        dAlTree_file=snakemake.input[1],
        log_file=snakemake.log[0],
        config=snakemake.config,
    )

    # Run step
    if steps.parameters["recombination"]:
        recomb_params = steps.get_params(["hyphySeuil", "hostfile"])
        steps.dAlTree = AnalysisFunc.gardRecomb(
            steps.Data, steps.dAlTree, **recomb_params
        )

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
    # Serialize dAllTree to disk
    steps.serialize_dAlTree(snakemake.output[1])
