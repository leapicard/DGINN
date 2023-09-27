"""
Script running the tree analysis step.
"""

import Steps
import LoadFileFunc
import AnalysisFunc

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
    # steps.Data = LoadFileFunc.phymlRecEntry(steps.Data)
    # LoadFileFunc.spTreeCheck(steps.Data, "tree", snakemake.config["duplication"])
    tree_parameters = steps.get_params(["phymlOpt"])
    steps.dAlTree = AnalysisFunc.phyMLTree(steps.Data, **tree_parameters)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])

    # Serialize data to disk
    steps.serialize_dAlTree(snakemake.output[1])
