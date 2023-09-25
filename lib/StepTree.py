"""
Script running the tree analysis step.
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
        debug=snakemake.config["debug"],
    )

    # Run step
    steps.Data = LoadFileFunc.phymlRecEntry(steps.Data)
    LoadFileFunc.spTreeCheck(steps.Data, "tree", snakemake.config["duplication"])
    tree_parameters = steps.get_params(["phymlOpt"])
    AnalysisFunc.phyMLTree(steps.Data, **tree_parameters)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
