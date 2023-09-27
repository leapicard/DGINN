"""
Script running the duplication analysis step.
"""

import Steps
import LoadFileFunc
import AnalysisFunc
import TreeFunc

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
    # steps.Data, steps.dAlTree = LoadFileFunc.duplPSEntry(steps.Data)

    # dTree = steps.dAlTree.pop(steps.Data.aln)
    # LoadFileFunc.spTreeCheck(steps.Data, "duplication", snakemake.config["duplication"])
    # steps.dAlTree[steps.Data.aln] = dTree

    check_params = steps.get_params(["nbspecies", "LBopt"])
    steps.dAlTree = AnalysisFunc.checkPhyMLTree(
        steps.Data, steps.dAlTree, **check_params
    )
    tree_params = steps.get_params(["nbspecies", "phymlOpt"])
    steps.dAlTree = TreeFunc.treeTreatment(steps.Data, steps.dAlTree, **tree_params)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
    # Serialize dAllTree to disk
    steps.serialize_dAlTree(snakemake.output[1])
