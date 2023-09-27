"""
Script running the positive selection analysis step.
"""

import Steps
import PosSelFunc

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
    if steps.parameters["positiveSelection"]:
        listArgsPosSel = []
        fAT = open(steps.Data.o + "files_list.txt", "w")
        for aln in steps.dAlTree:
            listArgs = [steps.Data, steps.parameters, aln, steps.dAlTree[aln]]
            listArgsPosSel.append(listArgs)
            if len(steps.dAlTree[aln]):  # only if it exists
                fAT.write(aln + "\t" + steps.dAlTree[aln])
                outDir = PosSelFunc.pspAnalysis(
                    steps.Data, steps.parameters, aln, steps.dAlTree[aln]
                )
                fAT.write("\t" + outDir + "\n")

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
    # Serialize dAllTree to disk
    steps.serialize_dAlTree(snakemake.output[1])
