"""
Script running the alignment analysis step.
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
    steps.Data = LoadFileFunc.prankEntry(steps.Data)
    LoadFileFunc.spTreeCheck(steps.Data, "alignment", snakemake.config["duplication"])
    AnalysisFunc.alnMafft(steps.Data)
    fasCov, nbOut = AnalysisFunc.covAln(
        steps.Data.aln, snakemake.config["mincov"], steps.Data.queryName, steps.Data.o
    )
    steps.Data.aln = AnalysisFunc.runPrank(fasCov, steps.Data.o)
    steps.Data.aln = AnalysisFunc.isoformAln(steps.Data.aln, steps.Data.o)

    # Serialize data to disk
    steps.serialize_data(snakemake.output[0])
