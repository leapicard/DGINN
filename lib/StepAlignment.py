"""
Script running the alignment analysis step.
"""

import Init, os
import AnalysisFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["output"] = snakemake.output
    config["input"] = str(snakemake.input)
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["step"] = snakemake.rule
    
    parameters = Init.paramDef(config)

    # Run step
    outMafft = AnalysisFunc.runMafft(parameters)
    
    fasCov, nbOut = AnalysisFunc.covAln(outMafft, parameters)
    
    outPrank = AnalysisFunc.runPrank(fasCov, parameters)

    if os.path.exists(outPrank):
      outAli = outPrank
    else:
      print("Prank did not run on file " + fasCov)
      outAli = fasCov

    outIso = AnalysisFunc.isoformAln(outAli, parameters)
    
    os.rename(outIso, str(snakemake.output))

