"""
Script running the alignment analysis step.
"""

import Init, os
import AnalysisFunc

if __name__ == "__main__":
    # Init and run analysis steps
    snakemake = globals()["snakemake"]

    config = snakemake.config
    config["queryName"] = str(snakemake.wildcards).split(":",1)[0]
    config["output"] = str(snakemake.output)
    cq = config["allquery"][config["queryName"]]
    if len(snakemake.input)>1 and cq!="void":
      config["input"] =  cq
    else:
      config["input"] = str(snakemake.input)
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
    
    os.rename(outIso, config["output"])

