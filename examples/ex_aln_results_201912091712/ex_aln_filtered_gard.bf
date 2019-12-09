inputRedirect = {};
inputRedirect["01"] = "ex_aln_results_201912091712/ex_aln_filtered.fasta";
inputRedirect["02"] = "010010";
inputRedirect["03"] = "None";
inputRedirect["04"] = "ex_aln_results_201912091712/ex_aln_filtered.gard";
ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GARD.bf", inputRedirect);
