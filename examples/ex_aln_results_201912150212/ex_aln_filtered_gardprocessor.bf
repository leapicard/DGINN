inputRedirect = {};
inputRedirect["01"] = "/home/lea/Documents/DGINN/examples/ex_aln_results_201912150212/ex_aln_filtered.fasta";
inputRedirect["02"] = "/home/lea/Documents/DGINN/examples/ex_aln_results_201912150212/ex_aln_filtered.gard_splits";
ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GARDProcessor.bf", inputRedirect);
