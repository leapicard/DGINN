# coding: utf8
import sys
import FastaResFunc, shutil
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, AlignIO
from collections import defaultdict
from itertools import chain
import pandas as pd

######Functions=============================================================================================================


def cmd(commandLine, choice, verbose=False, stdout = None):
    """
    Function executing a command line in a bash terminal.

    @param1 commandLine: String corresponding to a bash command line
    @param2 choice: Boolean determining whether the command is executed within the shell
    """
    if not stdout:
      if verbose:
        out = None
      else:
        out = subprocess.PIPE
    else:
      if stdout == subprocess.PIPE:
        out = stdout
      else:
        out = open(stdout, "w")

    if not choice:  ## if shell==False, split command & arguments. 
      commandLine = shlex.split(commandLine)

    try:
        run = subprocess.run(commandLine, shell=choice, stdout=out, stderr=out)
    except subprocess.CalledProcessError as err:
      sys.stderr.write(str(err))

        
    if stdout and stdout != subprocess.PIPE:
      out.close()
    return(run.stdout)
  
######ORF===================================================================================================================

def getORFs(parameters):
    """
    Function to find Open Reading Frames within the sequence of each gene and select the longest one.

    @param1 catFile: Path
    @param2 geneName: Gene name
    @param3 geneDir: Gene directory
    @return outORF: Path to the file containing the longest ORFs
    """

    catFile = parameters["input"]
    queryName = parameters["queryName"]
    geneDir = parameters["outdir"]
    
    outORFraw = geneDir + "/"+queryName+"_allORFs.fasta"
    logger = logging.getLogger("main.orf")

    logger.debug(
        "getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(
            catFile, outORFraw
        )
    )
    cmd(
        "getorf -sequence {:s} -outseq {:s} -table 0 -find 3 -noreverse".format(
            catFile, outORFraw
        ),
        False,
    )

    dId2ORFs = defaultdict(list)
    f = SeqIO.parse(open(outORFraw), "fasta")
    for fasta in f:
        fname, fseq = fasta.id, str(fasta.seq)
        if len(fname.split("_")) > 2:
            fname2 = "_".join(fname.split("_")[0:-1])
        else:
            fname2 = fname.split("_")[0]
        dId2ORFs[fname2].append(fseq)

    dId2Longest = {}
    for k, v in dId2ORFs.items():
        dId2Longest[k] = max(v, key=len)

    # delete duplicate sequences
    dRev = {}
    for k, v in dId2Longest.items():
        dRev.setdefault(v, set()).add(k)

    AllDupl = [values for key, values in dRev.items() if len(values) > 1]
    n = 0
    for dupl in AllDupl:
        species = set([x.split("_")[0] for x in dupl])

        for sp in species:
            if queryName in dupl:
                firstOcc = queryName
            else:
                lOcc = [x for x in dupl if sp in x]

                if len(lOcc) > 0:
                    firstOcc = lOcc[0]
                else:
                    firstOcc = str(lOcc)

            dupl.remove(firstOcc)

        for i in dupl:
            dId2Longest.pop(i, None)
            n += 1
            logger.debug("Deleted sequence {:s} (duplicate)".format(i))

    logger.info("Deleted {} sequences as duplicates".format(n))

    outORF = outORFraw.replace("allORFs.fasta", "orf.fasta")

    with open(outORF, "w") as outO:
        outO.write(FastaResFunc.dict2fasta(dId2Longest))
        outO.close()

    logger.info("Extracted longest ORFs: {:s}".format(outORF))

    return outORF



#######=================================================================================================================
######PRANK=============================================================================================================


def runPrank(ORFs, parameters):
    """
    Function to run PRANK software for codon alignment (Löytynoja, 2014).

    @param1 ORFs: Path
    @param2 geneDir: Gene Directory
    @return outPrank: Path to Prank results file
    """
    logger = logging.getLogger("main.alignment")
    logger.info("Started Prank codon alignment")
    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    outPrank = outdir + "/" + queryName + "_prank"

#    logger.debug("prank -d={:s} -o={:s} -codon -F".format(ORFs, outPrank))

    cmd("prank -d={:s} -o={:s} -codon -F".format(ORFs, outPrank), False)

    logger.info("Finished Prank codon alignment: {:s}.best.fas".format(outPrank))

    if os.path.exists(outPrank + ".fas"):
      shutil.copyfile(outPrank + ".fas", outPrank + ".best.fas")

    return outPrank + ".best.fas"

  
def runMacse(ORFs, parameters):
    """
    Function to run MACSE software for codon alignment.

    @return Path to Macse results file
    """
    logger = logging.getLogger("main.alignment")
    logger.info("Started Macse codon alignment")
    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    outFile = outdir + "/" + queryName + "_macse"

    fout = open(outFile+".out","w")

    subprocess.run("java -jar /opt/bin/macse.jar -prog refineAlignment -align {:s} -out_NT {:s}.best.fas".format(ORFs, outFile),shell=True, stdout=fout)

    fout.close()
    
    logger.info("Finished Macse codon alignment: {:s}.best.fas".format(outFile))

    return outFile + ".best.fas"


def runMafft(parameters):
    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    logger = logging.getLogger("main.alignment")
    logger.info("Started Mafft nucleotide alignment")

    
    ORFs=parameters["input"]
    cmdmafft = "mafft --auto --quiet {}".format(ORFs)
    outMafft = outdir+"/"+ queryName + "_mafft.fasta"
    cmd(cmdmafft,False,stdout=outMafft)
    # with open(outMafft, "w") as outM:
    #     run = subprocess.run(
    #         lCmd, shell=False, check=True, stdout=outM, stderr=subprocess.PIPE
    #     )

    logger.info("Finished Mafft nucleotide alignment: {:s}".format(outMafft))

    return outMafft


def covAln(aln, parameters):
    """
    Function to discard sequences from alignment according to coverage to query.

    @param1 aln: Path to next alignment
    @param2 parameters: 
    @return outCov: Path to file of conserved sequences
    """

    cov = parameters["mincov"]
    queryName = parameters["queryName"]

    dId2Seq = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(open(aln), "fasta")}
    logger = logging.getLogger("main.alignment")

    if queryName in dId2Seq:
        logger.info(
            "Discarding sequences with less than {:d}% coverage of query.".format(cov)
        )
        outCov = outdir+"/"+ queryName + "_mincov.fasta"
        
        nbOut = 0
        lIndexes = [pos for pos, char in enumerate(dId2Seq[queryName]) if char != "-"]

        dKeep = {}
        for ID, seq in dId2Seq.items():
            seqPos = [seq[x] for x in lIndexes]
            seqCov = (len(seqPos) - seqPos.count("-")) / len(seqPos) * 100

            if seqCov > cov:
                dKeep[ID] = seq

        nbOut = len(dId2Seq) - len(dKeep)

        with open(outCov, "w") as outC:
            outC.write(FastaResFunc.dict2fasta(dKeep))
            outC.close()
        logger.info("Discarded {:d} sequences".format(nbOut))

        return (outCov, nbOut)

    else:
        logger.warning(
            "Provided query name not found in the alignment, skipping coverage check."
        )
        return (aln, 0)


def isoformAln(aln, parameters):
    """Function to cluster isoforms according to the alignment. Return the
    overall coverage of these isoforms.

        Isoforms are from the same species (recognized through keyword
        xxxXxx at the beginning of their name) and same letters or
        indels at same positions in alignment.

        @param1 aln: Path to alignment
        @param2 o: output directory
        @return outAln: Path to file of resulting alignment

    """

    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    logger = logging.getLogger("main.alignment")
    logger.info("Clustering isoforms.")

    dRem = {}  # for remaining sequences
    dId2Seq = {}  # for remaining sequences
    laln = 0  # alignement length
    for fasta in SeqIO.parse(open(aln), "fasta"):
        post = fasta.id.find("_")
        if post != -1:  # regular format
            sp = fasta.id[:post]
            tag = fasta.id[post + 1 :]
            if not sp in dId2Seq:
                dId2Seq[sp] = {}
            dId2Seq[sp][tag] = str(fasta.seq)
            if laln == 0:
                laln = len(fasta.seq)
        else:
            dRem[fasta.id] = str(fasta.seq)

    outCov = outdir + "/"  + queryName + "_clustiso.fasta"
    clustok = False  # flag to check if a cluster has occured
    for sp, dtagseq in dId2Seq.items():
        lclust = [list(dtagseq)]  # list of clusters of tags to be split
        for pos in range(laln):
            lclust2 = []
            for clust in lclust:
                dlet = {tag: dtagseq[tag][pos] for tag in clust}
                llet = set([x for x in dlet.values() if x != "-"])
                if len(llet) <= 1:  # one letter at most, keep all
                    lclust2.append(clust)
                    continue
                else:
                    for x in llet:
                        lclust2.append([tag for tag in clust if dlet[tag] == x])
                    lind = [
                        tag for tag in clust if dlet[tag] == "-"
                    ]  # conservative, do not know wether to merge, may be improved
                    if len(lind) != 0:
                        lclust2.append(lind)
            lclust = lclust2

        # now merge sequences in each cluster
        for clust in lclust:
            if len(clust) == 1:
                dRem[sp + "_" + clust[0]] = dtagseq[clust[0]]
            else:
                clustok = True
                ntag = clust[-1] + "_clust"
                logger.info(
                    "Clustered sequences "
                    + sp
                    + "_"
                    + (", %s_" % (sp)).join(clust)
                    + " into %s_" % (sp)
                    + ntag
                )
                nseq = "".join(
                    [max([dtagseq[tag][pos] for tag in clust]) for pos in range(laln)]
                )
                dRem[sp + "_" + ntag] = nseq

    if clustok:
        with open(outCov, "w") as outC:
            outC.write(FastaResFunc.dict2fasta(dRem))
            outC.close()

        return outCov
    else:
        return aln



#######=================================================================================================================
######PhyML=============================================================================================================


def runPhyML(parameters):
    """
    Function converting fasta file to phylip and running PhyML.

    @param1 aln: Path
    @param2 geneDir: Gene directory
    @return outPhy: Path to PhyML results file
    """
    # convert to Phylip format and replace eventual "!" symbols (relic from using MACSE)

    geneDir = parameters["outdir"]
    aln = parameters["input"]
    queryName = parameters["queryName"]
    # origin = os.getcwd()
    # os.chdir(geneDir)
    outPhy = geneDir + "/" + queryName + "_tree.phylip"
    # aln = aln.split("/")[-1]
    tmp = geneDir + "/" + queryName + "tree.tmp"

    logger = logging.getLogger("main.tree")
    with open(aln, "r") as aln2:
        laln = aln2.read().replace("!", "N")
        aln2.close()
        with open(tmp, "w") as temp:
            temp.write(laln)
            temp.close()

    input_handle = open(tmp, "r")
    output_handle = open(outPhy, "w")

    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "phylip-relaxed")

    output_handle.close()
    input_handle.close()
    os.remove(tmp)

    phymlOpt = parameters["phymlOpt"]
    # PhyML
    if phymlOpt != "":
        try:
            opt = phymlOpt.split("ALN ")[1]
            logger.debug("phyml -i {:s} {}".format(outPhy, opt))
            cmd("phyml --quiet -i {:s} {}".format(outPhy, opt), False)
        except:
            logger.info(
                "PhyML couldn't run with the provided info {}, running with default options.".format(
                    phymlOpt
                )
            )
            cmd("phyml -i {:s} -v e -b -2".format(outPhy), False)
    else:
        logger.debug("phyml -i {:s} -v e -b -2".format(outPhy))
        cmd("phyml -i {:s} -v e -b -2".format(outPhy), False, False)

    return outPhy+"_phyml_tree.txt"

#######=================================================================================================================
###### IqTree =============================================================================================================


def runIqTree(parameters):
    """
    Function converting fasta file to phylip and running PhyML.

    @param1 aln: Path
    @param2 geneDir: Gene directory
    @return outPhy: Path to PhyML results file
    """
    # convert to Phylip format and replace eventual "!" symbols (relic from using MACSE)

    geneDir = parameters["outdir"]
    aln = parameters["input"]
    queryName = parameters["queryName"]

    cmd("iqtree2 --quiet -s {:s}".format(aln), False)

    return aln+".treefile"



def splitTree(parameters, step="duplication"):
  """
  Split the gene tree in several sub-trees, following several methods.

  1- Cutlongbranches
  2- Reconciliation

  @output The list of new subalignment files
  """

  nbspecies=parameters["nbspecies"]

  aln = parameters["input"].split()[0].strip()
  tree = parameters["input"].split()[1].strip()
  
  logger = logging.getLogger(".".join(["main", step]))

  dSubAln = cutLongBranches(parameters, aln, tree, nbspecies, logger)

  lquery=[]
  
  if parameters["sptree"]!="":
    sptree = parameters["sptree"]

    for query, aln in dSubAln.items():
      logger.info("Running Treerecs for " + query)
      recTree = runTreerecs(query, aln, tree, sptree, outdir, logger)

      if recTree:
        lquery += treeParsing(query, aln, recTree, nbspecies, outdir, logger)
      else:
        lquery.append(query)
          
  else:
    lquery=list(dSubAln.keys())

  return(lquery)


#######=================================================================================================================

######GARD==============================================================================================================


def runGARD(parameters):
    """
    Function creating the batch file to run GARD (Kosakovsky Pond et al., 2006).

    @param1 aln: Path
    @param2 o: Ouput directory
    @return gardRes: Path to GARD output file
    """

    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    aln = parameters["input"]
    hostfile = parameters.get("hostfile","")
    
    gardRes = outdir + "/" + queryName + ".gard"
    gardJson = outdir + "/" + queryName + ".GARD.json"
    outGard = gardRes + "_out"
    errGard = gardRes + "_err"

    # for hyphy 2.5
    if hostfile == "" or hostfile is None:
        cmd = "mpirun -np 4 HYPHYMPI GARD --alignment {:s} --output {:s} --output-lf {:s}".format(
            aln, gardJson, gardRes
        )
    else:
        cmd = "mpirun -np 4 -hostfile {:s} HYPHYMPI GARD --alignment {:s} --output {:s} --output-lf {:s}".format(
            hostfile, aln, gardJson, gardRes
        )

    lCmd = shlex.split(cmd)
    with open(outGard, "w") as o, open(errGard, "w") as e:
        runGARD = subprocess.run(lCmd, shell=False, check=True, stdout=o, stderr=e)
        o.close()
        e.close()
    logger.debug(cmd)
    logger.info(gardJson)
    return gardJson


def procGARD(gardRes, aln):
    """
    Function creating the batch file to run GARD processor for KH test (Kosakovsky Pond et al., 2006).

    @param1 gardRes: Path to GARD output file
    @param2 aln: Path to alignment file
    @return otGardProc: Path to GARDprocessor output file
    """
    batchFile = gardRes.split(".")[0] + "_gardprocessor.bf"
    splitsFile = gardRes + "_splits"
    outGardProc = gardRes + "_outproc"
    errGardProc = gardRes + "_errproc"
    logger = logging.getLogger("main.recombination")

    with open(batchFile, "w") as bf:
        bf.write("inputRedirect = {};\n")
        bf.write('inputRedirect["01"] = "{:s}";\n'.format(aln))
        bf.write('inputRedirect["02"] = "{:s}";\n'.format(splitsFile))
        bf.write(
            'ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "GARDProcessor.bf", inputRedirect);\n'.format()
        )
    # logger.debug("Batch file: {:s}".format(batchFile))

    cmd = "hyphy {:s}".format(batchFile)
    logger.info(cmd)
    logger.info(os.system(cmd))
    logger.info("====================")
    lCmd = shlex.split(cmd)
    with open(outGardProc, "w") as o, open(errGardProc, "w") as e:
        runGARD = subprocess.run(lCmd, shell=False, check=True, stdout=o, stderr=e)
        o.close()
        e.close()
    return outGardProc


def parseGard(kh, parameters):
    """
    Function returning the cut fragments following GARD analysis and identification of significant breakpoints.

    @param1 kh: Path to GARD.json output file
    @return lOutFrag: List of Path (Fragments in fasta files)
    """
    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    aln = parameters["input"]
    
    lBP = []
    f = open(kh, "r")
    lLine = f.readline()
    while lLine:
        if lLine.find('breakpoints') != -1:
            lLine = f.readline()
            lLine = lLine[lLine.find("[") + 1 : lLine.find("]")]
            if len(lLine)==0:
              break
            lBP = list(map(int, lLine.split(",")))
            break
        lLine = f.readline()
    f.close()
    index = 0

    # If there're breakpoint(s), cut sequence in subsequences according to breakpoints
    if len(lBP) > 0:
        dFname2Fseq = {}
        for fasta in SeqIO.parse(open(aln), "fasta"):
            dFname2Fseq[fasta.id] = str(fasta.seq)

        # looking for a multiple of 3 (number of letter) (subsequence ends on or after the breakpoint)
        nbSeq = len(dFname2Fseq)
        lenSeq = len(dFname2Fseq[list(dFname2Fseq.keys())[0]])
        lPos = [0]
        deb=0
        dFrag=[]
        for bp in lBP:
            dec=False

            while (bp-deb) % 3 != 0:
                bp += 1
                dec=True
            nF={}
            for  name,seq in dFname2Fseq.items():
              nseq=seq[deb:bp]
              if len(nseq.replace("-",""))!=0: # Non empty sequence
                     nF[name]=nseq
            dFrag.append(nF)
            deb=bp + [0,-3][dec]  # dec to reput codon if broken
            lPos.append(bp)

        # Adding subsequences that start at the last breakpoint to the end 
        dFrag += [{name:seq[deb:] for  name,seq in dFname2Fseq.items()}]

        lBP = lPos + [lenSeq]
        lQuerFrag = []
        lOutFrag = []
        index = 0
        for x in range(len(lBP)-1):
            extension = "_{:d}-{:d}".format(lBP[x]+1, lBP[x+1])

            name = queryName + extension
            outFrag = outdir + "/" + name + "_orf.fasta"
            
            with open(outFrag, "w") as outF:
                outF.write(FastaResFunc.dict2fasta(dFrag[x]))
                outF.close()
            lQuerFrag.append(name)
            lOutFrag.append(outFrag)
        return [lQuerFrag, lOutFrag]
    else:
        return [[queryName],[aln]]


#######=================================================================================================================



