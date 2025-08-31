# coding: utf8
import sys
import FastaResFunc, shutil
import logging, subprocess, shlex, os, ete3
from Bio import SeqIO, Align, AlignIO
from collections import defaultdict
from itertools import chain
import pandas as pd
import re

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
        # only keep sequences with limited number of Ns
        if fseq.count('N') < 10 and fseq.count('n') < 10:
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
        species = set(["_".join(x.split("_")[:2]) for x in dupl])

        for sp in species:
            if queryName in dupl:
                firstOcc = queryName
            else:
                lOcc = [x for x in dupl if sp in x]

                if len(lOcc) > 0:
                    firstOcc = lOcc[0]
                else:
                    firstOcc = str(lOcc)

            if firstOcc:
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


def runMafft(ORFs, parameters, verbose=True):
    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    logger = logging.getLogger("main.alignment")
    if verbose:
      logger.info("Started Mafft nucleotide alignment")
    
    cmdmafft = "mafft --auto --quiet {}".format(ORFs)
    outMafft = outdir+"/"+ queryName + "_mafft.fasta"
    cmd(cmdmafft,False,stdout=outMafft)
    # with open(outMafft, "w") as outM:
    #     run = subprocess.run(
    #         lCmd, shell=False, check=True, stdout=outM, stderr=subprocess.PIPE
    #     )

    if verbose:
      logger.info("Finished Mafft nucleotide alignment: {:s}".format(outMafft))

    return outMafft


def covAln(aln, parameters):
    """
    Function to discard sequences from alignment according to coverage on sequence query.

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

def isoformMafft(ORFs, parameters):
    """
    Cluster isoforms of the same species using MAFFT and return the resulting sequences file.
    """

    import os
    import logging
    from Bio import SeqIO

    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    logger = logging.getLogger("main.alignment")
    logger.info("Clustering isoforms.")

    dRem = {}       # for sequences that are not clustered
    dId2Seq = {}    # species -> {tag -> sequence}
    laln = 0        # alignment length

    # Step 1: Read sequences and organize by species and tag
    for fasta in SeqIO.parse(open(ORFs), "fasta"):
        parts = fasta.id.split("_")
        if len(parts) > 2:  # normal format: Species_species_tag...
            sp = "_".join(parts[:2])
            tag = "_".join(parts[2:])
            if sp not in dId2Seq:
                dId2Seq[sp] = {}
            dId2Seq[sp][tag] = str(fasta.seq)
            if laln == 0:
                laln = len(fasta.seq)
        else:
            dRem[fasta.id] = str(fasta.seq)

    outCov = os.path.join(outdir, queryName + "_clustiso.fasta")
    clustok = False

    # Step 2: Process each species
    for sp, dtagseq in dId2Seq.items():
        if len(dtagseq) == 1:
            # Only one sequence: keep as is
            single_tag = list(dtagseq.keys())[0]
            dRem[sp + "_" + single_tag] = dtagseq[single_tag]
            continue

        # --- Original tags (before MAFFT) ---
        ltags = list(dtagseq.keys())

        # Map tags to temporary IDs for MAFFT
        tag2tmp = {tag: f"seq{i}" for i, tag in enumerate(ltags)}
        tmp2tag = {v: k for k, v in tag2tmp.items()}

        tmpfseq = os.path.join(outdir, sp + ".fasta")
        with open(tmpfseq, "w") as outC:
            for tag, seq in dtagseq.items():
                outC.write(f">{tag2tmp[tag]}\n{seq}\n")

        # Step 3: Run MAFFT
        parameters["queryName"] = sp + "_" + queryName
        outMafft = runMafft(tmpfseq, parameters, False)
        parameters["queryName"] = queryName

        # --- Map aligned sequences back to original tags ---
        dtagseq_aligned = {}
        laln = 0
        for fasta in SeqIO.parse(open(outMafft), "fasta"):
            original_tag = tmp2tag[fasta.id]
            dtagseq_aligned[original_tag] = str(fasta.seq)
            if laln == 0:
                laln = len(fasta.seq)
        dtagseq = dtagseq_aligned

        # --- FIX: Refresh ltags so they match new dtagseq keys ---
        ltags = list(dtagseq.keys())

        # --- Safe cleanup of temporary files ---
        for f in (tmpfseq, outMafft):
            if f and os.path.exists(f):
                os.remove(f)

        # Step 4: Cluster sequences
        lclust = [[ltags[0]]]  # initialize clusters
        for itag in range(1, len(ltags)):
            tag = ltags[itag]
            seq1 = dtagseq[tag]
            for clust in lclust:
                for tag2 in clust:
                    seq2 = dtagseq[tag2]
                    dist = len([
                        pos for pos in range(laln)
                        if seq1[pos] != seq2[pos] and seq1[pos] != "-" and seq2[pos] != "-"
                    ])
                    if dist != 0:  # sequences differ
                        break
                if dist == 0:  # sequence matches cluster
                    clust.append(tag)
                    tag = ""
                    break
            if tag != "":
                lclust.append([tag])

        # Step 5: Merge sequences in each cluster
        for clust in lclust:
            if len(clust) == 1:
                dRem[sp + "_" + clust[0]] = dtagseq[clust[0]]
            else:
                clustok = True
                ntag = clust[-1] + "_clust"
                logger.info(
                    "Clustered sequences " +
                    sp + "_" +
                    (", %s_" % (sp)).join(clust) +
                    " into %s_" % (sp) + ntag
                )
                nseq = "".join(
                    [max([dtagseq[tag][pos] for tag in clust]) for pos in range(laln)]
                )
                dRem[sp + "_" + ntag] = nseq

    # Step 6: Write final sequences if clustering occurred
    if clustok:
        with open(outCov, "w") as outC:
            outC.write(FastaResFunc.dict2fasta(dRem))
        logger.info("%d remaining sequences" % len(dRem))
        return outCov
    else:
        return ORFs
    
def isoformAln(aln, parameters):
    """Function to cluster isoforms according to the alignment. Return the
    clustering of these isoforms.

        Isoforms are from the same species (recognized through 
        Species_species at the beginning of their name) and same letters or
        indels at same positions in alignment.

        @param1 aln: Path to alignment
        @param2 parameters: dict of parameters
        @return Path to file of resulting alignment

    """

    outdir = parameters["outdir"]
    queryName = parameters["queryName"]
    logger = logging.getLogger("main.alignment")
    logger.info("Clustering isoforms.")

    dRem = {}  # for remaining sequences
    dId2Seq = {}  # for remaining sequences
    laln = 0  # alignement length
    for fasta in SeqIO.parse(open(aln), "fasta"):
        post = len(fasta.id.split("_"))
        if post > 2:  # regular format
            sp = "_".join(fasta.id.split("_")[:2])
            tag = "_".join(fasta.id.split("_")[2:])
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
        ltags=list(dtagseq.keys())
        lclust = [[ltags[0]]]  # list of clusters of tags
        for itag in range(1,len(ltags)):
          tag=ltags[itag]
          seq1=dtagseq[tag]
          for clust in lclust:
            for tag2 in clust:
              seq2=dtagseq[tag2]
              dist=len([pos for pos in range(laln) if seq1[pos]!=seq2[pos] and seq1[pos]!="-" and seq2[pos]!="-"])
              if dist!=0: # not in the cluster
                break
            if dist==0: # final success -> in clust
              clust.append(tag)
              tag=""
              break
          if tag!="": # not put in a cluster
            lclust.append([tag])

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
        logger.info("%d remaining sequences"%(len(dRem)))
        return outCov
    else:
        return aln

      

#######=================================================================================================================
######PhyML=============================================================================================================

def fasta2phylip(aln, outPhy, format):
    tmp = aln+".tmp"
    with open(aln, "r") as aln2:
        laln = aln2.read().replace("!", "N")
        aln2.close()
    with open(tmp, "w") as temp:
      temp.write(laln)
      temp.close()

    input_handle = open(tmp, "r")
    output_handle = open(outPhy, "w")

    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, format)

    output_handle.close()
    input_handle.close()
    os.remove(tmp)
  

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
    outPhy = geneDir + "/" + queryName + ".phy"
    # aln = aln.split("/")[-1]
    tmp = geneDir + "/" + queryName + "tree.tmp"

    fasta2phylip(aln, outPhy, format="phylip-relaxed")
    
    logger = logging.getLogger("main.tree")
    logger.info("Run PhyML builder.")

    phymlOpt = parameters["phymlOpt"]
    # PhyML
    if phymlOpt != "":
        try:
            opt = phymlOpt.split("ALN ")[1]
            logger.debug("phyml --quiet -i {:s} {}".format(outPhy, opt))
            cmd("phyml --quiet -i {:s} {}".format(outPhy, opt), False)
        except:
            logger.info(
                "PhyML couldn't run with the provided info {}, running with default options.".format(
                    phymlOpt
                )
            )
            cmd("phyml --quiet -i {:s} -v e -b -2".format(outPhy), False)
    else:
        logger.debug("phyml --quiet -i {:s} -v e -b -2".format(outPhy))
        cmd("phyml --quiet -i {:s} -v e -b -2".format(outPhy), False, False)

    return outPhy+"_phyml_tree.txt"

#######=================================================================================================================
###### IqTree =============================================================================================================


def runIqTree(parameters):
    geneDir = parameters["outdir"]
    aln = parameters["input"]
    queryName = parameters["queryName"]

    logger = logging.getLogger("main.tree")
    logger.info("Run IqTree builder.")
    logger.debug("iqtree2 -redo --quiet -s {:s}".format(aln))
    cmd("iqtree2 -redo --quiet -s {:s}".format(aln), False)

    return aln+".treefile"


#######=================================================================================================================

######GARD==============================================================================================================

def runPhymlMulti(parameters):
    geneDir = parameters["outdir"]
    aln = parameters["input"]
    queryName = parameters["queryName"]
    outPhy = geneDir + "/" + queryName + ".phy"
    tmp = geneDir + "/" + queryName + "tree.tmp"

    logger = logging.getLogger("main.recombination")
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

    #phymlOpt = parameters["phymlOpt"]
    # PhyML
    if False:#phymlOpt != "":
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
    else:
        logger.info("phyml_multi %s 0 i 1 0 HKY e e 4 e BIONJ y y y 2"%(outPhy))
        subprocess.run("phyml_multi %s 0 i 1 0 HKY e e 4 e BIONJ y y y 2"%(outPhy),shell=True)

    ################
    ### SARMENT HMM

    outf = outPhy + "_phyml_siteLks.txt"
    fs = open(outPhy + "_phyml_lk.txt","r")
    autocor = 0.99
    for l in fs.readlines():
      if l.startswith("Autocorrelation"):
        p=l.find(":")
        autocor = float(l[p+1:])
        break
    fs.close()
  
    fout=open(os.path.join(outdir,"outpart.txt"),"w")
    subprocess.run("python3 /usr/bin/SARMENT/PartitioningHMM.py %s %f"%(outf,autocor),shell=True,stdout=fout)
    fout.close()

    l = fout.readline()
    while l:
      if l.find("Forward")!=-1:
        l = fout.readline()
        break
      l = fout.readline()

    pr=re.compile(r"<(\d+)-\d+>")
    lbeg=list(map(int,pr.findall(l)))

    threshold = 20
    ## avoid too short segments

    d=0
    lbeg2=[]
    for x in lbeg:
      if x-d>=threshold:
        lbeg2.append(x)
        d=x
    return(lbeg2)
        
  
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



