#################### MODULES ####################
## Python modules
import argparse, sys, os, numpy, pandas
from Bio import SeqIO, AlignIO
from collections import defaultdict, OrderedDict, Counter
from time import localtime, strftime
from statistics import median
from scipy import stats
import glob

#################################################
## Functions

def dict2fasta(dico):
    txtoutput = ""
    for key, value in dico.items():
        txtoutput += ">{:s}\n{:s}\n".format(str(key),str(value))
    return(txtoutput)

def check(arg):	
    if not os.path.exists(arg):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: arg does not exist, please enter valid argument.
        raise argparse.ArgumentTypeError("\"{0}\" does not exist, please enter valid argument.".format(arg))
    
    return os.path.abspath(arg)

def LRT(ll1, ll2, df):
    """
    Calculates likelihood ratio test between two models.
    :params ll1, ll2: likelihood of the two models studied
    :param df: degrees of freedom of difference between the two models
    """
    stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
    LR = max(2*(ll2-ll1),0)
    p = stats.chisqprob(LR, df)
    return(LR, p)

def ResBusted(baseName, posDir):
    WG = glob.glob(posDir+"/busted/busted.out")
    posSel = ""
    d = OrderedDict({"BUSTED_sel":"na","BUSTED_pv":"na"})

    if len(WG)>= 1 and os.path.exists(WG[0]):
        with open(WG[0], "r") as wg:
            try:
                p = wg.read().split("**")[-2].replace(" ", "").split("=")[1]

                if float(p) < 0.05:
                    posSel = "Y"
                else:
                    posSel = "N"
                d["BUSTED_sel"] = posSel
                d["BUSTED_pv"] = str(p)
            except:
                d = {}

    return(d)

def ResMeme(baseName, posDir, pv):
    BS = glob.glob(posDir+"/meme/meme.out")
    d = OrderedDict({"MEME_NbSites":"0", "MEME_PSS":"na"})

    if len(BS)>=1 and os.path.exists(BS[0]):
        try:
            bs = open(BS[0], "r")
            dPS={}
            for BSLine in bs.readlines():
                          if not BSLine.startswith("|"):
                            continue
                          Bspl= list(map(lambda s:s.strip(),BSLine.split("|")))
                          if not Bspl[1].isdigit():
                            continue
                          dPS[int(Bspl[1])]=float(Bspl[7].split("=")[1])
                        
            dPSok={k:v for k,v in dPS.items() if v<=pv}
            nbSites = len(dPSok)

            if nbSites > 0:
                d["MEME_NbSites"] = str(nbSites)
                d["MEME_PSS"] = ",".join(map(str,dPSok.keys()))
            bs.close()
        except:
            next

    """
    if len(d) == 0:
        d["MEME_PSS"] = "na"
        d["MEME_NbSites"] = 0
    """
    return(d)

def ResBppExtract(models, dLogLlh, dSAres, posDir, baseName, pr):
    [model1, model2] = models.split(" ")
    method = "Bpp{:s}{:s}".format(model1, model2)
    d = OrderedDict({"ps":"na", "pv":"na", "NbSites":"na", "wPS":"na", "PSS":"na"})
    if isinstance(dLogLlh[model1], float) and isinstance(dLogLlh[model2], float):
        try:
            if model1  in dLogLlh and model2 in dLogLlh:
                LR, p = LRT(dLogLlh[model1], dLogLlh[model2], 2)
                #if p < 0.05 and os.path.exists(dSAres[model2][0]):
                if len(dSAres[model2]) and os.path.exists(dSAres[model2][0]):
                    df = pandas.read_csv(dSAres[model2][0], sep='\t')
                    val = [float(x.split("=")[1]) if x[:2]=="Pr" else 0 for x in df.columns ]
                    ival = [i for i in range(len(val))  if val[i]>1 ]
                    if len(ival)>0:
                      wPS = sum([val[i] for i in ival])/len(ival)
                    else:
                                          wPS=0
                    lRespp= df[df.iloc[:,ival].sum(axis=1)>pr].iloc[:,0].tolist()
                    lResw = df[df.iloc[:,-1]>1].iloc[:,0].tolist()
                    #if wPS<5:
                    lRes = [x for x in lResw if x in lRespp]
                    #else:
                    #  lRes = lResw
                    #lRes = list(map(lambda x:x+1,lRes))
                    if len(lRes)!=0:
                            lResFinal = str("{}".format(lRes).replace("[", "").replace("]", ""))
                    else:
                            lResFinal = "na"
                    posSel = "NY"[p<0.05]
                    nbPSS = str(len(lRes))
                    wPS=str(wPS)
                else:
                    posSel = "N"
                    lResFinal = "na"
                    nbPSS = "0"
                    wPS = "0"
                        
                d["ps"] = posSel
                d["pv"] = str(p)
                d["NbSites"] = nbPSS
                d["PSS"] = lResFinal
                d["wPS"] = wPS
                    
        except:
            next
    return(d)

def ResBpp(baseName, posDir, pr):
    dSA = {}
    dSAres = {}
    dLogLlh = {}
    lModels = ["M1","M2","M7", "M8a", "M8", "M10", "DFP07_0", "DFP07"]
    lcpl=[("M1","M2"),("M7","M8"), ("M7","M10"), ("M8a","M8"), ("DFP07_0","DFP07")]
    res={}
    for suff in "GC":
          dSA[suff]={}
          dSAres[suff]={}
          dLogLlh[suff]={}
          for model in lModels:
              dSA[suff][model] = glob.glob(posDir+"/bpp_site/*"+model+"_"+suff+".params")
              dSAres[suff][model] = glob.glob(posDir+"/bpp_site/*_results_"+model+"_"+suff+".log")
              if len(dSA[suff][model])>=1 and os.path.exists(dSA[suff][model][0]):
                  with open(dSA[suff][model][0], "r") as params:
                      dLogLlh[suff][model] = float(params.readline().strip().split("= ")[-1])

              else:
                  dLogLlh[suff][model] = "na"
          res[suff]={}
          for (m1,m2) in lcpl:
            res[suff][(m1,m2)] = ResBppExtract(m1+" "+m2, dLogLlh[suff], dSAres[suff], posDir, baseName, pr)
    return(res)


def ResPamlExtract(models, dModelLlh, dModelFile, pr):
    [model1, model2] = models.split(" ")
    method = "codeml{:s}:{:s}".format(model1, model2)
    d = OrderedDict({method:"na", method+"_pv":"na", method+"_NbSites":"na", method+"_PSS":"na"})

    if isinstance(dModelLlh[model1], float) and isinstance(dModelLlh[model2], float):
        LR, p = LRT(dModelLlh[model1], dModelLlh[model2], 2)

        if p < 0.05:
            with open(dModelFile[model2], "r") as modFile:
                content = modFile.read()
                res = content.split("BEB")[1].split("The grid")[0].split("SE for w")[-1].split("\n")
                res = list(filter(None, res))
                PSS = [int(line.strip().split(" ")[0]) for line in res if float(list(filter(None, line.split(" ")))[2].replace("*", "")) > pr]
            
            posSel = "Y"
            nbPSS = str(len(PSS))
            lResFinal = str("{}".format(PSS).replace("[", "").replace("]", ""))
        
        else:
            posSel = "N"
            lResFinal = "na"
            nbPSS = "0"

        d[method] = posSel
        d[method+"_pv"] = str(p)
        d[method+"_NbSites"] = nbPSS
        d[method+"_PSS"] = lResFinal
    
    return(d)

def ResPaml(posDir, pr):
    PAML = posDir+"/paml_site/C/"
    lModels = ["M1","M2","M7","M8","M8a"]
    lcpl=[("M1","M2"),("M7","M8"), ("M8a","M8")]
    dModelLlh = OrderedDict({model:"na" for model in lModels})
    dModelFile = {model:"na" for model in lModels}
    #dModelOmega = {model:"na" for model in lModels}

    if os.path.exists(PAML):
        pamlDirs = {dI:os.path.join(PAML, dI) for dI in sorted(os.listdir(PAML)) if os.path.isdir(os.path.join(PAML, dI))}
        for pDir in pamlDirs:
            for model in lModels:
                if model in pamlDirs[pDir].split("/"):
                    dModelFile[model] = pamlDirs[pDir]+"/out"
                    if os.path.exists(pamlDirs[pDir]+"/rst1"):
                        with open(pamlDirs[pDir]+"/rst1", "r") as modelFile:
                            line = modelFile.read().strip()
                            try:
                                dModelLlh[model] = float(line.split("\t")[-1])
                                #dModelOmega[model] = float(line.split("\t")[-4])
                            except ValueError:
                                dModelLlh[model] = "na"
        
        for model in lModels:
            if model not in dModelLlh.keys():
                dModelLlh[model] = "na"

    res={}
    for (m1,m2) in lcpl:
          res[(m1,m2)] = ResPamlExtract(m1+" "+m2, dModelLlh, dModelFile, pr)
    return(res)
    
def getCov(fAln):
    lCov = []
    try:
            aln = AlignIO.read(open(fAln, "r"), "fasta")
    except ValueError:
                print("Alignment does not fit "  + fAln)
                return
    nbSeq = len(aln)
    alnLen = aln.get_alignment_length()
    for i in range(0, alnLen):
        nb = [record.seq[i] for record in aln].count("-")
        lCov.append((nbSeq-nb)/nbSeq)

    return(lCov[1::3])

