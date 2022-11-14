import logging, fastares_test
from Bio import SeqIO
import json, sys

def isoformAln(aln, o):
        """Function to cluster isoforms according to the alignment. Return the
overall coverage of these isoforms.

	Isoforms are from the same species (recognized through keyword
	xxxXxx at the beginning of their name) and same letters or
	indels at same positions in alignment.

	@param1 aln: Path to alignment
	@param2 o: output directory
	@return outAln: Path to file of resulting alignment

        """

        logger = logging.getLogger("main.alignment")
        logger.info("Clustering isoforms.")

        dRem={} #for remaining sequences
        dId2Seq={} #for remaining sequences
        laln=0 #alignement length
        for fasta in SeqIO.parse(open(aln),'fasta'):
                post=fasta.id.find("_")
                if post!=-1: #regular format
                        sp=fasta.id[:post]
                        tag=fasta.id[post+1:]
                        if not sp in dId2Seq:
                                dId2Seq[sp]={}
                        dId2Seq[sp][tag]=str(fasta.seq)
                        if laln==0:
                                laln=len(fasta.seq)
                else:
                        dRem[fasta.id]=str(fasta.seq)

        
        outCov = o
        clustok=False #flag to check if a cluster has occured
        for sp,dtagseq in dId2Seq.items():
                lclust=[list(dtagseq)] #list of clusters of tags to be split
                for pos in range(laln):
                        lclust2=[]
                        for clust in lclust:
                                dlet={tag:dtagseq[tag][pos] for tag in clust}
                                llet=set([x for x in dlet.values() if x!="-"])
                                if len(llet)<=1: #one letter at most, keep all
                                        lclust2.append(clust)
                                        continue
                                else:
                                        for x in llet:
                                                lclust2.append([tag for tag in clust if dlet[tag]==x])
                                        lind=[tag for tag in clust if dlet[tag]=="-"] #conservative, do not know wether to merge, may be improved
                                        if len(lind)!=0:
                                                lclust2.append(lind)
                        lclust=lclust2
                                        
                #now merge sequences in each cluster
                for clust in lclust:
                        if len(clust)==1:
                                dRem[sp+"_"+clust[0]]=dtagseq[clust[0]]
                        else:
                                clustok=True
                                ntag=clust[-1]+"_clust"
                                logger.info("Clustered sequences " + sp+"_" + (", %s_"%(sp)).join(clust) + " into %s_"%(sp)+ntag)
                                nseq="".join([max([dtagseq[tag][pos] for tag in clust]) for pos in range(laln)])
                                dRem[sp+"_"+ntag]=nseq

        if clustok:
            with open(outCov, "w") as outC:
                outC.write(fastares_test.dict2fasta(dRem))
                outC.close()
	
                return(outCov)
        else:
                return(aln)

if __name__ == "__main__" :
        
        with open(sys.argv[1], 'r') as config_in :
                config_dict = json.load(config_in)
                

        config_dict["data"]["aln"] = isoformAln(config_dict["data"]["aln"], sys.argv[2])

        with open(sys.argv[1], 'w') as config_out :
                json.dump(config_dict, config_out)
