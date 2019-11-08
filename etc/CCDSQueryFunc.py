import pandas, os, requests, argparse

def dict2fasta(dico):
	"""
	Function converting a dictionary associating gene names to their sequences into a writable Fasta format.

	@param dico: Dictionary associating gene names (keys) to their CCDS fasta sequence (values)
	@return txtoutput: Fasta formatted string of the dictionary
	"""

	txtoutput = ""
	for key, value in dico.items():
		txtoutput += ">{:s}\n{:s}\n".format(str(key),str(value))

	return(txtoutput)

def check(arg):	
	"""
	Function checking whether a path to a file exists.
	
	@param arg: Path to a file
	@return os.path.abspath(arg): Normalized/absolute version of the path
	"""

	if not os.path.exists(arg):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: arg does not exist, please enter valid argument.
		raise argparse.ArgumentTypeError("\"{0}\" does not exist, please enter valid path.".format(arg))
	
	return os.path.abspath(arg)

def getCCDS(genesFile, species, spName):
	"""
	Function downloading the input list of genes' CCDS sequences from the ENSEMBL server, given their CCDS accession.
	Fonction qui permet de récupérer les séquences des gènes à partir du nom du gène.

	@param obj: A basicData object
	@return1 allFasta: Path to the Fasta file containing all the genes' CCDS sequences
	@return2 dId2Seq: Dictionary associating gene identifiers (keys) to their fasta sequences (values)
	"""
	
	dId2Seq = {} # key = seqID (format speSpe|gene|CCDS012345) value = fasta sequence

	server = "http://rest.ensembl.org"
	species = species.lower().replace(" ", "_")
	
	df = pandas.read_csv(genesFile, sep='\t') # import gene list as pandas dataframe
	df = df.drop_duplicates(subset='Approved symbol', keep='first') # keep only first instance of each gene
	
	for index, row in df.iterrows():
		geneName = str(row['Approved symbol'])
		geneCCDS = str(row['CCDS accession'])
		seqID = "{:s}_{:s}_{:s}".format(spName, geneName, geneCCDS)
		#print(seqID)
		
		# download fasta sequence from ENSEMBL server
		add = 1
		while add <= 3:
			ext = "/sequence/id/{:s}.{:d}?content-type=text/x-fasta;species={:s};object_type=transcript;db_type=otherfeatures;type=cds".format(geneCCDS, add, species)	#db_type=otherfeatures;type=cds
			link = server+ext
			print(link)
			r = requests.get(link, headers={"Content-Type" : "text/x-fasta"})
			text = "".join(r.text.split("\n")[1:])

			if not r.ok or text == "":
				add += 1
			else:
				add = 10

		if add == 10:
			#print(geneName)
			seqID = seqID.split(".")[0]
			dId2Seq[seqID] = text
			#print(geneName)
			#print(len(text))
		else:
			print("Couldn't download sequence for {:s}".format(geneName))
	
	for key, value in dId2Seq.items():
		out = ".".join(genesFile.split(".")[:-1])+"_"+key.split("_")[1]+"_CCDS.fasta"
		with open(out, "w") as fasta:
			if key != "" or value != "":
				fasta.write(">"+key+"\n"+value)
		
	print("Wrote CCDS sequences to their respective files for {:d} genes.".format(len(dId2Seq.keys())))
