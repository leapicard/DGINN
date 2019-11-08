import os, logging, json

def treerecsInfo(path, duplication):
	with open(path, "r") as recontree:
		nbDup = recontree.readline().split(",")[1].split("=")[1]
	head = "#======Treerecs=====#\n"
	data = "Duplication: {:d}\nNb (>8): {:d}\n".format(nbDup, len(duplication))

	return head+data

def gardInfo(lGard):

	head = "#=====GARD=====#\n"
	text = ""
	for kh in lGard:
		lBP = []
		with open(kh, "r") as f:
			lLine = f.readlines()
		
		index = 0
		while lLine[index].startswith("Breakpoint") != True:
			index += 1

		nbBP = 0
		if lLine[index+1] != "":
			index += 1
			while lLine[index].startswith(" "):
				line = [float(item.strip()) for item in lLine[index].split("|")]
				if line[2] < pvalue and line[4] < pvalue:
					lBP.append(int(line[0]))
				nbBP += 1
				index += 1

		text += ">{:s} :\nNb BreakPoint: {:d}\nBreakPoints significatives: {:s}\n".format(kh.split("/")[-1].split(".")[0], nbBP, ",".join(lBP))

	return head+text

def bustedInfo(dico):
	head = "#======BUSTED=======#\n"
	text = ""
	for aln in dico:
		path = aln+".BUSTED.json"
		with open(path, "w") as busted:
			busted = busted.readlines()

		i=0
		while "p-value" not in busted[i]:
			i+=1
		text += ">{:s}\np-value: {:d}\n".format(aln, busted[i].split(":")[1])

	return head+text

def memeInfo(dico):
	head = "#======MEME=======#\n"
	text = ""
	for aln in dico:
		path = aln+".MEME.json"
		with open(path, "r") as busted:
			busted = json.load(busted)

		print(busted["MLE"])





