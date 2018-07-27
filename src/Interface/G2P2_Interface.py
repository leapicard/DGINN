#!/usr/bin/python3
# -*- coding: utf-8 -*-
from tkinter import messagebox, filedialog, StringVar, font
from tkinter import *
from time import localtime, strftime
import os, shlex, subprocess

class mainGUI:

	def __init__(self, master):

		self.master = master
		master.title("G2P Interface")

		menu = Menu(self.master)
		self.master.config(menu=menu)

		filename = StringVar()
		filename.set("No Path Found ")
		fileout = StringVar()
		fileout.set("None")
		dataBase = StringVar()
		dataBase.set("Obligatory")
		logger = StringVar()
		logger.set("")
		evalue = StringVar()
		coverage = StringVar()
		identity = StringVar()
		CheckVar = IntVar()
		CheckVar.set(False)
		queries = StringVar()
		queries.set("")
		apikey = StringVar()
		apikey.set("")


		file = Menu(menu)
		file.add_command(label="Open", command= lambda: self.openParamFile(filename, fileout, dataBase, logger, evalue, coverage, identity, CheckVar, queries, apikey, lStep))
		file.add_command(label="Exit", command=self.client_exit)
		menu.add_cascade(label="File", menu=file)


		helps = Menu(menu)
		helps.add_command(label="Help", command= lambda: messagebox.showinfo("Info", "G2P's Graphic interface by Léa Picard \nand Quentin Ganivet."))
		menu.add_cascade(label="Helps", menu=helps)



		infile = Button(master, text="In File")
		infile.grid(row=1, column=1)
		infile.configure(command= lambda: self.openfile(filename))

		
		infilePath = Entry(master, textvariable=filename)
		infilePath.grid(row=1, column=2, columnspan=8)


		text = Label(master, text= ">")
		text.grid(row=1, column=10)

		out = Button(master, text="Out Directory")
		out.grid(row=1, column=12)
		out.configure(command= lambda: self.openDir(fileout))

		outfile = Entry(master, textvariable=fileout)
		outfile.grid(row=1, column=13, columnspan=8)

		db = Label(master, text="Database :")
		db.grid(row=2, column=1, sticky=W)

		dbEntry = Entry(master, textvariable=dataBase)
		dbEntry.grid(row=2, column=2, columnspan=2)
		f = font.Font(db, db.cget("font"))
		f.configure(underline=True)
		db.configure(font=f)

		log = Label(master, text="Logger's Name :")
		log.grid(row=3, column=1, sticky=W)
		f = font.Font(log, log.cget("font"))
		f.configure(underline=True)
		log.configure(font=f)

		logEntry = Entry(master, textvariable=logger)
		logEntry.grid(row=3, column=2, columnspan=2)


		ev = Label(master, text="E-value :")
		ev.grid(row=4, column=1, sticky=W)
		f = font.Font(ev, ev.cget("font"))
		f.configure(underline=True)
		ev.configure(font=f)

		evEntry = Entry(master, textvariable=evalue)
		evEntry.grid(row=4, column=2, columnspan=2)

		cov = Label(master, text="Coverage :")
		cov.grid(row=5, column=1, sticky=W)
		f = font.Font(cov, cov.cget("font"))
		f.configure(underline=True)
		cov.configure(font=f)

		covEntry = Entry(master, textvariable=coverage)
		covEntry.grid(row=5, column=2, columnspan=2)

		ID = Label(master, text="Identity :")
		ID.grid(row=6, column=1, sticky=W)
		f = font.Font(ID, ID.cget("font"))
		f.configure(underline=True)
		ID.configure(font=f)

		idEntry = Entry(master, textvariable=identity)
		idEntry.grid(row=6, column=2, columnspan=2)


		step = Label(master, text="First Step :")
		step.grid(row=2, column=12, sticky=W)
		f = font.Font(step, step.cget("font"))
		f.configure(underline=True)
		step.configure(font=f)

		lStep = Listbox(master, height=3)
		lStep.grid(row=2, column=13, rowspan=3, columnspan=2)
		lStep.insert(END, "Blast")
		lStep.insert(END, "Accession")
		lStep.insert(END, "Fasta")
		lStep.insert(END, "ORFs")
		lStep.insert(END, "PRANK")
		lStep.insert(END, "PhyML")
		lStep.insert(END, "Tree")
		sb = Scrollbar(master, orient=VERTICAL)
		sb.grid(row=2, column=15, rowspan=3, sticky=NS, pady=5)
		sb.configure(command=lStep.yview)
		lStep.configure(yscrollcommand=sb.set, exportselection=False)

		remote = Label(master, text="Remote :")
		remote.grid(row=5, column=12, sticky=W)
		f = font.Font(remote, remote.cget("font"))
		f.configure(underline=True)
		remote.configure(font=f)

		Radiobutton(master, text="True", variable=CheckVar, value=True).grid(row=5, column=13)
		Radiobutton(master, text="False", variable=CheckVar, value=False).grid(row=5, column=14)

		query = Label(master, text="Query :")
		query.grid(row=6, column=12, sticky=W)
		f = font.Font(query, query.cget("font"))
		f.configure(underline=True)
		query.configure(font=f)

		queryEntry = Entry(master, textvariable=queries)
		queryEntry.grid(row=6, column=13, columnspan=2)

		apik = Label(master, text="API Key :")
		apik.grid(row=7, column=12, sticky=W)
		f = font.Font(apik, apik.cget("font"))
		f.configure(underline=True)
		apik.configure(font=f)

		apikEntry = Entry(master, textvariable=apikey)
		apikEntry.grid(row=7, column=13, columnspan=2)		


		Start = Button(master, text="Begin")
		Start.grid(row=9, column=10)
		Start.configure(command= lambda: self.start(filename.get(), fileout.get(), dataBase.get(), logger.get(), evalue.get(), coverage.get(), identity.get(), CheckVar.get(), queries.get(), apikey.get(), lStep.curselection()))



	def openParamFile(self, filename, fileout, dataBase, logger, evalue, coverage, identity, CheckVar, queries, apikey, lStep):
		path = filedialog.askopenfilename(parent=root)

		lParams = ["infile", "blastdb", "outdir", "logfile", "evalue", "mincov", "perc_id", "step_id", "remote", "entry_query", "API_Key"]
		dStep = {"blast":0, "accession":1, "fasta":2, "orf":3, "prank":4, "phyml":5, "tree":6}
		dParam = {}

		with open(path, "r") as param:
			for line in param:
				temp = line.split(":")
				dParam[temp[0]] = temp[1].strip("\n")
				if temp[0] not in lParams:
					messagebox.showerror("ERROR", "Parameters "+temp[0]+" doesn't exist !")
					break

		filename.set(dParam["infile"])
		fileout.set(dParam["outdir"])
		dataBase.set(dParam["blastdb"])
		logger.set(dParam["logfile"])
		evalue.set(dParam["evalue"])
		coverage.set(dParam["mincov"])
		identity.set(dParam["perc_id"])

		if dParam["remote"] == "":
			CheckVar.set(False)
		elif dParam["remote"].lower() != "true" and dParam["remote"].lower() != "false" :
			messagebox.showerror("ERROR", "Remote option need to be equal to True or False !")
		else:
			CheckVar.set(dParam["remote"])
		queries.set(dParam["entry_query"])
		apikey.set(dParam["API_Key"])

		if dParam["step_id"] not in dStep and dParam["step_id"] != "":
			messagebox.showerror("ERROR", "The "+dParam["step_id"]+" step doesn't exist !")
		elif dParam["step_id"] == "":
			lStep.selection_set(0)
		else:
			lStep.selection_set(dStep[dParam["step_id"]])


	def client_exit(self):
		exit()


	def openfile(self, filename):
		path = filedialog.askopenfilename(parent=root)
		filename.set(path)

	def openDir(self, fileout):
		path = filedialog.askdirectory(parent=root)
		fileout.set(path)


	def start(self,filename, fileout, dataBase, logger, evalue, coverage, identity, CheckVar, queries, apikey, step):
		args = {"in":filename, "out":fileout, "db":dataBase, "log":logger, "ev":evalue, "cov":coverage, "id":identity, "rmt":CheckVar, "query":queries, "apik":apikey, "step":step}
		lSteps = ["blast", "accession", "fasta", "orf", "prank", "phyml", "tree"]

		if os.path.exists(args["in"]) == False:
			messagebox.showerror("ERROR", "The source file doesn't exist !")

		elif args["db"] == "" or args["db"] == "Obligatory !!":
			messagebox.showerror("ERROR", "No DataBase specified !")

		else:
			if os.path.exists(args["out"]) == False:
				messagebox.showwarning("WARNING !", "Out directory doesn't exist, we will use the default directory....")

			if args["step"] == ():
				messagebox.showwarning("WARNING !", "No Step precised, \nwe will execute the pipeline from the Blast Step....")
				args["step"] = "blast"
			else:
				args["step"] = lSteps[args["step"][0]]

			commandLine = "pyhton3 G2P.py"
			for i in args:
				if args[i] != "":
					if args[i] == 1:
						args[i] = "False"
					elif args[i] == 0:
						args[i] = "True"

					commandLine += " -"+i+" "+args[i] 
			
			print(commandLine)

			lCmd = shlex.split(commandLine)
			try:
				run = subprocess.run(lCmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			except:
				messagebox.showerror("ERROR", "the pipeline did not run properly !...")

root = Tk()

root.geometry("650x200")

app = mainGUI(root)
root.mainloop()