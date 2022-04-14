import os, subprocess
import sys, glob

### Parser
##
## If there is no value argument for the option, if default equals
## True or False, the parser returns the negative of the default.

def parseargv(argv, arg, name, default=""):
  try:
    posi=argv.index(arg)
  except ValueError:
    return default
  else:  
    if posi==len(argv)-1 or argv[posi+1][0]=="-":
      if default!=True and default!=False:
        print("Missing ", name, " name after ", arg, " option.")
        sys.exit(0)
      else:
        return not default
    else:
      return(argv[posi+1])


def main(argv):
  if parseargv(argv,"-h","",False):
    print("""Usage:  python3  multi_dginn.py dataname -p parameters [-i image][-v][-j jobs]

    where:

       dataname: name of the file where the input names are stored per
                 line (aka used as argument of --infile option of
                 DGINN.py).

       -p parameters: name of the DGINN parameters file (aka used with
                      option -p in DGINN.py).

       -i image: name of the docker image used (default lpicard/dginn).

       -j jobs: number of jobs used in parallel (default 4)

       -v  verbosity of the commands.   (default False)
      """)
    return

  
  if len(argv)<2:
    print("Missing file name of data.")
    sys.exit(0)
    
  ### data file name
  fn=argv[1]
  
  ldata=[]
  try:
    fdata=open(fn,"r")
    ldata=list(map(lambda x:x.strip(),fdata.readlines()))
    fdata.close()
  except:
    print("Unknown data file " + fn)
    sys.exit(0)
                
  
  
  ## Optional image name (default lpicard/dginn)
  image=parseargv(argv,"-i","image","lpicard/dginn")
  
  print("Image name: ", image)
  
  ## Optional image name (default lgueguen/dginn)
  parameters=parseargv(argv,"-p","parameters","")
  
  if parameters=="":
    print("Missing parameter file.")
    sys.exit(0)
  
  if not os.path.exists(parameters):
    print("Unknown parameter file ",parameters)
    sys.exit(0)
    
  ## Number of parallel running jobs (default 4)
  njobs = int(parseargv(argv,"-j","number of jobs",4))
  
  ## Verbosity
  verbose = parseargv(argv,"-v","verbose",False)
  
  ################
  
  ### List running containers 
  lps=subprocess.run(["docker","ps"],stdout=subprocess.PIPE).stdout.decode("utf-8").strip().split("\n")[1:]
  run=[ps.split()[-1] for ps in lps]
  lrun=len(run)
  
  print("Number of docker running processes: ", lrun)
  
  if lrun>=njobs:
    print("Running jobs at max ", njobs)
    sys.exit(1)
  
  PWD=os.getcwd()

  for data in ldata:
      running=False
      pref=data[:data.rfind(".")]
      lpref=glob.glob(pref+"_DGINN_*.log")
      if len(lpref)>=1:
          if pref in run:
              running=True
          lpref.sort(key = lambda t: -os.stat(t).st_mtime) ## sort to get most recent ones
          f=open(lpref[0],"r")
          l = list(f.readlines())[-1].strip()
          print("sr"[running] + "at " + l)
          if running or l.find("Finished DGINN analyses")!=-1:
            print("- " + pref)
            continue
      print("+ " + pref)
      lcmd=["docker","run","--rm","--name=%s"%pref, "-u 1000:1000", "-e HOME=.", "-v %s:%s"%(PWD,PWD),"-w %s"%(PWD),str(image),"-p", parameters,"--infile %s"%data,"--outdir %s_Res"%(pref)]
      if verbose:
        print(" ".join(lcmd))
      res=subprocess.Popen(" ".join(lcmd)+ " &> %s.out &"%(pref),shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
      lrun+=1
      if lrun>=njobs:
        print("Running jobs at max ", njobs)
        sys.exit(1)


if __name__=="__main__":
  argv=sys.argv
  main(argv)
  

