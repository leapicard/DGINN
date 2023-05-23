import glob
import os, sys

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
    print("""Usage:  python3 recup_to_parse.py [-o outfile]

    where:

    -o outfile: output file that will parsed by parseResults.py
    
    """)

    return

  ## Optional output file (default toparse.txt)
  output=parseargv(argv,"-o","output","toparse.txt")

  
  lf = [x for x in glob.glob("*.log") if x.find("_DGINN")!=-1]
  lf.sort(key = lambda t: -os.stat(t).st_mtime) ## sort to get most recent ones

  llog=[]
  lpref=[]

  for f in lf:
    pD=f.rfind("_DGINN")
    pref=f[:pD]
    if not pref in lpref:
      flag=0
      fin=open(f,"r")
      for l in fin.readlines():
        if l.find("positiveSel")!=-1 and (l.find("Output dir")!=-1 or l.find("Alignement:")!=-1):
          if flag==0:
            lpref.append(pref)
            llog.append(f)
          flag=1
        if l.find("Finished DGINN")!=-1:
          if flag:
            flag=2
      fin.close()
      print("-~+"[flag] + " " + pref)
      
  fout=open(output,"w")
  for log in llog:
    fin=open(log,"r")
    flag=0
    for l in fin.readlines():
      if l.find("positiveSel")!=-1 and (l.find("Output dir")!=-1 or l.find("Alignement:")!=-1):
        ll=l.split()
        if flag==0:
          fout.write(ll[-1]+"\t")
          flag=1
        elif flag==1:
          fout.write(ll[-1]+"\n")
          flag=2
  
    fin.close()
    
  fout.close()

if __name__=="__main__":
  argv=sys.argv
  main(argv)
