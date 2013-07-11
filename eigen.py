import sys
sys.path.append("./PCA")
import pca



#calculate % ancestry by blobbling sides...

def calc_perc(intsr):
  fi=open(instr+".evec").readlines()
  tau=[]
  ind=[]
  
  for lin in fi:
    if lin.split()[-1] in ['100','101','104','113','124','127','129','214','218']:
      tau.append(float(lin.split()[1]))
    if lin.split()[-1] in ['506','507','505','503']:
      ind.append(float(lin.split()[1]))
  tavg=sum(tau)/len(tau)
  iavg=sum(ind)/len(ind)
  perc_dict={}
  for lin in fi:
    perc_dict[lin.split()[0]]=(float(lin.split()[1])-iavg)/(tavg-iavg)