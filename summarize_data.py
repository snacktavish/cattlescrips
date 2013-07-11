##cattle disambiguation

import glob
import re
import numpy

##NEED TO SUMMARIZE pca info in this document too!!!
## Associate with each indidv.
##What about 55k Data?
##rerun chrm by chrm?
## make 55k PCA figure!


##fis=glob.glob("../RawData_PNAS/*.txt")

##cattle=pca.count_missing(fis)

class cow(object): 
    """Class docstring."""
    def __init__(self, name):
        """Method docstring."""
        self.name=name
        self.breednum = name[0:3]
        self.breedname = "n/a"
        self.struct3kperc=None
        self.kperc=None

cattle=[]
inds=open("../ids_full_miss.csv",'r').readlines()
ind_info={}
for ind in inds[1:]:
    cattle.append(cow(ind.split(',')[0]))
    ind_info[ind.split(',')[0]]=ind.split(',')[1:]

for ind in cattle:
  ind.breednum = ind_info[ind.name][2]
  ind.region = ind_info[ind.name][5]

  


breeds={}
for item in cattle:
  if item.breednum not in breeds:
    breeds[item.breednum]=1
  else:
    breeds[item.breednum]=breeds[item.breednum]+1


pattern=re.compile(r"\d\d\d:")  
patt=re.compile(r" ?\d+ \d{9}")

cows=dict((item.name,[])for item in cattle)
for fi in glob.glob("struct_nov/*_2?_f"):
  fii=open(fi)
  for lin in fii:
    if lin.startswith('100:'):
      lii=lin.split()
      if lii[1]>lii[2]:
         tau=1
      else:
         tau=2
    if patt.search(lin):
#      print(lin)
      lii=lin.split()
      try:
         cows[lii[1]].append(float(lii[4+tau]))
      except: pass


pattern=re.compile(r"\d\d\d:")  
patt=re.compile(r" ?\d+ \d{9}")

cows3=dict((item.name,[])for item in cattle)
for fi in glob.glob("../struct_nov/*_3?_f"):
  fii=open(fi)
  for lin in fii:
    if lin.startswith('100:'):
      lii=list(float(item) for item in lin.split()[1:-1])
      tau=lii.index(max(lii))
    if lin.startswith('500:'):
      lii=list(float(item) for item in lin.split()[1:-1])
      afr=lii.index(max(lii))
    if lin.startswith('507:'):
      lii=list(float(item) for item in lin.split()[1:-1])
      ind=lii.index(max(lii))
    if patt.search(lin):
#      print(lin)
      lii=lin.split()
      try:
          cows3[lii[1]].append(float(lii[5+afr]))
#         cows3[lii[1]].append((float(lii[5+afr]),float(lii[5+ind]),float(lii[5+tau])))
      except: pass

cowverage={}
for it in cows:
   try:
       cowverage[it]=sum(cows[it])/len(cows[it])
   except:
       cowverage[it]=None

cowverage3={}
for it in cows3:
   try:
       cowverage3[it]=sum(cows3[it])/len(cows3[it])
   except:
       cowverage3[it]=None


for ind in cattle:
#    ind.st3k=cowverage[ind.name]
    ind.st3k3=cowverage3[ind.name]


#calculate % ancestry by blobbling sides...
#TREATS Angus, Hereford,Shorthorn, and Red Angus as pure taurus
# INDICUS Sahiwal, Gir,GUzera


def calc_perc(instr):
  fi=open(instr+".evec").readlines()
  tau=[]
  ind=[]
  for lin in fi:
    if lin.split()[-1] in ['100','101','124','129']:
      tau.append(float(lin.split()[1]))
    if lin.split()[-1] in ['506','507','505','503']:
      ind.append(float(lin.split()[1]))
  tavg=sum(tau)/len(tau)
  iavg=sum(ind)/len(ind)
  perc_dict={}
  for lin in fi:
    perc_dict[lin.split()[0]]=(float(lin.split()[1])-iavg)/(tavg-iavg)
  return perc_dict



perc_dict3k=calc_perc("eig/3k_breed2")
#perc_dict55k=calc_perc("eig/50k_breed2")

perc_dict55k=calc_perc("eig/BREEDS_F")


for ind in cattle:
  try:
    ind.pca55k=perc_dict55k[ind.name]
  except:
    ind.pca55k=None
  try:
    ind.pca3k=perc_dict3k[ind.name]
  except:
    ind.pca3k=None
    


#missdict3k=dict((ki,[]) for ki in breeds)
#missdict55k=dict((ki,[]) for ki in breeds)

st3k=dict((ki,[]) for ki in breeds)
perc55k=dict((ki,[]) for ki in breeds)
perc3k=dict((ki,[]) for ki in breeds)

st3k3=dict((ki,[]) for ki in breeds)

for cow in cattle:
#    missdict3k[cow.breednum].append(cow.tkmiss)
#    missdict55k[cow.breednum].append(cow.ffkmiss)
    st3k[cow.breednum].append(cow.st3k)
    perc55k[cow.breednum].append(cow.pca55k)
    perc3k[cow.breednum].append(cow.pca3k)
    st3k3[cow.breednum].append(cow.st3k3)

oufi=open("summaryAFR.csv",'w')
oufi=open("summaryFINAL.csv",'w')
oufi.write("breed, number, num3kstruct, num55k, numpca3k, struct 3k avg, stdv,avg pca 55k, sd 55,avg pca 3k, sd 3 \n")
for ki in breeds:
#for ki in ['504']:
#  avg=sum(missdict3k[ki])/len(missdict3k[ki])#
#  ffs=[]
#  for item in missdict55k[ki]:
#    if item: 
#       ffs.append(item)
#  avg55=sum(ffs)/len(ffs)
  sts3=[]
  for item in st3k3[ki]:
    if item:
      sts3.append(item) 
  avgst3k3=sum(sts3)/len(sts3)
  sts=[]
  for item in st3k[ki]:
    if item:
      sts.append(item) 
  avgst3k=sum(sts)/len(sts)
  ff=[]
  for item in perc55k[ki]:
    if item: 
       ff.append(item)
  avgpca55=sum(ff)/len(ff)
  ff3=[]
  for item in perc3k[ki]:
    if item: 
       ff3.append(item)
  avgpca3=sum(ff3)/len(ff3)
#  oufi.write(",".join([ki, str(breeds[ki]),str(avg)[:4],str(avg55)[:4],str(len(ffs)),str(avgst3k)[:5],str(numpy.std(st3k[ki]))[:5],str(avgpca55)[:5],str(numpy.std(ff))[:5],str(avgpca3)[:5],str(numpy.std(ff3))[:5],"\n"]))
  oufi.write(",".join([ki, str(breeds[ki]), str(len(sts)),str(len(ff)),str(len(ff3)),str(avgst3k)[:5],str(numpy.std(sts))[:5],str(avgpca55)[:5],str(numpy.std(ff))[:5],str(avgpca3)[:5],str(numpy.std(ff3))[:5],"\n"]))


oufi.close()

oufi=open("all_ind_info.csv",'w')
oufi.write("id number, breed code, Structure 1.8k, PCA 1.8k, PCA 55k, rawdata file, chip, region, country, in 55k?\n")
for item in cattle:
    lin=[item.name, item.breednum, str(item.st3k), str(item.pca3k), str(item.pca55k), ind_info[item.name][0][1:-1],ind_info[item.name][3][1:-1],ind_info[item.name][4][1:-1],ind_info[item.name][5][1:-1],ind_info[item.name][-1][1:-2].strip()]
    oufi.write(','.join(lin)+'\n') 


oufi.close()



'''
oufi=open("newworld.csv",'w')
nwd=['121','128','111','900','901']
oufi.write(','.join(['breednum','name','pca3k','pca55k','st3k'])+'\n')
for cow in cattle:
  if cow.breednum in nwd:
    oufi.write(','.join([cow.breednum,cow.name,str(cow.pca3k),str(cow.pca55k),str(cow.st3k)])+'\n')


oufi.close()



test=[]
for item in breeds: test.append(breeds[item])



pattern=re.compile(r"\d\d\d:")  
patt=re.compile(r" ?\d+ \d{9}")

cows=dict((item.name,[])for item in cattle)
for fi in glob.glob("../May16_3k/*2?_f"):
  fii=open(fi)
  for lin in fii:
    if lin.startswith('100:'):
      lii=lin.split()
      if lii[1]>lii[2]:
         tau=1
      else:
         tau=2
    if patt.search(lin):
      lii=lin.split()
      cows[lin.split()[1]].append(float(lin.split()[4+tau]))


cowverage=dict((it,sum(cows[it])/len(cows[it]))for it in cows)


for cow in cattle:
    cow.st3k=cowverage[cow.name]



missdict3k=dict((ki,[]) for ki in breeds)
missdict55k=dict((ki,[]) for ki in breeds)
st3k=dict((ki,[]) for ki in breeds)



for cow in cattle:
    missdict3k[cow.breednum].append(cow.tkmiss)
    missdict55k[cow.breednum].append(cow.ffkmiss)
    st3k[cow.breednum].append(cow.st3k)



for ki in breeds:
  avg=sum(missdict3k[ki])/len(missdict3k[ki])
  avg55=sum(missdict3k[ki])/len(missdict3k[ki])
  avgst3k=sum(st3k[ki])/len(st3k[ki])
  print(ki+' # '+ str(breeds[ki])+", 3kmissing = "+str(avg)[:4]+", 55kmissing = "+str(avg55)[:4], "struct 3k average = "+str(avgst3k)[:4])



(ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=pca.make_dicts()
(rev_dict3k, forw_dict3k, locdict3k)=pca.make3k_revdict()
raw_3ksnp=forw_dict3k.keys()
raw_55ksnp=forw_dict.keys()


nam_chrm_dict={}
fil=open("SNP_Map.txt", 'r')
fi=fil.readlines()
fi=fi[1:]


for lin in fi:
  nam_chrm_dict[lin.split('\t')[1].strip('"')]=(lin.split('\t')[2],lin.split('\t')[3])


outp={}
for kii in ibmclocdict:
  outp[kii]=[]
  nam=ibmclocdict[kii]
  if nam == '':
      nam='*'
  outp[kii].append(nam)
  if nam in locdict:
    outp[kii].append(locdict[nam])
  else:
    outp[kii].append("*")
  if nam in locdict3k:
    outp[kii].append(locdict3k[nam])
  else:
    outp[kii].append("*")
  if nam in nam_chrm_dict:
    outp[kii].append(nam_chrm_dict[nam][0])
    outp[kii].append(nam_chrm_dict[nam][1])
  else:
    outp[kii].append("*")
    outp[kii].append("*")




fi=open("all_markers_chrm.csv",'w')
fi.write(" RawData #,bovinesnp50_name, SNP_Map(umd3), SNP_MAP3k(umd3), Chrm, Pos\n")
for kii in outp:
  fi.write(kii+',')
  fi.write(",".join(outp[kii]))
  fi.write('\n')


fi.close()


chrms=dict((str(ki),0) for ki in range(1,30))
chrms['"X"']=0
for item in outp:
  try:
    chrms[outp[item][3]]=chrms[outp[item][3]]+1
  except:
    print(outp[item][3])

fi=open("chrm_count.csv",'w')
for item in chrms:
  fi.write(','.join([item,str(chrms[item])])+'\n')


fi.close()

class breed(object):     
    def __init__(self, name, tkmiss=None, ffkmiss=None, tkstruct_tau=None,members=set() ):
        """Method docstring."""
        self.name=name
        self.tkmiss = arg1
        self.ffkmiss = arg2
        self.breednum = name[0:3]
        self.breedname = "n/a"
        tkstruct_tau=None
        self.members=members
'''
