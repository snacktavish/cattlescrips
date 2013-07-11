import glob
from __future__ import division



snps={}
for line in (open("RawData_PNAS/50ksnp.csv").readlines()):
   snps.[line.split(',')[0]]={}


inds={}
for fi in glob.glob("RawData_PNAS/*.txt"):
    nam=fi.split('_')[-1].split('.')[-2]
    print(nam)
    inds[nam]={}


indlist=set()
infi=open("ids_full_miss.csv").readlines()
for lin in infi:
   lii=lin.split(",")
   if lii[-1].strip()[2:-1]=="Y":
      try:
        nam=lii[2]
        ind=lii[0]
        indlist.add(ind)
        inds[nam][ind]={}
        inds[nam][ind]['1']=0
        inds[nam][ind]['2']=0
        inds[nam][ind]['3']=0
        inds[nam][ind]['10']=0
      except:
        print lii
counts
counts={}
for fi in glob.glob("RawData_PNAS/*.txt"):

for fi in fis[25:]:

for fi in ['RawData_PNAS/genos_504.txt']:
    nam=fi.split('_')[-1].split('.')[-2]  
    counts[nam]={}  
    print("made empty ind dict")
    snps={}
    for line in (open("RawData_PNAS/50ksnp.csv").readlines()):
         snps[line.split(',')[0]]=set()
    print("made empty snp dict")
    for lin in open(fi).readlines():
        lii=lin[1:-2].split(',')
        indnum=lii[1].strip()
        if lii[0].strip() in snps.keys() and indnum in indlist:
                inds[nam][indnum][lii[-1].strip()]+=1
                if lii[-1] not in snps[lii[0]]:
                      snps[lii[0]].add(lii[-1].strip())
    counts[nam]['invars']=0
    print("ok adding up!")
    c1=[]
    c2=[]
    c3=[]
    c10=[]
    for indi in inds[nam]:
        c1.append(inds[nam][indi]['1'])
        c2.append(inds[nam][indi]['2'])
        c3.append(inds[nam][indi]['3'])
        c10.append(inds[nam][indi]['10'])
    counts[nam]['1']=(sum(c1),numpy.std(c1))
    counts[nam]['2']=(sum(c2),numpy.std(c2))
    counts[nam]['3']=(sum(c3),numpy.std(c3))
    counts[nam]['10']=(sum(c10),numpy.std(c10))
    counts[nam]["num inds"]=len(inds[nam])
    for item in snps:
       if len(snps[item]&set(['1','2','3'])) < 2:
            counts[nam]['invars']+=1 



oufi=open("counts.csv",'w')
header=['breed','invars','avg hets', 'het s.d.', 'num inds']

oufi.write(",".join(header)+'\n')
for item in counts:
 try:
   oulin=[item]
   oulin.append(str(counts[item]['invars']))
   oulin.append(str(counts[item]['3'][0]))
   oulin.append(str(counts[item]['3'][1]))
   oulin.append(str(counts[item]['num inds']))
   oufi.write(','.join(oulin)+'\n')
 except:
   print(item)

oufi.close()
print(item,(counts[item]['3']/counts[item]['inds'])/len(snps))



