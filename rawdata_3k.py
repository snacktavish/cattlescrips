# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

##LGH struct ms
import glob
import operator
import os
import glob
import pickle

# <codecell>

##Make location dictionaries
ibmclocdict={} ##takes an ibmc number (AS IN RAW DATA FILES)  and outputs SNP name
ibmcrevdict={}
fil=open("ibmc_markers_041709.csv", 'r')
for lin in fil:
        row=lin.split(',')
        num=row[1]
        nam=row[9]
        ibmclocdict[num]=nam
        ibmcrevdict[nam]=num


fil.close()

# seperate out markers by chormosome
markers=open("../UMC_snp50_markers_121007.csv").readlines()
mark_dict={} #takes chrm and returns all markers
mark_chrm_dict={}#takes marker name and returns snp50 number, chromosome, and position.
mark_num_dict={}
mark_nam_dict={}
for lin in markers:
    lii=lin.split(',')
    chrm = lii[6]
    num = lii[0]
    nam = lii[1]
    try:
            pos=int(lii[7])
    except:
            pos=0
    mark_chrm_dict[nam]=(num,chrm, pos)
    mark_num_dict[nam]=num #translates to snp number
    mark_nam_dict[num]=nam #translates back to SNP name
    if chrm not in mark_dict:
    	mark_dict[chrm]=[]
    mark_dict[chrm].append((nam,pos))


#this should make it so that there is a dictionary where the key is 1-30, and the values are a list of names and positions in order
for item in mark_dict.keys():
      mark_dict[item].sort(key=operator.itemgetter(1))


# <codecell>

##Read in marker files
nope=open("../snp50_removed.csv").readlines()
remove=set()
for lin in nope:
  if lin.split(',')[0][0] not in ['4','5']:
     remove.add(lin.split(',')[2])


snpnames=set()
for chrm in list(range(1,30)):
  for mark in mark_dict[str(chrm)]:
    if mark[0] not in remove:
      snpnames.add(mark[0])

# <codecell>

##Read in list of indivduals
indfi=open("../ids_full.csv").readlines()
cownames={}
inds=[]
for lin in indfi[1:]:
  inds.append(lin.split(",")[0])
  try: 
     cownames[lin.split(',')[7].strip()[1:-1]]=lin.split(",")[0]
  except: pass

del cownames['x']


genos={}
for ind in inds:
  genos[ind]={}


def gen(X,Y): #Decode's AB in genotypes
    if X=='A' and Y =='A': return '1'
    if X=='B' and Y =='B': return '2'
    if X=='A' and Y =='B': return '3'
    if X=='B' and Y =='A': return '3'

    

#NOw, for my newest trick,
#add in all the final report genos.
fis=glob.glob("../final_reports/*.txt")
for fi in fis:
  fil=open(fi).readlines()
  assert(fil[9].split("\t")[0]=="SNP Name")
  for lin in fil[10:]:
      lii=lin.split("\t")
      if lii[0] in snpnames:
        try:
           genos[cownames[lii[1]]][lii[0]]=gen(lii[6],lii[7])
        except:
            pass



fi="../final_reports/Texas Longhorn 21oct2009_FinalReport AB.txt2"
fil=open(fi).readlines()
assert(fil[9].split("\t")[0]=="SNP Name")
for lin in fil[10:]:
      lii=lin.split("\t")
      if lii[0] in snpnames:
        try:
           genos[cownames[lii[1]]][lii[0]]=gen(lii[2],lii[3])
        except:
            pass



fi=open("all_ind_info.csv").readlines()
hybs=[]
nohybs=[]
for lin in fi:
     if lin.split(',')[5]=="'genos_999111Hyb.txt'":
             hybs.append(lin.split(',')[0])
     if lin.split(',')[5]=="'genos_999111.txt'":
             nohybs.append(lin.split(',')[0]) 
 

fi=open('genos_999111Hyb.txt','w')
for item in hybs:
	for snp in genos[item]:
		if genos[item][snp] in ['1','2','3'] and snp in ibmcrevdict:
			fi.write("[%s, %s, %s]\n"%(ibmcrevdict[snp],item,genos[item][snp]))

fi.close()

fi=open('genos_999111.txt','w')
for item in nohybs:
	for snp in genos[item]:
		if genos[item][snp] in ['1','2','3'] and snp in ibmcrevdict:
			fi.write("[%s, %s, %s]\n"%(ibmcrevdict[snp],item,genos[item][snp]))

fi.close()
