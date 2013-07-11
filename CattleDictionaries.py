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

def create_dictionaries(ibmc_csv = "./PCA/ibmc_markers_041709.csv", snp_50="UMC_snp50_markers_121007.csv"):
    ##Make location dictionaries
    ibmclocdict={} ##takes an ibmc number (AS IN RAW DATA FILES)  and outputs SNP name
    fil=open(ibmc_csv, 'r')
    for lin in fil:
        row=lin.split(',')
        num=row[1]
        nam=row[9]
        ibmclocdict[num]=nam
    fil.close()
    # seperate out markers by chormosome 
    markers=open(snp_50).readlines()
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
    return(ibmclocdict,mark_num_dict,mark_nam_dict,mark_dict)

# <codecell>

ibmclocdict,mark_num_dict,mark_nam_dict,mark_dict = create_dictionaries()

# <codecell>

#remove SNPs suggested by Taylor lab
nope=open("snp50_removed.csv").readlines()
remove=set()
for lin in nope:
  if lin.split(',')[0][0] not in ['4','5']:
     remove.add(lin.split(',')[2])

# <codecell>

snpnames=set()
for chrm in list(range(1,30)):
  for mark in mark_dict[str(chrm)]:
    if mark[0] not in remove:
      snpnames.add(mark[0])

# <codecell>

##Add in func here to add new indivs to full id file?

# <codecell>

def make_cownames(id_file="ids_full_chin.csv"):
    indfi=open(id_file).readlines()
    cownames={}
    inds=[]
    for lin in indfi[1:]:
        inds.append(lin.split(",")[0])
        try: 
            cownames[lin.split(',')[7].strip()]=lin.split(",")[0]
        except: pass
    del cownames['x']
    return(inds,cownames)

# <codecell>

inds,cownames=make_cownames()

# <codecell>

def gen(X,Y): #Decode's AB in genotypes
    if X=='A' and Y =='A': return '1'
    if X=='B' and Y =='B': return '2'
    if X=='A' and Y =='B': return '3'
    if X=='B' and Y =='A': return '3'
    else: return '10'

