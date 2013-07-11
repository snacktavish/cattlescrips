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
rev_num_dict={}
fil=open("ibmc_markers_041709.csv", 'r')
for lin in fil:
        row=lin.split(',')
        num=row[1]
        nam=row[9]
        ibmclocdict[num]=nam
        rev_num_dict[nam]=num

fil.close()


# <codecell>


def gen(X,Y): #Decode's AB in genotypes
    if X=='A' and Y =='A': return '1'
    if X=='B' and Y =='B': return '2'
    if X=='A' and Y =='B': return '3'
    if X=='B' and Y =='A': return '3'


##Read in list of indivduals
indfi=open("ids_full_miss.csv").readlines()
cownames={}
inds=[]
for lin in indfi[1:]:
  inds.append(lin.split(",")[0])
  try: 
     cownames[lin.split(',')[7].strip()[1:-1]]=lin.split(",")[0]
  except: pass


genos={}
for ind in inds:
  genos[ind]={}


#NOw, for my newest trick,
#add in all the final report genos.
prob_inds=set()
fis=glob.glob("final_reports/*.txt")
for fi in fis:
  fil=open(fi).readlines()
  assert(fil[9].split("\t")[0]=="SNP Name")
  for lin in fil[10:]:
      lii=lin.split("\t")
      try:
           genos[cownames[lii[1]]][lii[0]]=gen(lii[6],lii[7])
      except:
            prob_inds.add(lii[1])


fi="final_reports/Texas Longhorn 21oct2009_FinalReport AB.txt2"
fil=open(fi).readlines()
assert(fil[9].split("\t")[0]=="SNP Name")
for lin in fil[10:]:
      lii=lin.split("\t")
      try:
           genos[cownames[lii[1]]][lii[0]]=gen(lii[2],lii[3])
      except:
            prob_inds.add(lii[1])




inds={}
for ind in genos:
  if len(genos[ind])>1:
    inds[ind]=set()


outfi=open('RawData_PNAS/genos_999111Hyb.txt', 'w')
for ind in inds:
  for snp in genos[ind]:
    try:
      num=mrev_num_dict[snp]
      inds[ind].add(num)
      outfi.write('['+str(num)+', '+ind+', '+str(lin[1])+']'+'\n')
    except:
      pass


outfi.close()


def(final)

infi=open('''/home/ejbm/Documents/LGH3K/Cattleman's TLR 07feb2011_FinalReport.txt''','r')

rawdat=[]
for lin in infi.readlines()[9:]:
  rawdat.append(lin.strip().split())


infi.close()
infi=open('/home/ejbm/Documents/LGH3K/10Mar/Bovine 3K 10mar2011/CTL 07Mar2011_FinalReport.txt', 'r')
for lin in infi.readlines()[9:]:
  rawdat.append(lin.strip().split())



nudat=[]
nudat.append(["snp name", "sample id", "allele1 - ab", "allele1 - ab"])
for lin in rawdat:
    if lin[1]=='282':
      nudat.append([lin[0],999999005,lin[9], lin[8]])
    if lin[1]=='Dublin':
      nudat.append([lin[0],999999006,lin[8], lin[7]])
    if lin[1]=='Herd':
      nudat.append([lin[0],999999007,lin[8], lin[7]])


locations=[]
locdict={}
fi=open("/home/ejbm/Documents/pyScrips/ibmc_markers_041709.csv", 'r')
for lin in fi:            
	row=lin.split(',')
	try:
		num=int(row[1])
	except:
		num=row[1]
	nam=row[0]
	try:
		chrm=int(row[37])
	except:
		chrm=row[37]
	try:
		pos=int(row[38])
	except:
		pos=row[38]
	locs=[num, nam, chrm, pos]
	locations.append(locs)
	locdict[nam]=num
	##'''if num in colsort: ##this is key- excludes SNPS that are not in the struct data file'''


outfi=open('new_999.txt', 'w')
for lin in nudat:
  try:
    num=locdict[lin[0]]
    if lin[2]==lin[3]=='A':
      geno='1'
    if lin[2]==lin[3]=='B':
      geno='2'
    if lin[2]!=lin[3]:
      geno='3'
    outfi.write('['+str(num)+', '+str(lin[1])+', '+geno+']'+'\n')
  except: pass


outfi.close()  
