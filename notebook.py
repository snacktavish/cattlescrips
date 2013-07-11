# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

##LGH struct ms
import glob

# <codecell>

##Make location dictionaries
ibmclocdict={} ##takes an ibmc number (AS IN RAW DATA FILES)  and outputs SNP name
    fil=open("ibmc_markers_041709.csv", 'r')
    for lin in fil:
        row=lin.split(',')
        num=row[1]
        nam=row[9]
        ibmclocdict[num]=nam
    fil.close()

# <codecell>

locations=open("../UMC_snp50_markers_121007.csv")
chrm_dict={}
chrm_n=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29","30","31","99",""]
for item in chrm_n:
        chrm_dict[item]={}


for lin in locations:
    try:
        lii=lin.split(",")
        nam=lii[1]
        chrm=lii[6]
        try:
            pos=int(lii[7])
        except:
            pos=0
        chrm_dict[chrm][pos]=nam
    except: print lin


for item in chrm_n:
        chrm_dict[item].keys().sort()





        

# <codecell>

##Read in marker files
nope=open("../snp50_removed.csv").readlines()
remove=set()
for lin in nope:
  if lin.split(',')[0][0] not in ['4','5']:
     remove.add(lin.split(',')[2])


autos=open("../autosomes.csv").readlines()
snpnames=set()
for item in autos:
  lii=item.split(',')
  if lii[1] not in remove:
    snpnames.add(lii[1])

# <codecell>

##Read in list of indivduals
indfi=open("../ids.csv").readlines()
inds=[]
for lin in indfi[1:]:
  inds.append(lin.split(",")[0])


genos={}
for ind in inds:
  genos[ind]={}

# <codecell>

##Run through raw data files for individuals to test mnssingness:
#fis=glob.glob("../RawData_PNAS/*.txt")
fis=["../RawData_PNAS/genos_116.txt"]
for fi in fis:
  fil=open(fi)
  for lin in fil:
      lii=lin[1:-2].split(',')
      snpnam=ibmclocdict[lii[0].strip()]
      if snpnam in snpnames:
          genos[lii[1].strip()][snpnam]=lii[2].strip()
    


# <codecell>

##Drop indivs missing >10%
miss_inds=[]
for ki in genos:
  if len(genos[ki]) < len(snpnames) * 0.9:
    miss_inds.append(ki)

for ind in miss_inds:
    del genos[ind]

# <codecell>

##Run through SNPs. Remove any missing >10%
snp_miss={}
for snp in snpnames:
    snp_miss[snp]=0


for item in genos:
  for ite in genos[item]:
    snp_miss[ite]=snp_miss[ite]+1

snp_removal=set()
for snp in snp_miss:
    if snp_miss[snp]<len(genos)*0.9:
        snp_removal.add(snp)

filtered_snps=snpnames-snp_removal
print("number of snps is %i" %len(filtered_snps))

# <codecell>

##make phase input files

# <codecell>

##phase data

# <codecell>

##run PCA

# <codecell>

##estimate indivdual ancestry

# <codecell>


