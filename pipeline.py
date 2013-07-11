##LGH struct ms
import glob

##Read in marker files
autos=open("../autosomes.csv").readlines()

nope=open("../snp50_removed.csv").readlines()
remove=set()
for lin in nope:
  if lin.split(',')[0][0] not in ['4','5']:
     remove.add(lin.split(',')[2])

snpnames=set()
for item in autos:
  lii=item.split(',')
  if lii[1] not in remove:
    snpnames.add(lii[1])


(ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=pca.make_dicts()


##Read in list of indivduals
indfi=open("../ids.csv").readlines()
inds=[]
for lin in indfi[1:]:
  inds.append(lin.split(",")[0])q


genos={}
for ind in inds:
  genos[ind]={}


##Run through raw data files for individuals to test mnssingness:

infiles=["../RawData_PNAS/genos_116.txt"]
fis=glob.glob("../RawData_PNAS/*.txt")
for fi in fis:
  fil=open(fi)
  for lin in fil:
      lii=lin[1:-2].split(',')
      snpnam=ibmclocdict[lii[0].strip()]
      if snpnam in snpnames:
          genos[lii[1].strip()][snpnam]=lii[2].strip()
    




class indiv(object):   
    """Class docstring."""
    def __init__(self, name, arg1=None, arg2=None):
        """Method docstring."""
        self.geno={}


##Drop indivs missing >10%
miss_inds=[]
for ki in genos:
  if len(genos[ki]) < len(snpnames) * 0.9:
    miss_inds.append(ki)


for ind in miss_inds:
    del genos[ind]


##Run through SNPs. Remove any missing >10%
snp_miss={}
for snp in snpnames:
    snp_miss[snp]=0


for item in genos:
  for ite in genos[item]:
    snp_miss[ite]=snp_miss[ite]+1

##Make structure input files.


##make phase input files



##phase data



##run PCA


##estimate indivdual ancestry