# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

##LGH struct ms
import glob
import operator
import os

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
##Run through raw data files for individuals to test missingness:
fis=glob.glob("../RawData_PNAS/*.txt")
##fis=["../RawData_PNAS/genos_101.txt","../RawData_PNAS/genos_503.txt"]
for fi in fis:
  fil=open(fi)
  for lin in fil:
      print('.',end='')
      lii=lin[1:-2].split(',')
      snpnam=ibmclocdict[lii[0].strip()]
      if snpnam in snpnames and lii[2].strip() != '10':
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
    snp_miss[ite]=snp_miss[ite]+1 #counts up how many indivs are misisng that snp

snp_removal=set()
for snp in snp_miss:
    if snp_miss[snp]<len(genos)*0.9:
        snp_removal.add(snp)

filtered_snps=snpnames-snp_removal
print("number of snps is %i" %len(filtered_snps))


##### REMOVE THEM FORM THE mark_dict?!
##****FIND A BETTER WAY TO REMOVE ITEMS FROM DICT
def keep(x):
    return x[0] in filtered_snps


for chrm in mark_dict:
	mark_dict[chrm]= [x for x in mark_dict[chrm] if keep(x)] #strips snpa with too much missing data out of marker dictionary


##make phase input files
for chrm in mark_dict:
        print("Hromosome %s" %chrm)
        oufi=open("../phase_inputs/phase_chrm%s.inp"%chrm,'w')
        oufi.write(str(len(genos))+"\n")
        oufi.write(str(len(mark_dict[chrm]))+"\n"+"P ")
        for mark in mark_dict[chrm]:
             oufi.write(str(mark[1])+' ')
        oufi.write('\n')
        oufi.write("S" * len(mark_dict[chrm])+"\n")#WTF ARE THESE SSSSs for?
        for ind in genos:
                oufi.write("# id %s"%ind+"\n")
                hap1=[]
                hap2=[]
                for mark in mark_dict[chrm]:
                    try:
                      gen=genos[ind][mark[0]]
                    except KeyError:
                      gen='?'
                    if gen=='1':
                       hap1.append('1')
                       hap2.append('1')
                    if gen=='2':
                       hap1.append('2')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                       hap2.append('2')
                    if gen=='3':
                       hap1.append('1')
                       hap2.append('2')
                    if gen=='?' or gen=='10':
                       hap1.append('?')
                       hap2.append('?')
                    if gen not in ['1','2','3','?','10']:
                       hap1.append('?')
                       hap2.append('?')
                assert(len(hap1)==len(hap2)==len(mark_dict[chrm]))
                oufi.write("".join(hap1)+'\n')
                oufi.write("".join(hap2)+'\n')

##Actually phase the data
##MAKE THIS HAPPEN ON THE CLUSTER

for chrm in mark_dict:
    os.system("/home/ejbm/Documents/Phase/fastPHASE/fastPHASE_Linux -o../phase_outputs/chrm%s ../phase_inputs/phase_chrm%s.inp"%(chrm,chrm))##RECHECK PHASING SETTINGS

##for chrm in ['2']:
##    os.system("/home/ejbm/Documents/Phase/fastPHASE/fastPHASE_Linux -T2 -o../phase_outputs/chrm%s ../phase_inputs/phase_chrm%s.inp"%(chrm,chrm))##RECHECK PHASING SETTINGS


##BOUNCE BACK TO LOCAL??
def ABtrans(A,B):
  if A and B in ['1','2']:
    if A==B=='1': return '1'
    if A==B=='2': return '2'
    if (A,B)==('1','2'): return '3'
    if (A,B)==('2','1'): return '3'
  else: return "error"

for chrm in mark_dict:
#for chrm in ['1']:
        oufi=open('../phased/Raw_phased_chrm%s.txt'%(chrm),'w')
        print(item)
        snps=[mark_num_dict[x[0]] for x in mark_dict[chrm]]     
        hg=open("../phase_outputs/chrm%s_hapguess_switch.out"%(chrm)).readlines() #opens fastphase outputfile  
        assert len(snps)==len(hg[22].split())#makes sure number of snps in head == number in lines
        for x,lin in enumerate(hg):#skipps all the headers and footers
           if lin.startswith('# id'):
                idi=lin.split()[2]
                for i,snp in enumerate(snps):
                    A=hg[x+1].split()[i]
                    B=hg[x+2].split()[i]
                    geno=ABtrans(A,B)
                    oufi.write("["+snp+','+idi+','+geno+']'+'\n')
           else:pass             



##run PCA


def raw_to_eig(infiles,outstr,group='region'):#infiles is a list of files    
    def trans(x):
        if x=='1':return '0'
        if x=='2':return '2'
        if x=='3':return '1'
        if x=='10':return '10'
    def ind_info(ind):
        ind_info=open('ids.csv','r').readlines()
        inddict={}
        for lin in ind_info:
            lii=lin.split(',')
            inddict[lii[0]]=lii[1:]
        return inddict[ind]
    if type(infiles) != list:
       print("infiles needs to be alist")
    snps=set()
    ids=set()
    snpfi=open("../eig/"+outstr+'.snp','w')
    goufi=open("../eig/"+outstr+'.geno','w')
    indfi=open("../eig/"+outstr+'.ind','w')
    for infile in infiles:
        goufi=open("../eig/"+outstr+'.geno','a')
        indfi=open("../eig/"+outstr+'.ind','a')
        print(infile)
        fi=open(infile)
        rawdat=[]
        for lin in fi:
          try:
                stri=str(str(lin)[1:-2]) #strips off brackets
                stri2=(str(stri)).split(',')
                snpnum=stri2[0].strip()
                snpnam=mark_nam_dict[snpnum]
                lst=[snpnum,stri2[1].strip(),trans(stri2[2].strip())]
                ids.add(stri2[1].strip())
                snps.add(snpnum)
                if (lst[2]!='10'):
                                rawdat.append(lst)
          except:pass
        for lin in rawdat:
                goufi.write(lin[0]+' '+lin[1]+' '+lin[2]+'\n')
        goufi.close()
    for item in ids:
        if group=='region':
               indfi.write(item+'\tF\t'+ind_info(item)[4].strip('"')+'\n')#by region      
        if group =='breed':
               indfi.write(item+'\tF\t'+ind_info(item)[2].strip('"')+'\n')#by breedcode
        if group =='country':
               indfi.write(item+'\tF\t'+ind_info(item)[5].strip('"')+'\n')
    indfi.close()
    for snp in snps:
            nam=mark_nam_dict[snp] 
            (num,chrm, pos)=mark_chrm_dict[nam]
            snpfi.write(" ".join([snp,chrm,'0.0',str(pos)])+'\n') #pulls outsnp number, chrm, dummy recombination distance, lcoaction
    snpfi.close()
    par=open("../eig/par.example",'r')
    opar=open('../eig/par.'+outstr,'w')
    for lin in par.readlines():
           opar.write(lin.replace('example',outstr))
    opar.close()
    runner="../../EIG4.2/bin/smartpca -p ../eig/par.%s"%outstr
    print(runner)
    os.system(runner)
    os.system("../../EIG4.2/bin/ploteig -p ../eig/par.%s"%outstr)

infs=['../phased/Raw_phased_chrm%s.txt'%(chrm) for chrm in ['1']]
fis=glob.glob("../phase_outputs/*.txt")
##estimate indivdual ancestry


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



##------------------------------------------------------------------------------------------
##SECTION OF THINGS THAT HAPPEN ONLY TO 3k DATA
three=open("../autosomes3k.csv").readlines()
snpnamesTK=set()
for item in autos:
  lii=item.split(',')
  if lii[1] not in remove:
    snpnamesTK.add(lii[1])

##Read in list of indivduals
genosTK={}
for ind in inds:
  genosTK[ind]={}

# <codecell>

##Run through raw data files for individuals to test missingness:
#fis=glob.glob("../RawData_PNAS/*.txt")
fis=["../RawData_PNAS/genos_116.txt"]
for fi in fis:
  fil=open(fi)
  for lin in fil:
      lii=lin[1:-2].split(',')
      snpnam=ibmclocdict[lii[0].strip()]
      if snpnam in snpnamesTK:
          genosTK[lii[1].strip()][snpnam]=lii[2].strip()
    


# <codecell>

##Drop indivs missing >10%
miss_indsTK=[]
for ki in genosTK:
  if len(genosTK[ki]) < len(snpnamesTK) * 0.9:
    miss_indsTK.append(ki)

for ind in miss_inds:
    del genosTK[ind]


##define way to make structure files from 3k data
stfile=open("Structureout.str",'w')
positions=[0]
snpnums=[]
for chrm in mark_dict:
	for i in range(len(mark_dict[chrm])-1):
            try:
                dist=mark_dict[chrm][i+1][1]-mark_dict[chrm][i][1]#calc distance between markers
                if dist ==0:
                	dist = -1            
            except:
                dist=-1
            header.append(dist)
            mark_dict[chrm][i+1].add(dist) ##HOW TO INSERT INTO TUPLE??!?
		

## FIIIIXXX MEEEEEE
def raw_to_struct(finame):
    (ibmclocdict, locdict, nam_dict,trans,chrm_dict,chrm_n,forw_dict,rev_dict,inddict,ABdict)=make_dicts()
    print(finame)
    fi=open(finame)
    rawdat=[]
    i=1
    for lin in fi:
            stri=str(str(lin)[1:-2])
            stri2=(str(stri)).split(',')
            try:
                #renames snpnum from ibmc numbers to SNP_Map numbers in rawdat
                snpnum=nam_dict[ibmclocdict[stri2[0].strip()]] 
                lst=[snpnum,stri2[1].strip(),stri2[2].strip()]
                rawdat.append(lst)                                                    
            except:
                i=i+1
    print("Skipped lines = %i" %i)
    fi.close()
    print('data read')
    nudat=rawdat
    #now nudat holds all the raw data recoded as [SNP_Map snpnumber,individual id number, genotype]
    snp_dict={}
    indiv_n=set()
    for item in nudat:#dictionarizes nudat
        indiv_n.add(item[1])#makes a lis of inidvdual ids
        snp_dict[str(item[0]),item[1]]=item[2]#tuple key of (SNP_Map #, ID #) calls genotype
    
    indiv_n=list(indiv_n)
    indiv_n.sort()
    for chrm in chrm_n: #goes through each chromosome SNP dictionary, sorts it, and then pulls the genotypes for that SNP.
        #here globals()[nm] calls the list called chrm1 then sorts it
        #return locals()[nm]
        chrm_dict[chrm]=sorted(chrm_dict[chrm], key=operator.itemgetter(3))
        chrm_dict[chrm]=sorted(chrm_dict[chrm], key=operator.itemgetter(2))
        header=[-1]
        print("sorted")
        for i in range(len(chrm_dict[chrm])-1):
            try:
                dist=chrm_dict[chrm][i+1][3]-chrm_dict[chrm_n][i][3]#calc distance between markers
            except:
                dist=-1
            header.append(dist)
            chrm_dict[chrm][i+1].append(dist)
        matri=[]
        for ite in indiv_n: #starts each line with id number
            matri.append([ite])
        snp_order=[]
        for snp in chrm_dict[chrm]:#lists snps in order
            snp_order.append(str(snp[0]))
        for ids in matri:
            idss=ids[0]#    
            for snp in snp_order:
                try:
                    ids.append(snp_dict[(snp,idss)])
                except KeyError: 
                    ids.append('10')
        newdat=matri
        allele1=copy.deepcopy(newdat[:])
        allele2=copy.deepcopy(newdat[:])
        for x in range(len(newdat)):
            for y in range(len(newdat[x])):
                if newdat[x][y]=='1':
                    allele1[x][y]='1'
                    allele2[x][y]='1'
                elif newdat[x][y]=='2':
                    allele1[x][y]='2'
                    allele2[x][y]='2'
                elif newdat[x][y]=='3':
                    allele1[x][y]='1'
                    allele2[x][y]='2'
                elif newdat[x][y]=='10':
                    allele1[x][y]='-9'
                    allele2[x][y]='-9'
                elif newdat[x][y]=='':
                    allele1[x][y]='-9'
                    allele2[x][y]='-9'
        structdat1=[]
        structdat1=allele1+allele2
        structdat=copy.deepcopy(structdat1)
        structdat.sort()
        for x in range(len(structdat)):
            structdat[x].insert(1,structdat[x][0][0:3])
        for x in range(len(structdat)):
            if (int(structdat[x][0])>121488000 and int(structdat[x][0])<121488999):
                structdat[x][1]='999'
        for x in range(len(structdat)):
            if (int(structdat[x][0])>=999999000 and int(structdat[x][0])<999999100):
                structdat[x][1]='998'
        outputfile='%s.txt1'%chrm
        if os.path.exists(outputfile): pass
        else:
          structdat.insert(0,copy.deepcopy(snp_order))
          structdat.insert(1, header)        
        fs=open(outputfile, 'a')
        for x in range(len(structdat)):
            for y in range(len(structdat[x])):
                fs.write(str(structdat[x][y]))
                fs.write(' ')
            fs.write('\n')
        fs.close()




#actually makes the structure files

#for infile in glob.glob( os.path.join(path, '*.txt') ):
#    print infile
#    raw_to_struct(infile)



