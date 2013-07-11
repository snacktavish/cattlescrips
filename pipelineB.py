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

for item in ['504202950','504202970','5042023190','514203200','504203220','504203230']:#problematic nelore individuals
	del genos[item]

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



pifi=open("genos_5.pkl","wb")
pickle.dump(genos,pifi)

'''
indsnum=[]
for ki in genos:
  indsnum.append((len(genos[ki]),ki))

TK=set()
TK=list(TK)
TK.sort()
for item in indsnum[:75]:
  TK.add(item[1])
'''


##STAGE 1 SNP REMOVAL
##Run through SNPs. Remove any missing in >70 of inds%
snp_miss={}
for snp in snpnames:
    snp_miss[snp]=0

for item in genos:
      for ite in genos[item]:
        snp_miss[ite]=snp_miss[ite]+1 #counts up how many indivs have misisng that snp

snp_removal=set()
for snp in snp_miss:
    if snp_miss[snp]<(len(genos))*0.3:
        snp_removal.add(snp)

filtered_snps=snpnames-snp_removal

print("number of snps is %i" %len(filtered_snps))

#STRIP BAD SNPS FROM GENOS
for ki in genos:
   for mark in snp_removal:
     try:
      del genos[ki][mark]
     except: pass


##NOW Drop indivs missing >10%
miss_inds=[]
for ki in genos:
  if len(genos[ki]) < len(filtered_snps) * 0.9:
    miss_inds.append(ki)



for ind in miss_inds:
    del genos[ind]


##STAGE 2 SNP REMOVAL
##Run through SNPs AGAIN. Remove any missing >10%
snp_miss2={}
for snp in filtered_snps:
    snp_miss2[snp]=0

for item in genos:
      for ite in genos[item]:
        snp_miss2[ite]=snp_miss2[ite]+1 #counts up how many indivs have that snp


snp_removal2=set()
for snp in snp_miss2:
    if snp_miss2[snp]<len(genos)*0.9:
        snp_removal2.add(snp)

final_snps=filtered_snps-snp_removal2

print("number of snps is %i" %len(final_snps))

##### REMOVE THEM FORM THE mark_dict?!
##****FIND A BETTER WAY TO REMOVE ITEMS FROM DICT
def keep(x):
    return x[0] in final_snps

for chrm in mark_dict:
      mark_dict[chrm]= [x for x in mark_dict[chrm] if keep(x)] #strips snpa with too much missing data out of marker dictionary

del mark_dict['']

del mark_dict['umd30_bta']

del mark_dict['31']

del mark_dict['99']

pifi=open("../phase_inputs/mark_dict.pkl","wb")
pickle.dump(mark_dict,pifi)
##make phase input files

def phase_inp(chrm):
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
        oufi.close()

for chrm in mark_dict:
	phase_inp(chrm)


##Actually phase the data
##MAKE THIS HAPPEN ON THE CLUSTER

for chrm in mark_dict:
    os.system("/home/ejbm/Documents/Phase/fastPHASE/fastPHASE_Linux -o../phase_outputs/chrm%s ../phase_inputs/phase_chrm%s.inp"%(chrm,chrm))##RECHECK PHASING SETTINGS

##for chrm in ['2']:
##    os.system("/home/ejbm/Documents/Phase/fastPHASE/fastPHASE_Linux -T2 -o../phase_outputs/chrm%s ../phase_inputs/phase_chrm%s.inp"%(chrm,chrm))##RECHECK PHASING SETTINGS


##BOUNCE BACK TO LOCAL??
##Translate phased genos to EIG
def ABtrans(A,B):
  if A and B in ['1','2']:
    if A==B=='1': return '1'
    if A==B=='2': return '2'
    if (A,B)==('1','2'): return '3'
    if (A,B)==('2','1'): return '3'
  else: return "error"


##Writes phased data back down to chrms
for chrm in mark_dict:
for chrm in ['25']:
        oufi=open('../phased/Raw_phased_chrm%s.txt'%(chrm),'w')
        snps=[mark_num_dict[x[0]] for x in mark_dict[chrm]] #gets the number for each snp from the name in mark dict
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




for chrm in mark_dict:
     hg=open("../phase_outputs/chrm%s_hapguess_switch.out"%(chrm)).readlines()


##run PCA

ind_skip=[]
for lin in open("../ind_rem.txt"):
    ind_skip.append(lin.strip())

ind_info=open('../ids_full_miss.csv','r').readlines()
inddict={}
for lin in ind_info:
            lii=lin.split(',')
            inddict[lii[0]]=lii[1:]



for item in inddict:
    if inddict[item][4]=="'Africa'":
       ind_skip.append(item)

def raw_to_eig(infiles,outstr,group='region'):#infiles is a list of files    
    try:
        len(mark_nam_dict)
    except:
      print("Create marker dictionary!")
    def trans(x):
        if x=='1':return '0'
        if x=='2':return '2'
        if x=='3':return '1'
        if x=='10':return '10'
    def ind_info(ind):
        ind_info=open('../ids_full_miss.csv','r').readlines()
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
                if stri2[1].strip() not in ind_skip:
                  ids.add(stri2[1].strip())
                  snps.add(snpnum)
                  if (lst[2]!='10'):
                                rawdat.append(lst)
          except:
                print(lin)
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
    os.system("perl ../../EIG4.2/bin/ploteig -i ../eig/%s.evec -p ../../EIG4.2/POPGEN/breeds.txt  -x -o %s.xtxt"%(outstr,outstr))


##infs=['../phased/Raw_phased_chrm%s.txt'%(chrm) for chrm in ['1']]
fis=glob.glob("../phased/*.txt")
##estimate indivdual ancestry

raw_to_eig(fis,"sub_ind_breed_F", group='breed')

raw_to_eig(fis,"BREEDS_F",group='breed')


snps=set()
for infi in fis:
   fi=open(infi).readlines()
   for lin in fi:
#          try:
                stri=str(str(lin)[1:-2]) #strips off brackets
                stri2=(str(stri)).split(',')
                snpnum=stri2[0].strip()
                snps.add(snpnum)


                snpnam=mark_nam_dict[snpnum]

def calc_perc(instr):
  fi=open(instr+".evec").readlines()
  tau=[]
  ind=[]
  for lin in fi:
    if lin.split()[0].startswith(tuple(['100', '101', '124', '129'])):
      tau.append(float(lin.split()[1]))
    if lin.split()[0].startswith(tuple(['506','507','505','503','504'])):
      ind.append(float(lin.split()[1]))
  tavg=sum(tau)/len(tau)
  iavg=sum(ind)/len(ind)
  perc_dict={}
  for lin in fi:
    perc_dict[lin.split()[0]]=(float(lin.split()[1])-iavg)/(tavg-iavg)
  return perc_dict




##------------------------------------------------------------------------------------------
##For 3k:

genos3k=pickle.load(open("genos_5.pkl","rb"))


##STAGE 1 SNP REMOVAL
##Run through SNPs. Remove any missing in >70 of inds%
genos={}
for item in genos3k:
    if item.startswith('999111'):
        genos[item]=genos3k[item]

snp_miss={}
for snp in snpnames:
    snp_miss[snp]=0

for item in genos:
      for ite in genos[item]:
        snp_miss[ite]=snp_miss[ite]+1 #counts up how many indivs have misisng that snp

snp_removal=set()
for snp in snp_miss:
    if snp_miss[snp]<(len(genos))*0.70:
        snp_removal.add(snp)

filtered_snps=snpnames-snp_removal

print("number of snps is %i" %len(filtered_snps))

#STRIP BAD SNPS FROM GENOS
for ki in genos3k:
   for mark in snp_removal:
     try:
      del genos3k[ki][mark]
     except: pass


##NOW Drop indivs missing >10%
miss_inds=[]
for ki in genos3k:
  if len(genos3k[ki]) < len(filtered_snps) * 0.9:
    miss_inds.append(ki)


for ind in miss_inds:
    del genos3k[ind]


##STAGE 2 SNP REMOVAL
##Run through SNPs AGAIN. Remove any missing >10%
snp_miss2={}
for snp in filtered_snps:
    snp_miss2[snp]=0

for item in genos3k:
      for ite in genos3k[item]:
        snp_miss2[ite]=snp_miss2[ite]+1 #counts up how many indivs have that snp


snp_removal2=set()
for snp in snp_miss2:
    if snp_miss2[snp]<len(genos)*0.9:
        snp_removal2.add(snp)

final_snps=filtered_snps-snp_removal2

print("number of snps is %i" %len(final_snps))

##### REMOVE THEM FORM THE mark_dict?!
##****FIND A BETTER WAY TO REMOVE ITEMS FROM DICT
def keep(x):
    return x[0] in final_snps

for chrm in mark_dict:
      mark_dict[chrm]= [x for x in mark_dict[chrm] if keep(x)] #strips snpa with too much missing data out of marker dictionary

del mark_dict['']

del mark_dict['umd30_bta']

del mark_dict['31']

del mark_dict['99']

del mark_dict['30']


pifi=open("../phase_inputs3k/mark_dict.pkl","wb")
pickle.dump(mark_dict,pifi)
##make phase input files


def phase_inp(chrm):
        print("Hromosome %s" %chrm)
        oufi=open("../phase_outputs3k/phase_chrm%s.inp"%chrm,'w')
        oufi.write(str(len(genos3k))+"\n")
        oufi.write(str(len(mark_dict[chrm]))+"\n"+"P ")
        for mark in mark_dict[chrm]:
             oufi.write(str(mark[1])+' ')
        oufi.write('\n')
        oufi.write("S" * len(mark_dict[chrm])+"\n")#WTF ARE THESE SSSSs for?
        for ind in genos3k:
                oufi.write("# id %s"%ind+"\n")
                hap1=[]
                hap2=[]
                for mark in mark_dict[chrm]:
                    try:
                      gen=genos3k[ind][mark[0]]
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
        oufi.close()

for chrm in mark_dict:
	phase_inp(chrm)


redo=[]
for chrm in mark_dict:
#for chrm in ['16','19','18']:
  try:
        oufi=open('../phased3k/Raw_phased3k_chrm%s.txt'%(chrm),'w')
        snps=[mark_num_dict[x[0]] for x in mark_dict[chrm]] #gets the number for each snp from the name in mark dict
        hg=open("../phase_outputs3k/chrm3k_%s_hapguess_switch.out"%(chrm)).readlines() #opens fastphase outputfile  
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
  except(AssertionError):
       redo.append(chrm)



## FIIIIXXX MEEEEEE

popnum={}
for ind in genos3k:
    if 999111000<int(ind)<999111999:
       if ind in ['999111012','999111039','999111047','999111050','999111058']:
          popnum[ind]='998'
       else:
          popnum[ind]='121'
    else:
       popnum[ind]=ind[:3]


oufi=open("struct3k.inp","w")
nams = []
pos = []
for chrm in mark_dict:
        for i,mark in enumerate(mark_dict[chrm]):
            nams.append(mark_num_dict[mark[0]])
            if i == 0:
               pos.append('-1')
            else:
               pos.append(str(mark_dict[chrm][i][1]-mark_dict[chrm][i-1][1]))


oufi.write(" ".join(nams)+"\n")
oufi.write(" ".join(pos)+"\n")


for ind in genos3k:
        hap1=[ind,popnum[ind]]
        hap2=[ind,popnum[ind]]
        for chrm in mark_dict:
            for mark in mark_dict[chrm]:
                    try: gen=genos3k[ind][mark[0]]
                    except: gen='10'
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
                       hap1.append('-9')
                       hap2.append('-9')
                    if gen not in ['1','2','3','?','10']:
                       hap1.append('-9')
                       hap2.append('-9')
        oufi.write(" ".join(hap1)+"\n")
        oufi.write(" ".join(hap2)+"\n")   


oufi.close()

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



