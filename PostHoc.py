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
fil=open("./PCA/ibmc_markers_041709.csv", 'r')
for lin in fil:
        row=lin.split(',')
        num=row[1]
        nam=row[9]
        ibmclocdict[num]=nam


fil.close()

# <codecell>

# seperate out markers by chormosome
markers=open("UMC_snp50_markers_121007.csv").readlines()
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

# <codecell>

#this should make it so that there is a dictionary where the key is 1-30, and the values are a list of names and positions in order
for item in mark_dict.keys():
      mark_dict[item].sort(key=operator.itemgetter(1))


# <codecell>

##Read in marker files
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

##Read in list of indivduals
indfi=open("ids_full_SPAN.csv").readlines()
cownames={}
inds=[]
for lin in indfi[1:]:
  inds.append(lin.split(",")[0])
  try: 
     cownames[lin.split(',')[7].strip()]=lin.split(",")[0]
  except: pass

del cownames['x']


# <codecell>

genos={}
for ind in inds:
  genos[ind]={}

# <codecell>

##Run through raw data files for individuals to test missingness:
fis=glob.glob("RawData_PNAS/*.txt")
##fis=["../RawData_PNAS/genos_101.txt","../RawData_PNAS/genos_503.txt"]
for fi in fis:
  fil=open(fi)
  for lin in fil:
      lii=lin[1:-2].split(',')
      snpnam=ibmclocdict[lii[0].strip()]
      if snpnam in snpnames and lii[2].strip() != '10':
          genos[lii[1].strip()][snpnam]=lii[2].strip()




# <codecell>

def gen(X,Y): #Decode's AB in genotypes
    if X=='A' and Y =='A': return '1'
    if X=='B' and Y =='B': return '2'
    if X=='A' and Y =='B': return '3'
    if X=='B' and Y =='A': return '3'

# <codecell>

for item in ['504202950','504202970','5042023190','514203200','504203220','504203230']:#problematic nelore individuals
    try:
        del genos[item]
    except: pass

# <codecell>

new_names=set()
fis=glob.glob("final_reports/*.txt")
for fi in fis:
  fil=open(fi).readlines()
  for lin in fil[10:]:
      lii=lin.split("\t")
      try:
            cownames[lii[1]]
      except:
                    new_names.add(lii[1]) 
            

# <codecell>

print(new_names)

# <codecell>

for item in new_names:
    print(item)

# <codecell>

##########################
#add in all the final report genos.
fis=glob.glob("final_reports/*.txt")
for fi in fis:
  fil=open(fi).readlines()
  assert(fil[9].split("\t")[0]=="SNP Name")
  for lin in fil[10:]:
      lii=lin.split("\t")
      if lii[0] in snpnames and float(lii[-1])>0.2:
        try:
           genos[cownames[lii[1]]][lii[0]]=gen(lii[6],lii[7])
        except:
            pass

# <codecell>

fi

# <codecell>

float(lii[-1])

# <codecell>

len(genos['888000001'])

# <codecell>

fi="final_reports/Texas Longhorn 21oct2009_FinalReport AB.txt2"
fil=open(fi).readlines()
assert(fil[9].split("\t")[0]=="SNP Name")
for lin in fil[10:]:
      lii=lin.split("\t")
      if lii[0] in snpnames:
        try:
           genos[cownames[lii[1]]][lii[0]]=gen(lii[2],lii[3])
        except:
            pass

# <codecell>

pifi=open("genos_span.pkl","wb")
pickle.dump(genos,pifi)

# <codecell>

genos3k=pickle.load(open("genos_span.pkl","rb"))

# <codecell>

#genos3k=genos

# <codecell>

##STAGE 1 SNP REMOVAL
##Run through SNPs. Narrows them down to ones common to the SPANISH samples.
sub_genos={}
for item in genos3k:
    if item.startswith('888000'):
        sub_genos[item]=genos3k[item]

# <codecell>

len(sub_genos)

# <codecell>

len(snpnames)

# <codecell>

snp_miss={}
for snp in snpnames:
    snp_miss[snp]=0

for item in sub_genos:
      for ite in sub_genos[item]:
        snp_miss[ite]=snp_miss[ite]+1 #counts up how many indivs have misisng that snp

snp_removal=set()
for snp in snp_miss:
    if snp_miss[snp]<(len(sub_genos))*0.70:
        snp_removal.add(snp)

filtered_snps=snpnames-snp_removal

print("number of snps is %i" %len(filtered_snps))


# <codecell>

#STRIP non 19K SNPS FROM GENOS
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

# <codecell>

def keep(x):
    return x[0] in final_snps

for chrm in mark_dict:
      mark_dict[chrm]= [x for x in mark_dict[chrm] if keep(x)] #strips snpa with too much missing data out of marker dictionary

del mark_dict['']

del mark_dict['umd30_bta']

del mark_dict['31']

del mark_dict['99']

del mark_dict['30']

# <codecell>

#pifi=open("mark_dictSPAN.pkl","wb")
#pickle.dump(mark_dict,pifi)

# <codecell>

mark_dict=pickle.load(open("mark_dictSPAN.pkl","rb"))

# <codecell>

len(mark_dict)

# <codecell>

def phase_inp(chrm):
        print("Hromosome %s" %chrm)
        oufi=open("phase_inputSPAN/phase_chrm%s.inp"%chrm,'w')
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

# <codecell>

for chrm in mark_dict:
	phase_inp(chrm)

# <codecell>

for chrm in mark_dict:
#for chrm in ['1']:
    os.system("./fastPHASE_MacOSX-Darwin -ophase_outputSPAN/chrm%s phase_inputSPAN/phase_chrm%s.inp"%(chrm,chrm))##RECHECK PHASING SETTINGS

# <codecell>

def ABtrans(A,B):
  if A and B in ['1','2']:
    if A==B=='1': return '1'
    if A==B=='2': return '2'
    if (A,B)==('1','2'): return '3'
    if (A,B)==('2','1'): return '3'
  else: return "error"

# <codecell>

for chrm in mark_dict:
#for chrm in ['1','24','25','26','27','20','21','22','23','28']:
        oufi=open('phasedSPAN/Raw_phased_chrm%s.txt'%(chrm),'w')
        snps=[mark_num_dict[x[0]] for x in mark_dict[chrm]] #gets the number for each snp from the name in mark dict
        hg=open("phase_outputSPAN/chrm%s_hapguess_switch.out"%(chrm)).readlines() #opens fastphase outputfile  
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

# <codecell>

len(snps)

# <codecell>

'''##run PCA
ind_skip=[]
for lin in open("ind_rem.txt"):
    ind_skip.append(lin.strip())

ind_info=open('ids_full_miss.csv','r').readlines()
inddict={}
for lin in ind_info:
            lii=lin.split(',')
            inddict[lii[0]]=lii[1:]
            
            
for item in inddict:
    if inddict[item][4]=="'Africa'":
       ind_skip.append(item)'''

# <codecell>

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
        ind_info=open('ids_fullSPAN.csv','r').readlines()
        inddict={}
        for lin in ind_info:
            lii=lin.split(',')
            inddict[lii[0]]=lii[1:]
        return inddict[ind]
    if type(infiles) != list:
       print("infiles needs to be alist")
    snps=set()
    ids=set()
    snpfi=open("eigSPAN/"+outstr+'.snp','w')
    goufi=open("eigSPAN/"+outstr+'.geno','w')
    indfi=open("eigSPAN/"+outstr+'.ind','w')
    for infile in infiles:
        goufi=open("eigSPAN/"+outstr+'.geno','a')
        indfi=open("eigSPAN/"+outstr+'.ind','a')
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
    par=open("eig/par.example",'r')
    opar=open('eigSPAN/par.'+outstr,'w')
    for lin in par.readlines():
           opar.write(lin.replace('example',outstr))
    opar.close()
#    runner="../EIG4.2/bin/smartpca -p eigSPAN/par.%s"%outstr
#    print(runner)
#    os.system(runner)
#    os.system("perl ../EIG4.2/bin/ploteig -i eigSPAN/%s.evec -p ../EIG4.2/POPGEN/breeds.txt  -x -o %s.xtxt"%(outstr,outstr))

# <codecell>

def trans(x):
        if x=='1':return '0'
        if x=='2':return '2'
        if x=='3':return '1'
        if x=='10':return '10'

ind_info=open('ids_full_SPAN.csv','r').readlines()
inddict={}
for lin in ind_info:
            lii=lin.split(',')
            inddict[lii[0]]=lii[1:]

# <codecell>

def raw_to_eig(infiles,outstr,group='region'):#infiles is a list of files    
    try:
        len(mark_nam_dict)
    except:
      print("Create marker dictionary!")
    try:
        len(inddict)
    except:
      print("Create ind info dictionary!")
    if type(infiles) != list:
       print("infiles needs to be alist")
    snps=set()
    ids=set()
    snpfi=open("eigSPAN/"+outstr+'.snp','w')
    goufi=open("eigSPAN/"+outstr+'.geno','w')
    indfi=open("eigSPAN/"+outstr+'.ind','w')
    for infile in infiles:
        goufi=open("eigSPAN/"+outstr+'.geno','a')
        indfi=open("eigSPAN/"+outstr+'.ind','a')
        print(infile)
        fi=open(infile)
        rawdat=[]
        for lin in fi:
                stri=str(str(lin)[1:-2]) #strips off brackets
                stri2=(str(stri)).split(',')
                snpnum=stri2[0].strip()
                snpnam=mark_nam_dict[snpnum]
                lst=[snpnum,stri2[1].strip(),trans(stri2[2].strip())]
                ids.add(stri2[1].strip())
                snps.add(snpnum)
                if (lst[2]!='10'):
                                rawdat.append(lst)
        for lin in rawdat:
                goufi.write(lin[0]+' '+lin[1]+' '+lin[2]+'\n')
        goufi.close()
    for item in ids:
        if group=='region':
               indfi.write(item+'\tF\t'+inddict[item][4].strip('"')+'\n')#by region      
        if group =='breed':
               indfi.write(item+'\tF\t'+inddict[item][2].strip('"')+'\n')#by breedcode
        if group =='country':
               indfi.write(item+'\tF\t'+inddict[item][5].strip('"')+'\n')
    indfi.close()
    for snp in snps:
            nam=mark_nam_dict[snp] 
            (num,chrm, pos)=mark_chrm_dict[nam]
            snpfi.write(" ".join([snp,chrm,'0.0',str(pos)])+'\n') #pulls outsnp number, chrm, dummy recombination distance, lcoaction
    snpfi.close()
    par=open("eig/par.example",'r')
    opar=open('eigSPAN/par.'+outstr,'w')
    for lin in par.readlines():
           opar.write(lin.replace('example',outstr))
    opar.close()

# <codecell>

outstr="test"
snps=set()
ids=set()
snpfi=open("eigSPAN/"+outstr+'.snp','w')
goufi=open("eigSPAN/"+outstr+'.geno','w')
indfi=open("eigSPAN/"+outstr+'.ind','w')
infile=fis[0]
goufi=open("eigSPAN/"+outstr+'.geno','a')
indfi=open("eigSPAN/"+outstr+'.ind','a')
print(infile)
fi=open(infile)
rawdat=[]
for lin in fi:
#          try:
                stri=str(str(lin)[1:-2]) #strips off brackets
                stri2=(str(stri)).split(',')
                snpnum=stri2[0].strip()
                snpnam=mark_nam_dict[snpnum]
                lst=[snpnum,stri2[1].strip(),trans(stri2[2].strip())]
#                if stri2[1].strip() not in ind_skip:
                ids.add(stri2[1].strip())
                snps.add(snpnum)
                if (lst[2]!='10'):
                                rawdat.append(lst)
#          except:
#                print(lin)#
#                 pass

# <codecell>

for lin in rawdat:
                goufi.write(lin[0]+' '+lin[1]+' '+lin[2]+'\n')

goufi.close()

# <codecell>

ind_skip

# <codecell>


##infs=['../phased/Raw_phased_chrm%s.txt'%(chrm) for chrm in ['1']]
fis=glob.glob("phasedSPAN/*.txt")
##estimate indivdual ancestry

raw_to_eig(fis,"regions", group='region')

#raw_to_eig(fis,"BREEDS_F",group='breed')



# <codecell>

-i DeprecationWarningsnps=set()
for infi in fis:
   fi=open(infi).readlines()
   for lin in fi:
#          try:
                stri=str(str(lin)[1:-2]) #strips off brackets
                stri2=(str(stri)).split(',')
                snpnum=stri2[0].strip()
                snps.add(snpnum)


                snpnam=mark_nam_dict[snpnum]

# <codecell>

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

# <codecell>

perc_dict=calc_perc("eigSPAN/BREEDS_F")

# <codecell>

for item in perc_dict.keys():
    if item.startswith('88'):
        print(item)

# <codecell>

perc_dict['888000025']

# <codecell>


