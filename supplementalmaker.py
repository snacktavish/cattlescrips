
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
    mark_chrm_dict[nam]=(num,chrm, str(pos))
    mark_num_dict[nam]=num #translates to snp number
    mark_nam_dict[num]=nam #translates back to SNP name

    
f3k=open("../LGH_struct_ms/eig/3k.snp",'r').readlines()
f50k=open("../LGH_struct_ms/eig/50K.snp",'r').readlines()

snp3=open("50ksnps.csv",'w')
snp50=open("18snps.csv",'w')

def csvify(infi,oufi):
  oufi.write("Snpnum,Name,Numcheck,chrm,pos\n")
  for lin in infi:
    snp=[lin.split()[0]]
    snp.append(mark_nam_dict[snp[0]])
    out=snp+list(mark_chrm_dict[snp[1]])
    oufi.write(','.join(out)+'\n')
  oufi.close()

