# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from CattleDictionaries import cownames, snpnames, ibmclocdict, gen


#This doesn't test to see if indiviuals are in genos yet...
def add_final_reports(genos,fis,cutoff=0.7):
    ""Witre a reasonable docstring!
    ""
    new_genos=genos.copy()
    assert type(fis) is list
    skips=set()
    for fi in fis:
        fil=open(fi).readlines()
        assert(fil[9].split("\t")[0]=="SNP Name")
        assert (fil[9].split("\t")[-3]=="GC Score")
    for lin in fil[10:]:
        lii=lin.split("\t")
        if lii[1] in cownames:
            cname=cownames[lii[1]]
            if cname in new_genos:
               pass
            else:
                new_genos[cname]={}
            if lii[0] in snpnames and float(lii[-3])>cutoff:
                  new_genos[cname][lii[0]]=gen(lii[6],lii[7])
        else:
           skips.add(lii[1])
    print(skips)
    return(new_genos)

# <codecell>

def create_genos(raw_data_folder="RawData_PNAS/*.txt"):
    genos={}
    for ind in inds:
        genos[ind]={}
    ##Run through raw data files for individuals to test missingness:
    fis=glob.glob(raw_data_folder)
    for fi in fis:
         fil=open(fi)
    for lin in fil:
        lii=lin[1:-2].split(',')
        snpnam=ibmclocdict[lii[0].strip()]
        if snpnam in snpnames and lii[2].strip() != '10':
            genos[lii[1].strip()][snpnam]=lii[2].strip()
    #remove mislabeled Nelore individuals
    for item in ['504202950','504202970','5042023190','514203200','504203220','504203230']:#problematic nelore individuals
        try:
            del genos[item]
        except: pass
    #add in all the final report genos.
    add_final_reports(genos)
    return(genos)

# <codecell>

def get_cownames(genos,fis):
    skips=set()
    for fi in fis:
        fil=open(fi).readlines()
        assert(fil[9].split("\t")[0]=="SNP Name")
    for lin in fil[10:]:
        lii=lin.split("\t")
        if lii[1] in cownames: pass
        else:
           skips.add(lii[1])
    return(skips)
    


