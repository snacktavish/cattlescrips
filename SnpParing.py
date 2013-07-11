# -*- coding: utf-8 -*-


from CattleDictionaries import snpnames


def snp_strip(genos, cutoff, preflist=None):
    #takes a genos base, strips SNPS with higher than %cutoff missingness.
    snp_miss={}
    for snp in snpnames:
        snp_miss[snp]=0
    sub_genos={}
    if preflist: 
      for item in genos:
         if item in preflist:
             sub_genos[item]=genos[item]
    else:
        sub_genos=genos.copy()
    for item in sub_genos:
      for ite in sub_genos[item]:
        snp_miss[ite]=snp_miss[ite]+1 #counts up how many indivs have that snp
    snp_removal=set()
    for snp in snp_miss:
         if snp_miss[snp]<(len(sub_genos))*cutoff:
            snp_removal.add(snp)
            filtered_snps=snpnames-snp_removal
    strip_genos={}
    for ind in genos:
        strip_genos[ind]={}
        for snp in filtered_snps:
            try:
              strip_genos[ind][snp]=genos[ind][snp]
            except: pass
    return(strip_genos,filtered_snps)
    print("number of snps is %i" %len(filtered_snps))



def ind_strip(genos, cutoff,num_snps):
    #returns a dictionray w/o indivs w more than %cutoff missingness
    keep_inds=[]
    strip_genos={}
    for ind in genos:
       if len(genos[ind]) > num_snps * cutoff:
            keep_inds.append(ind)
    for ind in keep_inds:
        strip_genos[ind]=genos[ind]
    return(strip_genos)
            


