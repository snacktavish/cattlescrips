# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import SnpParing
import pickle
from CattleDictionaries import cownames
import MakeGenos


#Load in already processed FULL GENOS FILE.
genos=pickle.load(open("genos_span.pkl","rb"))
print("Loaded pickfile")


#MakeGenos.get_cownames(genos,["china/Bovine GGP HD 28mar2013/Penn State-Lei GGP HD 28mar2013_FinalReport.txt"])

genoB=MakeGenos.add_final_reports(genos,["china/Bovine GGP HD 28mar2013/Penn State-Lei GGP HD 28mar2013_FinalReport.txt"])
print("geno B len = %i" (len(genoB)))
pickB=open("genoCHIN",rb)
pickle.dump(genoB,pickB)

prefstrs=[str(val) for val in range(901000000,902000000)]

genoC=SnpParing.snp_strip(genoB, cutoff=0.3, preflist=prefstrs)
print("geno C len = %i" (len(genoC)))


genoD=SnpParing.ind_strip(genoC,cutoff=0.8)
print("geno D len = %i" (len(genoD)))

genoE=SnpParing.snp_strip(genoB, cutoff=0.8, preflist=None)
print("geno E len = %i" (len(genoE)))
pickE=open("genoCHIN_strip388",rb)
pickle.dump(genoE,pickE)

pickout=open(stripped_chin)