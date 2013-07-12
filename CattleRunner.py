# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import SnpParing
import pickle
from CattleDictionaries import cownames, mark_dict
import MakeGenos
import formats


#Load in already processed FULL GENOS FILE.
genos=pickle.load(open("genos_span.pkl","rb"))
print("Loaded pickfile")


#MakeGenos.get_cownames(genos,["china/Bovine GGP HD 28mar2013/Penn State-Lei GGP HD 28mar2013_FinalReport.txt"])

genoB=MakeGenos.add_final_reports(genos,["china/Bovine GGP HD 28mar2013/Penn State-Lei GGP HD 28mar2013_FinalReport.txt"])
print("geno B len = %i" (len(genoB)))
pickB=open("genoCHIN",'wb')
pickle.dump(genoB,pickB)


prefstrs=[str(val) for val in range(901000000,902001000)]

genoC,filtered_snps=SnpParing.snp_strip(genoB, cutoff=0.3, preflist=prefstrs)
print("geno C len = %i" %len(genoC))


genoD=SnpParing.ind_strip(genoC,cutoff=0.8,len(filtered_snps))
print("geno D len = %i" %len(genoD))

genoE,filt_snpsE=SnpParing.snp_strip(genoD, cutoff=0.8, preflist=None)
print("geno E len = %i" %len(genoE))
pickE=open("genoCHIN_strip388",'wb')
pickle.dump(genoE,pickE)


#    runner="../EIG4.2/bin/smartpca -p eigSPAN/par.%s"%outstr
#    print(runner)
#    os.system(runner)
#    os.system("perl ../EIG4.2/bin/ploteig -i eigSPAN/%s.evec -p ../EIG4.2/POPGEN/breeds.txt  -x -o %s.xtxt"%(outstr,outstr))

# <codecell>



formats.phase_inp(genoE,mark_dict,filt_snpsE,"china/phase_inp/chin")

#PHASE IT ESLSEWHWE!
#scp -r china/phase_inp/ ebm447@phylocluster.ccbb.utexas.edu:/home/ebm447/china/

