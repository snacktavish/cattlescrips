
fis=glob.glob("../final_reports/*.txt*")
newids=set()
for fi in fis:
  fil=open(fi).readlines()
  assert(fil[9].split("\t")[0]=="SNP Name")
  for lin in fil[10:]:
      lii=lin.split("\t")
      newids.add(lii[1])

def idnum(Z):
	dicto={}
	for i,item in enumerate(newids):
		dicto[item]=str(i+999111000)
	return(dicto[Z])

for item in list(newids):
	print(item+", "+idnum(item))


for item in newids:
	genos[idnum(item)]={}


indfi=open("../ids_full.csv").readlines()
oufi=open("../ids_full_miss.csv","w")
oufi.write(indfi[0].strip()+", in_55k\n")
for lin in indfi[1:]:
  ind=(lin.split(",")[0])
  if ind in miss_inds:
    oufi.write(lin.strip()+", N\n")   
  else:
    oufi.write(lin.strip()+", Y\n")   

