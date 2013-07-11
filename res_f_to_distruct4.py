##res_f_to_distruct.py
import os, glob
import operator

##this path should be where all of your results files are
path = '/home/ejbm/Documents/LGH_struct_ms/struct_nov/'
start=len(path)
os.chdir(path)
print(path)
##this for loop iterates through all the files in that folder

sortindqs=[]
sortdat=open("struct_out_2a_f").readlines()
n=sortdat[(sortdat.index('Run parameters:\n')+1)].split()[0]
i=(sortdat.index('Inferred ancestry of individuals:\n')+1)
while i<((sortdat.index('Inferred ancestry of individuals:\n')+2)+int(n)):
		i=i+1
		if "999 :" in sortdat[i]:
			print(sortdat[i])
			lii=sortdat[i].split()
			lii[3]='121'
			sortindqs.append(lii)
		elif "504203230" not in sortdat[i] and "504202950" not in sortdat[i] and "121692020" not in sortdat[i]:	
			sortindqs.append(sortdat[i].split())


sortindqs=sortindqs[:-1]
sortindqs.sort(key=operator.itemgetter(5))
order={}
for i, lin in enumerate(sortindqs):
	order[lin[0]]=i


for infile in glob.glob( os.path.join(path, '*_f') ):
	os.chdir(path)
##the -6 in the next line assumes that your file names end in.txt_f.  It names each folder with the distruct files the name of your results files. Prepare for many folders!
	nam=infile.split('/')[-1][:-2]
	print(nam)
##	nam='Link_10K_2a'
##infile='/home/emilyjane/Documents/BigTest/Links/samp10K/runs/link10K_1_6.txt_f'
	fi=open(infile)
	rawdat=fi.readlines()
	K=rawdat[(rawdat.index('Run parameters:\n')+3)].split()[0]
	loci=rawdat[(rawdat.index('Run parameters:\n')+2)].split()[0]
	n=rawdat[(rawdat.index('Run parameters:\n')+1)].split()[0]
	pops=59 #set number of predefined populations
	##obviously this next line leads to your distruct folder
	os.chdir('/home/ejbm/Documents/distruct1.1')
	dirname = nam
	if not os.path.isdir("./" + dirname + "/"):
		os.mkdir("./" + dirname + "/")
	os.chdir("./" + dirname + "/")
	fss=open("popq.popq", 'w')
	ind=1
	for i, lin in enumerate(rawdat):
		try:
			if lin.split()[0]=='Given':
				ind=i
			else: pass
		except: pass
	i=(ind+2)
	while i<(ind+2+1+pops):
		if rawdat[i].startswith('504:'):
			fss.write(' '.join(rawdat[i].split()[:-1])+" 60\n")
		elif rawdat[i].startswith('999:'): pass
		else:
			fss.write(rawdat[i])
		i=i+1
	fss.close()
	fs=open("indivq.indivq", 'w')
	##ditto with the K*3 + 60. I think it should be K*3 + 48 + your number of predefined populations
	i=(rawdat.index('Inferred ancestry of individuals:\n')+1)
	##replace 563 with your number of individuals
	indqs=[]
	while i<((rawdat.index('Inferred ancestry of individuals:\n')+2)+int(n)):
		i=i+1
		if "999 :" in rawdat[i]:
			print(rawdat[i])
			lii=rawdat[i].split()
			lii[3]='121'
			indqs.append(lii)
		elif "504203230" not in rawdat[i] and "504202950" not in rawdat[i] and "121692020" not in rawdat[i]:	
			indqs.append(rawdat[i].split())
	indqs=indqs[:-1]
	for item in indqs:
		item.append(order[item[0]])
	
	indqs.sort(key=operator.itemgetter(-1))			
	for item in indqs:
		fs.write(' '.join(item[:-1])+'\n')
	fs.close()
	##os.chdir('/home/emilyjane/Documents/distruct1.1')
	##so for some stupid reason it was getting mad when I tried to put the .ps file inthe same folder as the input files. So they all just go to the distruct folder, named whatever your results file was named
	##Chqange num pops in darwparams
	os.system('/home/ejbm/Documents/distruct1.1/distructLinux1.1 -d /home/ejbm/Documents/distruct1.1/drawparams  -K %i -N %i  -p popq.popq -i indivq.indivq -o %s.ps'%(int(K), int(n)-3, nam))
