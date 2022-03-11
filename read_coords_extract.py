import os
import glob

#grab files
files = glob.glob("filtered/*_Yeast_sorted_filtered.bam")

#variables
Coords=["chr6:54167-55910" ,"chr13:830805-832885","chr2:277252-281281","chr12:288718-290775",
"chr4:461985-463960","chr16:81790-83777","chr2:274971-276775","chr13:170559-172140","chr15:549702-551697",
"chr8:292077-294889","chr2:803932-805889","chr2:797990-801076","chr9:35814-37913",
"chr13:98851-100933","chr10:72875-74320","chr1:140772-142742","chr13:886030-887537","HO:1492-4089"]

Gene_Names=["ACT1","CAT8","GAL1","GAL2","GAL3","GAL4","GAL7","GAL80",
"GCY1","HXT1","MAL31","MAL33","SUC2","TUB1","RPS14B","EFB1","GAS1","HO_GAL1"]

for i in range(0,len(files)):
	for n in range(0,len(Coords)):
			command = "".join(["samtools view ",files[i]," ",Coords[n],' | awk',''' '{print $4 "\t" length($10) "\t" $5 "\t" $9}' > ''',files[i].split('/')[0], '/', 'txt/',
	Gene_Names[n],"starts_S",files[i].split('/')[-1].split('S')[1].split('_Yeast')[0],".txt"])
			print(command)
			os.system(command)
			#        	command = "".join(["samtools view ",files[i]," ",Coords[n],' | awk',''' '{print $4 "\t" length($10) "\t" $5 "\t" $9}' > ''',files[i].split('/')[0], '/', 'txt/',
			#	Gene_Names[n],"starts_",files[i].split('/')[-1].split('_')[2],".txt"])
			# updated on 22/1/25 to split on 'S' and then '_Yeast' (and not just '_' as before) to work with files that start with SXX and not with LIBXXXX

# Now extract All chromosomes
Coords=["chr1:1-230218","chr2:1-813184","chr3:1-316620","chr4:1-1531933","chr5:1-576874","chr6:1-270161","chr7:1-1090940",
"chr8:1-562643","chr9:1-439888","chr10:1-745751","chr11:1-666816","chr12:1-1078177","chr13:1-924431","chr14:1-784333",
"chr15:1-1091291","chr16:1-948066","chrm:1-85779","HO:1-5816"]


Gene_Names=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10",
"Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chrm","ChrHO"]

for i in range(0,len(files)):
	for n in range(0,len(Coords)):
		command = "".join(["samtools view ",files[i]," ",Coords[n],' | awk',''' '{print $4 "\t" length($10) "\t" $5 "\t" $9}' > ''',files[i].split('/')[0], '/', 'txt/AllChr/',
		Gene_Names[n],"starts_S",files[i].split('/')[-1].split('S')[1].split('_Yeast')[0],".txt"])
		#print(command)
		os.system(command)
		# command = "".join(["samtools view ",files[i]," ",Coords[n],' | awk',''' '{print $4 "\t" length($10) "\t" $5 "\t" $9}' > ''',files[i].split('/')[0], '/', 'txt/AllChr/',
		# Gene_Names[n],"starts_",files[i].split('/')[-1].split('_')[2],".txt"])
		# updated on 22/1/25 to split on 'S' and then '_Yeast' (and not just '_' as before) to work with files that start with SXX and not with LIBXXXX
