

# Read the filtered SNP file.
 

###########################################################
# 1. create list from SNP file.
###########################################################
inFile1 = open("data/high_impact_snps.csv", 'r') 

geneList = []
first = True
for line1 in inFile1:
	cols1 = line1.replace('"','').rstrip().split(",")
	geneList.append(cols1[8])
	if first == True:
		header1 = cols1
		first = False
genelist = list(set(geneList))
print(len(geneList), "genes found across 100+ samples")

#print (geneList)
inFile1.close()


###########################################################
# 2. Read annotation file
###########################################################
inFile2 = open("data/rattus_terms_rdo_fixed.csv", 'r')
previousLine = ""

# "data/corrected_GO.csv"
annotationDict = {}
first = True
for line2 in inFile2:
	#print(line2)
	cols2 = line2.rstrip().split(";")
	if cols2 == previousLine:
		print("duplicate removed!")
		pass
	previousLine = cols2
	
	if first == True:
		header2 = cols2
		first = False
	
	#print (cols2)
	goGene = cols2[1]
	# Only keep relevant annotations. in dict with lists
	if goGene in genelist:
		if goGene in annotationDict.keys():
			annotationDict[goGene].append(cols2)
		else:
			annotationDict[goGene] = []
			annotationDict[goGene].append(cols2)
	else:
		pass
	
print(len(annotationDict.keys()), "genes with an annotation")
print(annotationDict.keys())
inFile2.close()


###########################################################
# 3. fold together the SNPs and annotations
###########################################################

inFile3 = open("data/high_impact_snps.csv", 'r') 

outFile = open("data/high_impact_annotated.csv","w")

outFile.write("\t".join(header1)+ "\t" + "\t".join(header2)+ "\n")

# Create new snp file with annotations included.
writtenLines = 0
for line3 in inFile3:
	cols3 = line3.replace('"','').rstrip().split(",")
	gene = cols3[8]
	
	if gene in annotationDict.keys():
		#print(gene, "found")
		for outList in annotationDict[gene]:
			#print (outList)
			outFile.write("\t".join(cols3) + "\t" + "\t".join(outList) + "\n")
			writtenLines += 1
			
print("All done writing: ", writtenLines, " lines")			
outFile.close()

			


