
import itertools
import operator
import sys
import os
import time

start = time.time()
print("lets get started")

quality_cutoff = 30

stepSize = 10000


chromlist = ['1','2','3','4','5','6','7','8','9','10','11', '12','13','14','15','16','17','18','19','20','X','Y']

for chromNumber in chromlist:
	#location = 'hs_founders_genotype/deepvariant_chr'+chromNumber+'_168_samples_mRatBN7.2_annot_error_removed_HS_founders.gvcf' #deepvariant_chr1_168_samples_mRatBN7.2_annot_error_removed.gvcf
	location = 'data/deepvariant_chr'+chromNumber+'_168_samples_mRatBN7.2_annot_error_removed.gvcf' #deepvariant_chr1_168_samples_mRatBN7.2_annot_error_removed.gvcf
	outFile = "results/"+"foundInAll_chr" + chromNumber + ".gvcf.csv"
	print("running", location)

	vcf_file = open(location, 'r') # dissabled in case of accidental activation

	for line in vcf_file:
			if not line.startswith("##"):
					##print line
					header = line.split()
					break
					
	#print (header[0:9]) 

	names = header[9:]
	out_names =  "\t".join(names) 

	founderlist = [	"ACI/N_JL.ILM","ACI/N_AP.ILM",
	"BN/N_JL.ILM","BN/N_AP.ILM",
	"BUF/N_JL.ILM","BUF/N_AP.ILM",
	"F344/N_JL.ILM","F344/N_AP.ILM",
	"M520/N_JL.ILM","M520/N_AP.ILM",
	"MR/N_JL.ILM","MR/N_AP.ILM",
	"WN/N_JL.ILM","WKY/N_AP.ILM",
	"WKY/N_JL.ILM","WN/N_AP.ILM"
		]
	callGroups = "\t".join( founderlist + ["all_low", "all_ref", "all_het", "all_alt", "all_complex"])

	## We need, range,  chrom, pos, id, ref, alt,     region_type,         impact,        nearby_gene, Out names
	out_header = "\t".join(header[0:5]) + "\t" + "region_type" + "\t" + "impact" + "\t" + "nearby_gene" + "\t" + callGroups + "\n"


	founder_full_file = open(outFile+".founder_analysis", "w")
	founder_full_file.write(out_header)
	
	id_format = "{chromo}_{pos}_{ref}_{alt}".format
	round_dict = {0.0 : "0", 1.0 : "1"}

	# create empty list for rolled means
	local_snp_numbers = [0]*len(names)
	current_position = 0

	for line_nr, line in enumerate(vcf_file, start=0):

		cols = line.rstrip().split()
		alt_alleles = cols[4].split(",")
		quality = cols[5]
		
		if int(float(quality)) < quality_cutoff:
			continue

		chromosome = cols[0]
		position = cols[1]
		allele_ids = cols[2]
		ref = cols[3]
		alt = cols[4]
		
		line_format = cols[8].split(":")
		genotype_index = line_format.index("GT") # index of the genotypes in the individual
		allelic_depth_index = line_format.index("AD") # index of the allelic depth estimation in the individual
		quality_index = line_format.index("GQ") # index of the allelic depth estimation in the individual
		RNC_index = line_format.index("RNC") 
		indiv_data = cols[9:]
		

		information = cols[7].split("|")
		region_type = information[1]
		impact = information[2]
		nearby_gene = information[3]
		
		# range for rolled only
		#range,  chrom, pos, id, ref, alt,     region_type,         impact,        nearby_gene, Out names
		infoLine = [chromosome, position, allele_ids, ref, alt,  region_type, impact, nearby_gene]
			
		gt_list = []
		is_alt_list = []
		ad_list = []
		gq_list = []
		rnc_list = []
		adjusted_gt_list = []	
		ad_str_list = []
		
		a_low = 0
		a_ref = 0
		a_het = 0
		a_alt = 0
		a_complex = 0
		
		
		
		outList = []
		inFounders = 0
		for indiv_index in range(0,len(indiv_data),1):
			indiv = indiv_data[indiv_index]
			#print(indiv)
			indiv_info = indiv.split(":")

			# Genotype calls
			indiv_gt = indiv_info[genotype_index]
			indiv_gq = int(indiv_info[quality_index].replace(".", "0") )
			#gt_list.append(indiv_gt)

			##split the call
			callA,callB =  str(indiv_gt).split('/')

			if names[indiv_index] in founderlist:
				if float(str(indiv_gq)) < quality_cutoff:	
					outList.append( "./." )
				else:
					outList.append(indiv_gt) 
					if indiv_gt != "0/0":
						inFounders += 1
				
			if float(str(indiv_gq)) < quality_cutoff:
				a_low +=1
			elif str(indiv_gt) == "0/0":
				a_ref +=1
			elif callA == callB:
				a_alt +=1
			elif callA == "0" and callB != "0":
				a_het +=1
			else:	
				a_complex +=1
		
		infoLine = [chromosome, position, allele_ids, ref, alt,  region_type, impact, nearby_gene]
		callLine = outList + [str(a_low), str(a_ref), str(a_het), str(a_alt), str(a_complex)]
		
		if inFounders > 0:
			founder_full_file.write("\t".join(infoLine) + "\t" + "\t".join(callLine) + "\n")

		stored_position = position
		
	end = time.time()
	print(end - start)

	founder_full_file.close()	
	
	print("all done")	
