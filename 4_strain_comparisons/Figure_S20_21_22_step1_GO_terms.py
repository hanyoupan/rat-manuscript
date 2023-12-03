
import itertools
import operator
import sys
import os
import time

start = time.time()
print("lets get started")




stepSize = 10000

chromlist = ['1','2','3','4','5','6','7','8','9','10','11', '12','13','14','15','16','17','18','19','20','X','Y']

for chromNumber in chromlist:
	location = 'data/deepvariant_chr'+chromNumber+'_163_samples_mRatBN7.2_annot_error_removed_GQ30_HIGH_IMPACT_only.gvcf' #deepvariant_chr1_168_samples_mRatBN7.2_annot_error_removed.gvcf

	print("running", location)

	vcf_file = open(location, 'r') # dissabled in case of accidental activation

	for line in vcf_file:
			if not line.startswith("##"):
					##print line
					header = line.split()
					break

	names = header[9:]
	out_names =  "\t".join(names) 

	## We need, range,  chrom, pos, id, ref, alt,     region_type,         impact,        nearby_gene, Out names
	out_header = "\t".join(header[0:5]) + "\t" + "region_type" + "\t" + "impact" + "\t" + "nearby_gene" + "\t" + out_names + "\n"
	out_header2 = "chromosome" + "\t" + "start" + "\t" + "end" + "\t" + out_names + "\n"
	
	gt_rolled_file = open(location+".rolled.gt", "w")
	gt_rolled_file.write(out_header2)
	
	gt_impact_file = open(location+".high.gt", "w")
	gt_impact_file.write(out_header)
	
	gq_impact_file = open(location+".high.gq", "w")
	gq_impact_file.write(out_header)
	
	id_format = "{chromo}_{pos}_{ref}_{alt}".format
	round_dict = {0.0 : "0", 1.0 : "1"}

	# create empty list for rolled means
	local_snp_numbers = [0]*len(names)
	current_position = 0

	for line_nr, line in enumerate(vcf_file, start=0):

		cols = line.rstrip().split()
		alt_alleles = cols[4].split(",")
		quality = cols[5]
		
		if int(float(quality)) < 30:
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
		for indiv in indiv_data:
			#print(indiv)
			indiv_info = indiv.split(":")

			# Genotype calls
			indiv_gt = indiv_info[genotype_index] # .replace(".", "0") was in between
			gt_list.append(indiv_gt)

			# Get if it is an alt call for the rolled overview.
			if str(indiv_gt) == "1/1" or str(indiv_gt) == "2/2":
				is_alt_list.append(1)
			else:
				is_alt_list.append(0)

			# Quality calls
			indiv_gq = indiv_info[quality_index].replace(".", "0")
			
			#print(indiv_gq)
			gq_list.append(str(indiv_gq))

			indiv_ad = list(map(int, indiv_info[allelic_depth_index].replace(".", "0").split(",")))
			ad_list.append(indiv_ad)
			ad_str_list.append(str(sum(indiv_ad)))

		infoLine = [chromosome, position, allele_ids, ref, alt,  region_type, impact, nearby_gene]
		
		############################################################################################
		# To roll the means, to find which positions for each strain have a lot of SNPs for the entire chromosome.
		############################################################################################
		
		position = int(position)
		while True:
			# If we are before the position, there is a big problem
			if position < current_position:
				print("error, we should not be able to be before the starting position!!!")
				1/0
			# If we meet the position, then we must combine the current local_snp numbers with the larger list and continue for a while
			elif position >= current_position and position < (current_position + stepSize):
				# Get the local calls
				zipped_lists = zip(local_snp_numbers, is_alt_list)
				local_snp_numbers = [x + y for (x, y) in zipped_lists]
				break
			# If we walked by the position, we must write the previous set, and start the next one.
			elif position > current_position:
				#print ("next frame!")
				string_ints = [str(int) for int in local_snp_numbers]
				gt_rolled_file.write(chromosome + "\t" + str(current_position) + "\t" + str(current_position+stepSize) + "\t" + "\t".join(string_ints) + "\n")
				local_snp_numbers = [0]*len(names)
				current_position += stepSize
				# Use a round number to jump to next step
				#current_position = round(position/stepSize)*stepSize
				#print (current_position)
			else:
				print("what do you mean else?")
				1/0

			
	
		############################################################################################
		# To extract high impact SNPs and write them to file.
		############################################################################################
		
		if impact == "HIGH": # or impact == "MODERATE"
			gt_impact_file.write("\t".join(infoLine) + "\t" +  "\t".join(gt_list) + "\n")
			gq_impact_file.write("\t".join(infoLine) + "\t" +  "\t".join(gq_list)+ "\n")
			stored_position = position
		
	end = time.time()
	print(end - start)
		
	gt_rolled_file.close()
	gt_impact_file.close()
	print("all done")	
