
import itertools
import operator
import sys
import os



chromlist = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20',"X","Y"]

	
for chromNumber in chromlist:
	# deepvariant_mRatNor1_chr12_36_samples.gvcf
	#location = 'data/deepvariant_mRatNor1_chr'+chromNumber+'_36_samples.gvcf' #mratbn7
	location = 'data/deepvariant_rn6_chr'+chromNumber+'_36_samples.gvcf' #rn6
	print("running", location)

	vcf_file = open(location, 'r') # dissabled in case of accidental activation

	for line in vcf_file:
			if not line.startswith("##"):
					##print line
					header = line.split()
					break
					
	#print header[9:]
	names = header[9:]
	out_header = "X" + "\t" + "\t".join(names) + "\n"
	#out_header2 = "\t" + "\t".join(names) +  "\t" + "call_qual" +  "\t" + "call_depth" +  "\t" + "call_custom" + "\tRNC_" + "\tRNC_".join(names) + "\n"

	gt_out_file = open(location+".gt", "w")
	gt_out_file.write(out_header)
	ad_out_file = open(location+".ad", "w")
	ad_out_file.write(out_header)
	gq_out_file = open(location+".gq", "w")
	gq_out_file.write(out_header)
	
	adj_gt_outfile = open(location+".adj", "w")
	adj_gt_outfile.write(out_header)
	

		
	id_format = "{chromo}_{pos}_{ref}_{alt}".format

	round_dict = {0.0 : "0", 1.0 : "1"}

	for line_nr, line in enumerate(vcf_file, start=0):

		cols = line.rstrip().split()
		alt_alleles = cols[4].split(",")
		quality = cols[5]
		if int(quality) < 30:
			continue

		allele_ids = cols[2]

		line_format = cols[8].split(":")
		genotype_index = line_format.index("GT") # index of the genotypes in the individual
		allelic_depth_index = line_format.index("AD") # index of the allelic depth estimation in the individual
		quality_index = line_format.index("GQ") # index of the allelic depth estimation in the individual
		RNC_index = line_format.index("RNC") 
		indiv_data = cols[9:]
		
		gt_list = []
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

			# Quality calls
			indiv_gq = indiv_info[quality_index].replace(".", "0")
			#print(indiv_gq)
			gq_list.append(int(indiv_gq))

			if int(indiv_gq) >= 30:
				adjusted_gt_list.append(indiv_gt)
			else:
				adjusted_gt_list.append("x/x")

			indiv_ad = list(map(int, indiv_info[allelic_depth_index].replace(".", "0").split(",")))
			ad_list.append(indiv_ad)
			ad_str_list.append(str(sum(indiv_ad)))


		gt_out_file.write(allele_ids+"\t")
		ad_out_file.write(allele_ids+"\t")
		gq_out_file.write(allele_ids+"\t")
		adj_gt_outfile.write(allele_ids+"\t")


		gt_out_file.writelines("\t".join(gt_list))
		ad_out_file.writelines("\t".join(ad_str_list))
		gq_out_file.writelines("\t".join(map(str, gq_list)))
		adj_gt_outfile.writelines("\t".join(adjusted_gt_list))

		gt_out_file.write("\n")
		ad_out_file.write("\n")
		gq_out_file.write("\n")
		adj_gt_outfile.write("\n")
	
	print("all done")	

	gt_out_file.close()
	ad_out_file.close()
	gq_out_file.close()
	adj_gt_outfile.close()
