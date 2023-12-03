
import itertools
import operator
import sys
import os

#cur_dir = "D:\\TriGene\\01_UniversityOfTennessee\\13_Deepvar_vs_gatk\\"
#os.chdir(cur_dir)
#location = [
#	'../09_GATK/GATK/GATK_filtered/joint_called_gvcfs_chr10.vcf.recode.vcf',
#	'../09_GATK/Deep_Variant/deepVariant_filtered/deepvariant_chr10_48_samples.gvcf.gz.recode.vcf'
#	
#][1]
#chromlist = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']

chromlist = ['X','Y']


def make_strain_call(gt_list, ad_list, gq_list, quality_cutoff=10, get_left = operator.itemgetter(0, 1, 2), get_right = operator.itemgetter(3, 4, 5)):

	summed_depths_per_sample = list(map(sum, ad_list))
	
	gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list), key=operator.itemgetter(0)))
	
	#gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list, summed_depths_per_sample), key= operator.itemgetter(3)))
	gq_ordered_gt_list_left,  gq_ordered_ad_list_left,  gq_ordered_gq_list_left,  summed_depths_per_sample_left  = zip(*sorted(zip(get_left( gt_list), get_left( ad_list), get_left( gq_list), get_left( summed_depths_per_sample)), key= operator.itemgetter(2)))
	gq_ordered_gt_list_right, gq_ordered_ad_list_right, gq_ordered_gq_list_right, summed_depths_per_sample_right = zip(*sorted(zip(get_right(gt_list), get_right(ad_list), get_right(gq_list), get_right(summed_depths_per_sample)), key= operator.itemgetter(2)))

	#gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list, summed_depths_per_sample), key= operator.itemgetter(3)))
	ad_ordered_gt_list_left,  ad_ordered_ad_list_left,  ad_ordered_gq_list_left,  summed_depths_per_sample_left  = zip(*sorted(zip(get_left( gt_list), get_left( ad_list), get_left( gq_list), get_left( summed_depths_per_sample)), key= operator.itemgetter(3)))
	ad_ordered_gt_list_right, ad_ordered_ad_list_right, ad_ordered_gq_list_right, summed_depths_per_sample_right = zip(*sorted(zip(get_right(gt_list), get_right(ad_list), get_right(gq_list), get_right(summed_depths_per_sample)), key= operator.itemgetter(3)))
	
	summed_depth_left   = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(*get_left( ad_list), fillvalue=0) ]
	summed_depth_right  = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(*get_right(ad_list), fillvalue=0) ]
	merged_summed_depth = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(summed_depth_left, summed_depth_right, fillvalue=0) ]
	
	call_left  = gq_ordered_gt_list_left[ -1]
	call_right = gq_ordered_gt_list_right[-1]
	merged_call_gq = call_left + "_" + call_right
	
	call_left  = ad_ordered_gt_list_left[ -1]
	call_right = ad_ordered_gt_list_right[-1]
	merged_call_ad = call_left + "_" + call_right
	
	merged_call = [ merged_call_gq, merged_call_ad ]
	
	call_left  = gq_ordered_gt_list_left[ -1]
	call_right = gq_ordered_gt_list_right[-1]
	merged_call_custom = "?/?"
	
	summed_quality_left = sum(get_left(gq_list))
	summed_quality_right = sum(get_right(gq_list))
	
	force_call_left  = ""
	force_call_right = ""
	if summed_quality_left == 0:
		force_call_left  = "./."
	if summed_quality_right == 0:
		force_call_right  = "./."
	
	if merged_call_ad == merged_call_gq:
		merged_call_custom = merged_call_ad
		merged_call.append(merged_call_custom)
	else:
		scale_factors = []
		scaled = []
		for idx, (ad, gq) in enumerate(zip(ad_list, gq_list), start=0):
			if idx < 3:
				summed_quality = summed_quality_left
			else:
				summed_quality = summed_quality_right
			
			scale_factor = (gq + 1.0) / (summed_quality + 1.0)
			tmp = [ depth * scale_factor for depth in ad ]
			scaled.append(tmp)
			scale_factors.append(scale_factor)
			
		#gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list, summed_depths_per_sample), key= operator.itemgetter(3)))
		custom_ordered_gt_list_left,  custom_value_left  = zip(*sorted(zip(get_left( gt_list), get_left(scaled)) , key= operator.itemgetter(1)))
		custom_ordered_gt_list_right, custom_value_right = zip(*sorted(zip(get_right(gt_list), get_right(scaled)), key= operator.itemgetter(1)))
		
		
		if force_call_left:
			call_left = force_call_left
		else:
			call_left  = custom_ordered_gt_list_left[ -1]
			
		if force_call_right:
			call_right = force_call_right
		else:
			call_right = custom_ordered_gt_list_right[-1]
		merged_call_custom = call_left + "_" + call_right
		merged_call.append(merged_call_custom)
	
	
	return merged_call
	


for chromNumber in chromlist:
	location = 'data/deepvariant_chr'+chromNumber+'_163_samples_mRatBN7.2_annot_error_removed_GQ30_HXBonly.gvcf' #deepvariant_wmiwli_chrY_6_samples.gvcf

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

	gt_out_file = open(location+".gt", "w")
	gt_out_file.write(out_header)
	ad_out_file = open(location+".ad", "w")
	ad_out_file.write(out_header)
	gq_out_file = open(location+".gq", "w")
	gq_out_file.write(out_header)
	
	adj_gt_outfile = open(location+".adj", "w")
	adj_gt_outfile.write(out_header)
	
	complex_out_file = open(location+".complex_out.txt", "w")
	complex_out_file.write("This is the file with complex calls that were removed",)
	
	#rnc_out_file = open(location+".rnc", "w")
	#rnc_out_file.write(out_header)
	
	id_format = "{chromo}_{pos}_{ref}_{alt}".format

	round_dict = {0.0 : "0", 1.0 : "1"}

	for line_nr, line in enumerate(vcf_file, start=0):

		cols = line.rstrip().split()
		alt_alleles = cols[4].split(",")
		quality = cols[5]
		if int(float(quality)) < 30:
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
