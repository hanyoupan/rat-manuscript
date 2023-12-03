

location = 'data/rattus_terms_rdo.txt' #deepvariant_chr1_168_samples_mRatBN7.2_annot_error_removed.gvcf
names(LabelsOfInterest)
print("running", location)

in_file = open(location, 'r') # dissabled in case of accidental activation
out_file = open("data/rattus_terms_rdo_parsed.txt", 'w') 

for line in in_file:
		if not line.startswith("#"):
				print(line)
				header = line.split()
				break

print(header)
print(len(header))

out_file.write("\t".join(header)+"\n")

for line_nr, line in enumerate(in_file, start=0):

	cols = line.rstrip().split("\t")

	outline = cols[0:16]

	if len(outline)!= 16:
		#print (outline)
		#print (len(outline))
		new_out = outline + ["."]*(16-len(outline))
		out_text = "\t".join(new_out)+"\n"
		if len(new_out) != 16:
			print(line_nr)
			print(new_out)
		
	else:
		out_text = "\t".join(outline)+"\n"
		
	if out_text.count("\t") != 15:
		print(line_nr)
		print (out_text)
	out_file.write(out_text)

out_file.close()
in_file.close()

















