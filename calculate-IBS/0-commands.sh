## filter-b-snp-q30-concat
# IBS
INPUT="../../../vcf/163-rats/filter-mnref1-b-snp-q30-concat/deepvariant_chr1-20_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz"
OUTPUT_FOLDER="chr1-20-mnref1-b-snp-q30"
OUTPUT_P=${OUTPUT_FOLDER}/${OUTPUT_FOLDER}
mkdir $OUTPUT_FOLDER
plink --vcf $INPUT --double-id --vcf-idspace-to _ --vcf-half-call m --keep-allele-order --distance square ibs --out $OUTPUT_P




