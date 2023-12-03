# Unzip (so we can zip later again with ), and gz
cd raw
parallel --dry-run gunzip {} ::: *gz
parallel gunzip {} ::: *gz

parallel --dry-run bgzip {} ::: *vcf
parallel bgzip {} ::: *vcf


# Index with tabix
chr_l=("1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "X"  "Y")
for i in "${chr_l[@]}"
do
INPUT="raw/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30.gvcf.gz"
tabix -p vcf $INPUT &
done


# Stats
chr_l=("1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "X"  "Y")
for i in "${chr_l[@]}"
do
INPUT="raw/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30.gvcf.gz"
OUTPUT="raw-stats/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30.gvcf.gz.txt"
bcftools stats $INPUT > $OUTPUT &
done


# Filter for b minimum non-reference allele count >=1
mkdir filter-mnref1
chr_l=("1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "X"  "Y")
for i in "${chr_l[@]}"
do
INPUT="raw/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30.gvcf.gz"
OUTPUT="filter-mnref1/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1.gvcf.gz"
bcftools view -c 1 -O z -o $OUTPUT $INPUT &
done

mkdir filter-mnref1-stats
for i in "${chr_l[@]}"
do
INPUT="filter-mnref1/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1.gvcf.gz"
OUTPUT="filter-mnref1-stats/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1.gvcf.gz.txt"
bcftools stats $INPUT > $OUTPUT &
done


# Filter for b snp q30, stats
mkdir -p filter-mnref1-b-snp-q30
chr_l=("1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "X"  "Y")
for i in "${chr_l[@]}"
do
INPUT="filter-mnref1/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1.gvcf.gz"
OUTPUT="filter-mnref1-b-snp-q30/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz"
bcftools view -m2 -M2 -V indels,mnps,ref,bnd,other -i 'QUAL>30' -O z -o $OUTPUT $INPUT &
done

mkdir -p filter-mnref1-b-snp-q30-stats
chr_l=("1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "X"  "Y")
for i in "${chr_l[@]}"
do
INPUT="filter-mnref1-b-snp-q30/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz"
OUTPUT="filter-mnref1-b-snp-q30-stats/deepvariant_chr"${i}"_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz.txt"
bcftools stats $INPUT > $OUTPUT &
done

# Concat 20 chromosomal VCF files into one, stats
mkdir -p filter-mnref1-b-snp-q30-concat 
cd filter-mnref1-b-snp-q30
ls -1 *gz | sort -n +0.15 > ../filter-mnref1-b-snp-q30-concat/concat.txt
# Manual remove 2 lines for the X and Y
OUTPUT="../filter-mnref1-b-snp-q30-concat/deepvariant_chr1-20_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz"
bcftools concat -f "../filter-mnref1-b-snp-q30-concat/concat.txt" -O z -o $OUTPUT
cd ../

mkdir -p filter-mnref1-b-snp-q30-concat-stats # Stats
INPUT="filter-mnref1-b-snp-q30-concat/deepvariant_chr1-20_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz"
OUTPUT="filter-mnref1-b-snp-q30-concat-stats/deepvariant_chr1-20_163_samples_mRatBN7.2_annot_error_removed_GQ30-mnref1-b-snp-q30.gvcf.gz.txt"
bcftools stats $INPUT > $OUTPUT




