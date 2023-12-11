#Rat linkage map, December 2023
#Pasi Rastas (c) 2023
#data is in rn6 coordinates, files:TN.gen, TN.sample, pedigree.txt

#header
awk '{print $1}' TN.sample|./transpose_tab >header.txt
#edit first two columns to CHR:POS and POS

#individual names
awk '(NR>2){print $1}' TN.sample >individuals.txt

#families with 3 offspring and one parent
awk '(NR==FNR){i[$1]}(NR!=FNR && ($1 in i)){print  $2,$3}' individuals.txt pedigree.txt |sort|uniq -c|awk '($1>=3 && $2!="NA"){print $2 "_" $3"\t"$2"\t0\t0\t1\t0";print $2 "_" $3"\t"$3"\t0\t0\t2\t0"}' >fam3.txt

grep -w -F -f individuals.txt fam3.txt|cut -f 1|uniq >fam3p.txt

#make the final pedigree
(cat fam3.txt; awk '{print $2 "_" $3"\t"$1"\t" $2 "\t" $3 "\t0"}' pedigree.txt)|grep -w -F -f fam3p.txt|sort -V >ped_final.txt

#transpose pedigree for Lep-MAP3
./transpose_tab ped_final.txt|awk '{print "CHR\tPOS\t"$0}' >ped_t.txt

#LM3 points to Lep-MAP3
export LM3=../workspace/Lep-MAP3/bin

#ParentCall2
for i in {1..20}
do
java -cp $LM3 ParentCall2 data=ped_t.txt halfSibs=1 posteriorFile=<(cat header.txt; grep -w -F chr$i TN.gen|awk -f convert.awk) removeNonInformative=1 |gzip >chr$i.call.gz &
done

#Filtering2
for i in {1..20}
do
zcat chr$i.call.gz|java -cp $LM3 Filtering2 data=- dataTolerance=0.01 |gzip >chr$i.filt.gz
done

#increase informativeness by binning markers, one marker per 10kb should be enough...
for i in {1..20}
do
zcat chr$i.filt.gz|cut -f 1,2|awk 'BEGIN{print "#binned markers"}(NR>7){if (prev==0 || $2-prev >= 10000) {++n;prev=$2}; print n}' >bin$i.txt
done

#Run OrderMarkers2 for each bin
for i in $(seq 20 -1 1)
do
echo "zcat chr$i.filt.gz | java -cp $LM3 OrderMarkers2 map=bin$i.txt data=- improveOrder=0 recombination1=0 recombination2=0 phasingIterations=3 outputPhasedData=4 hyperPhaser=1 2>/dev/null |gzip >binned$i.gz"
done|parallel --jobs 4

#order2data.awk
for i in  {1..20}
do
zcat binned$i.gz | awk -vchr=$i -f order2data.awk|gzip >bdata$i.gz
done

#SNP names
for i in {1..20}; do zcat chr$i.filt.gz|cut -f 1,2|awk '(NR>=7)' >snps$i.txt; done
for i in {1..20};do awk '(NR==FNR){s[NR-1]=$0}(NR!=FNR){print s[$1]}' snps$i.txt <(zcat bdata$i.gz) >bin_snps$i.txt; done
for i in {1..20}; do cat bin_snps$i.txt; done >bin_snps_all.txt

#get pedigree in the correct order
java -cp $LM3 OrderMarkers2 data=ped_t.txt outputPedigree=ped_t2.txt 1

#SeparateChromosomes2
for i in {1..20}
do
(cat ped_t2.txt;zcat bdata$i.gz)|java -cp $LM3 SeparateChromosomes2 data=- lodLimit=90 numThreads=16 >map$i.90.txt
done

#JoinSingles2All
for i in  {1..20}; do awk -vn=$i '(NR>1){if ($1==1) print n; else print 0}' map$i.90.txt; done >map1.chr
for i in  {1..20}; do awk '(NR>1){if ($1>0) print 1; else print 0}' map$i.90.txt; done >map.mask
(cat ped_t2.txt;for i in {1..20}; do zcat bdata$i.gz; done)|java -cp lm3/ JoinSingles2All data=- lodLimit=80 lodDifference=10 numThreads=16 map=map1.chr mask=map.mask >map_js.txt

#run OrderMarkers2
for i in {1..20}; do echo "(cat ped_t2.txt;for j in {1..20}; do zcat bdata\$j.gz; done) | java -cp $LM3 OrderMarkers2 map=map_js.txt chromosome=$i data=- phasingIterations=3 numMergeIterations=1 numThreads=4 >orders/order$i.txt 2>orders/order$i.err"; done >order_commands.txt
parallel --jobs 8 < order_commands.txt

#with different scale
for i in {1..20}; do echo "(cat ped_t2.txt;for j in {1..20}; do zcat bdata\$j.gz; done) | java -cp $LM3 OrderMarkers2 map=map_js.txt chromosome=$i data=- phasingIterations=3 numMergeIterations=1 numThreads=4 scale=0.5M/N >orders2/order$i.txt 2>orders2/order$i.err"; done >order_commands2.txt
parallel --jobs 8 < order_commands2.txt

#map to rn6 coordinates
for i in {1..20}; do awk -vn=$i '(NR!=FNR && /^[^#]/){print m[$1]"\t"$2"\t"$3}(NR==FNR){m[NR]=substr($1,1,index($1,":")-1) "\t"$2"\t"n}' bin_snps_all.txt orders/order$i.txt >orders/order$i.mapped; done

for i in {1..20}; do awk -vn=$i '(NR!=FNR && /^[^#]/){print m[$1]"\t"$2"\t"$3}(NR==FNR){m[NR]=substr($1,1,index($1,":")-1) "\t"$2"\t"n}' bin_snps_all.txt orders2/order$i.txt >orders2/order$i.mapped; done

#take chrs 1,2,4,8,11 and 12 from orders, remaining from orders2

(for i in 1 2 4 8 11 12
do
        cat orders/order$i.mapped
done
for i in 3 5 6 7 9 10 {13..20}
do
        cat orders2/order$i.mapped
done)>map_cm_v6.txt

#liftover to v7 genome
sed -e 's/\t/-/g' bin_snps_all.txt >bin_snps_all.bed
#liftover, using first https://genome.ucsc.edu/cgi-bin/hgLiftOver
#output are hglft_genome_2fe11_93a5d0.bed and hglft_genome_2fe11_93a5d0.err

awk '(pass==1){e[$1]}(pass==2){d[++n]=$1}(pass==3){if ($1 in e) print ""; else print d[++m]}' pass=1 hglft_genome_2fe11_93a5d0.err pass=2 hglft_genome_2fe11_93a5d0.bed pass=3 bin_snps_all.bed >snps_mapped.bed

paste bin_snps_all.bed snps_mapped.bed|sed -e 's/-/\t/g' -e 's/:[0-9]*//g'|awk '(NF==4)' >mapping2v7.txt

awk -vOFS="\t" '(NR==FNR){m[$1"\t"$2]=$3"\t"$4}(NR!=FNR && ($1"\t"$2 in m)){$1=m[$1"\t"$2];$2="";print}' mapping2v7.txt map_cm_v6.txt >map_cm_v7.txt


#maps with header
awk 'BEGIN{print "CONTIG\tPOS\tCHR\tMALE_POS\tFEMALE_POS"}1' map_cm_v6.txt >map_cm_v6_withheader.txt

awk 'BEGIN{print "CONTIG\tPOS\tCHR\tMALE_POS\tFEMALE_POS"}1' map_cm_v7.txt >map_cm_v7_withheader.txt

