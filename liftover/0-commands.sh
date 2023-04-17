# Liftover, es1k, 6 to 7
input_bed="/Users/hanyou/hanyou-c/projects/2022-09-02-li-lab/2022-10-10-rat-manuscript/liftover/bed-files/es1k-rn6.bed"
chain_file="/Users/hanyou/hanyou-c/projects/2022-09-02-li-lab/2022-10-10-rat-manuscript/liftover/chain-files/rn6ToRn7.over.chain.gz"
outout_bed="/Users/hanyou/hanyou-c/projects/2022-09-02-li-lab/2022-10-10-rat-manuscript/liftover/bed-files/es1k-rn6-l2rn7.bed"
time liftover/liftover_mac ${input_bed} ${chain_file} ${outout_bed} ${outout_bed}".unlifted"