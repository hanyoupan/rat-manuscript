
# Liftover ES1K to from Rnor6_0 to mRatBN7.2 (rn7) 

The chain_file is obtained from https://hgdownload.soe.ucsc.edu/goldenPath/rn6/liftOver/


```
input_bed="bed-files/es1k-rn6.bed"
chain_file="rn6ToRn7.over.chain.gz"
outout_bed="es1k-rn6-l2rn7.bed"
time liftover/liftover_mac ${input_bed} ${chain_file} ${outout_bed} ${outout_bed}".unlifted"
```
