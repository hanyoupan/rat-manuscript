# A revamped rat reference genome improves the discovery of genetic diversity in laboratory rats 

The manuscript is available from [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.13.536694v2). This repository hosts the code linked to the paper currently under review by Cell Genomics.

This repository is organized into the following folders:

## 1_compare_rn6_v_mRatBN7.2

We compared Rnor_6.0 and mRatBN7.2 using dotplots, genetic map, comparing mapping stats of whole genome sequencing data, and compared eQTL and pQTL mappings. The scripts for these analysis are placed in their own folders.  

## 2_mapping_of_whole_genome_sequencing_data

This folder contains bash commands for analyzing the WGS data.

## 3_liftover 

This folder contains script for liftover analysis between Rnor_6.0 and mRatBN7.2

## 4_strain_comparisons

This script contains scripts for comparing genetic variants between inbred strains of rats. These analysis are based on mapping of whole genome sequencing data to mRatBN7.2. The VCF files obtained from joint variant calling is available from [Zenodo](https://zenodo.org/records/10398344).   

