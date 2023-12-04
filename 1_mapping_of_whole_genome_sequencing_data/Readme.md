## Source of WGS data
SRA IDs are listed in Table S9 of the manuscript. The following script downloads data from SRA (using two SRA IDs as an example)  

```

for i in ERR224465 SRR7755593  ; do
    echo "downloading $i"
    ~/sra_tools/sratoolkit.2.10.7-centos_linux64/bin/fasterq-dump $i  -e 10
    sleep 5;
done

```

## Mapping linked-read data using longranger

```
#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH --partition campus-long
#SBATCH --qos campus-long 
#SBATCH --time 120:00:00
#SBATCH -N 1
#SBATCH -n 48 
#SBATCH -p long 
#SBATCH -q long 
#SBATCH -o /lustre/isaac/proj/UTHSC0013/longranger_commands_mRatNor1/lraci.log
#SBATCH -e /lustre/isaac/proj/UTHSC0013/longranger_commands_mRatNor1/lraci.err
#SBATCH -J lracin 

DIR="/lustre/isaac/proj/UTHSC0013/rat_fastas/95-39697"
ID=`cat $DIR/SampleName`
export PATH=/lustre/isaac/proj/UTHSC0013/longranger-2.2.2/:$PATH
cd /lustre/isaac/proj/UTHSC0013/longranger_results_mRatNor1/
longranger  wgs --id  $ID --vcmode gatk:/lustre/isaac/proj/UTHSC0013/gatk-4.0.6.0/gatk-package-4.0.6.0-local.jar --reference /lustre/isaac/proj/UTHSC0013/mRatNor1/refdata-mRatNor1_1.curated_primary.20200520/ --fastqs $DIR  --localcores=40 --localmem=180

```

## Mapping short Illumina data using bwa

```
#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH --partition campus  
#SBATCH --qos campus 
#SBATCH --time 24:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -o /lustre/isaac/scratch/hchen3/log/ACI-N_umich.log
#SBATCH -e /lustre/isaac/scratch/hchen3/log/ACI-N_umich.err
#SBATCH -J ACI-N_umich 

module load bwa
Datadir=/lustre/isaac/proj/UTHSC0013/rat_fastas/data_from_sra/fastq/
R1=SRR7755593_1.fastq.gz
R2=SRR7755593_2.fastq.gz
Sample=ACI-N_umich
Out=${Sample}_markdup_bn72_withM.bam

Ref="/lustre/isaac/proj/UTHSC0013/mRatBN7.1/rn7.fa"
THREADS=48
cd $Datadir 
TMPDIR=$SCRATCHDIR/${Sample}_temp

bwa mem -t ${THREADS} $Ref $R1  $R2  > ${Sample}.sam
sambamba view -S ${Sample}.sam -f bam -t ${THREADS} -o ${Sample}.unsorted.bam
sambamba sort ${Sample}.unsorted.bam -t ${THREADS} --tmpdir ${TMPDIR} -o ${Sample}.bam
java -jar /lustre/isaac/proj/UTHSC0013/bin/picard_v2.25.2.jar MarkDuplicates I=${Sample}.bam  O=$Out M=${Out}_metrics.txt
sambamba index -t ${THREADS} ${Out}

```

## Variant discovery using Deepvariant

We do this for one chromosome at a time. 
```
#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH --partition campus  
#SBATCH --qos campus 
#SBATCH --time 24:00:00
#SBATCH -N 1
#SBATCH -n 24 
#SBATCH -o $SCRATCHDIR/log/singularity_chr1.log
#SBATCH -e $SCRATCHDIR/log/singularity_chr1.err
#SBATCH -J cmchr1 

export SINGULARITY_CACHEDIR=$SCRATCHDIR 
export REF="/lustre/isaac/proj/UTHSC0013/mRatNor1/refdata-mRatNor1_1.curated_primary.20200520/fasta/genome.fa"

#Pull the image.
#singularity pull docker://google/deepvariant:"${BIN_VERSION}"

cd /lustre/isaac/proj/UTHSC0013/deepvar_singularity

for SAMPLENAME in  `cat ./samples`; do 
    if [ -d $SCRATCHDIR/tmp_chr1_${SAMPLENAME} ] ; then 
        rm  $SCRATCHDIR/tmp_chr1_${SAMPLENAME}/*
    else 
        mkdir $SCRATCHDIR/tmp_chr1_${SAMPLENAME}
    fi
    if [ ! -f "./deepvar140_results/deepvariant_${SAMPLENAME}_deepvar140_chr1.g.vcf.gz" ] ; then 
        singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
              deepvariant_1.4.0.sif \
              /opt/deepvariant/bin/run_deepvariant \
             --model_type=WGS \
             --ref="${REF}" \
             --reads="../longranger_results_mRatNor1/${SAMPLENAME}/outs/phased_possorted_bam.bam" \
             --regions="chr1" \
             --intermediate_results_dir="$SCRATCHDIR/tmp_chr1_${SAMPLENAME}" \
             --output_vcf="./deepvar140_results/deepvariant_${SAMPLENAME}_deepvar140_chr1.vcf.gz" \
            --output_gvcf="./deepvar140_results/deepvariant_${SAMPLENAME}_deepvar140_chr1.g.vcf.gz" \
             --num_shards=20
    fi
done
```

##  Joint variant calling using DLNexus

```
#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH --partition campus
#SBATCH --qos campus
#SBATCH --time 24:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -o /lustre/isaac/scratch/hchen3/log/glnexuslog
#SBATCH -e /lustre/isaac/scratch/hchen3/log/glnexuserr
#SBATCH -J glnexus 

cd /lustre/isaac/proj/UTHSC0013/glnexus/

for i in `seq 1 20; echo X; echo Y`; do
        cd chr$i
        rm -rf ./GLnexus.DB
        Samples=`ls deepvariant_*g.vcf.gz |wc -l`
        echo "$Samples g.vcfs found"
        if [ $Samples == 168 ] ; then
                if [ ! -e deepvariant140_168_rats_${i}.bcf ] ; then
                        echo "run glnexus"
                        glnexus_cli141 --config DeepVariant deep*.g.vcf.gz  > deepvariant140_168_rats_chr${i}.bcf
                        bcftools view deepvariant140_168_rats_chr${i}.bcf | gzip - >  deepvariant140_168_rats_chr${i}.gvcf.gz
                else
                        echo "bcf file already exist"
                fi
        else
                echo "sample count is not $Samples, skip "
        fi
        cd ..
done

```
