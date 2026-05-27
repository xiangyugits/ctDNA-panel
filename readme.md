## 下载参考基因组

```sh
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
cat chr{1..22}.fa chrX.fa chrY.fa chrM.fa > hg19.fa
rm -rf chr*.fa chromFa.tar.gz

module load samtools

srun -p q_cn samtools faidx hg19.fa


wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
# 重命名
mv 00-All.vcf.gz dbsnp_hg19.vcf.gz


```



## test

```sh
conda activate smk-ctdna

srun -p q_cn seqtk sample -s42 ../202611/202611_R1.fq.gz 10000 | gzip > test_R1.fq.gz
srun -p q_cn seqtk sample -s42 ../202611/202611_R2.fq.gz 10000 | gzip > test_R2.fq.gz

# ~/.config/snakemake/slurm/config.yaml

conda activate smk-ctdna

module load bwa-mem
module load samtools
module load umi_tools
module load multiqc
module load fastqc/0.11.9
module load fastp
module load GATK4
module load manta
module load bcftools

module load tabix
srun -p q_cn bgzip loci.bed
srun -p q_cn tabix -p bed loci.bed
srun -p q_cn gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict



snakemake -j10 \
--snakefile /home/zhangli_lab/zhouxiangyu/DATA/workflow/ctDNA-panel/Snakefile \
--configfile config.yaml \
--profile slurm -np

```



## run 
```sh
conda activate smk-ctdna

module load bwa-mem
module load samtools
module load umi_tools
module load multiqc
module load fastqc/0.11.9
module load fastp
module load GATK4
module load manta
module load bcftools


snakemake -j20 \
--snakefile /home/zhangli_lab/zhouxiangyu/DATA/workflow/ctDNA-panel/Snakefile \
--configfile /home/zhangli_lab/zhouxiangyu/DATA/workflow/ctDNA-panel/config.yaml \
--profile slurm -np


```


