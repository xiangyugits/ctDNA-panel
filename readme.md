# ctDNA Panel Analysis Pipeline

一个完整的循环肿瘤 DNA（ctDNA）分析流程，使用 Snakemake 框架整合多个生物信息学工具，实现从原始测序数据到变异注释的端到端分析。

## 📋 项目概述

该流程针对 ctDNA 面板数据进行多层次变异检测和分析，包括：

- **质量控制（QC）**：Fastp 质控，MultiQC 报告
- **序列比对**：BWA-MEM2 比对，UMI 处理
- **SNV 检测**：MuTect2 体细胞单核苷酸变异检测
- **CNV 检测**：CNVKit 拷贝数变异检测
- **SV 检测**：Delly 和 Manta 结构变异检测
- **变异注释**：ANNOVAR 功能注释

## 🔄 分析流程

```
┌─────────────────────────────────────────────────────────────────┐
│                   ctDNA Panel 测序数据                           │
└──────────────────────────────┬──────────────────────────────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        │                      │                      │
    ┌───▼────┐          ┌──────▼──────┐         ┌────▼────┐
    │ Fastp  │          │  FastQC     │         │ MultiQC │
    │ 质控   │          │ 质量评估    │         │ 报告    │
    └───┬────┘          └──────┬──────┘         └────┬────┘
        │                      │                     │
        └──────────────────────┼─────────────────────┘
                               │
                    ┌──────────▼──────────┐
                    │  清洁 FastQ 数据    │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │  序列比对 (Mapping) │
                    │   BWA-MEM2          │
                    └──────────┬──────────┘
                               │
                               │
                    ┌──────────▼──────────┐
                    │  UMI 处理 umi-tools │
                    │  去重/错误纠正       │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        │                      │                      │
    ┌───▼────────┐      ┌──────▼──────┐         ┌────▼────────┐
    │ SNV 检测   │      │ CNV 检测    │         │ SV 检测     │
    │ MuTect2    │      │ CNVKit      │         │ Delly/Manta │
    └───┬────────┘      └──────┬──────┘         └────┬────────┘
        │                      │                     │
    ┌───▼────────┐      ┌──────▼──────┐         ┌────▼────────┐
    │ VCF 过滤   │      │ 结果过滤    │         │ VCF 过滤    │
    │ VAF/TLOD   │      │ log2 阈值   │         │ MAPQ/支持数 │
    └───┬────────┘      └──────┬──────┘         └────┬────────┘
        │                      │                     │
        └──────────────────────┼─────────────────────┘
                               │
                    ┌──────────▼──────────┐
                    │  ANNOVAR 注释       │
                    │ (RefGene/ClinVar)   │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        │                      │                      │
    ┌───▼────────┐      ┌──────▼──────┐         ┌────▼────────┐
    │ 注释 SNV   │      │ 注释 CNV    │         │ 注释 SV     │
    └───┬────────┘      └──────┬──────┘         └────┬────────┘
        │                      │                     │
        └──────────────────────┼─────────────────────┘
                               │
                    ┌──────────▼──────────┐
                    │  生成分析报告       │
                    │  (HTML/汇总表)      │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │  最终输出            │
                    │  VCF/CNS/HTML       │
                    └─────────────────────┘
```
## 🔧 环境配置


### 系统要求

- Linux/Unix 系统
- Python 3.7+
- Conda 或 Mamba（推荐）
- Snakemake >= 6.0

### 必需工具

| 工具 | 用途 |
|------|------|
| fastp | FastQ 质控 |
| bwa-mem2 | 序列比对 |
| samtools | BAM 处理 |
| umi-tools | UMI 处理 |
| GATK4 | MuTect2 SNV 检测 |
| bcftools | VCF 处理 |
| cnvkit | CNV 检测 |
| delly | SV 检测 |
| manta | SV 检测 |
| annovar | 变异注释 |
| multiqc | QC 报告汇总 |
| fastqc | FastQ 质量评估 |



### 安装依赖

```bash
# 创建 conda 环境
conda create -n smk-ctdna -c bioconda -c conda-forge snakemake

# 激活环境
conda activate smk-ctdna

# 安装其他依赖（详见 envs/environment.yaml）
conda env update -f envs/environment.yaml -n smk-ctdna
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


module load tabix
srun -p q_cn bgzip loci.bed
srun -p q_cn tabix -p bed loci.bed
srun -p q_cn gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict

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











