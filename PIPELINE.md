# ctDNA Panel Analysis Pipeline - 详细分析流程

## 📊 完整分析工作流

```
┌─────────────────────────────────────────────────────────────────┐
│                   ctDNA Panel 测序数据                           │
│              (FastQ 文件：R1.fq.gz, R2.fq.gz)                   │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                    ┌──────────▼──────────┐
                    │   质量控制（QC）     │
                    └──────────┬──────────┘
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
                    ┌──────────▼──────────┐
                    │  SAM → BAM 转换    │
                    │  排序和索引         │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │  UMI 处理 (fgbio)   │
                    │  去重/错误纠正      │
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

---

## 🔍 详细模块说明

### 1️⃣ 质量控制（QC）模块

#### 步骤 1.1：Fastp 质控
```bash
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz \
      -o clean_R1.fq.gz -O clean_R2.fq.gz \
      --cut_front --cut_tail \
      --length_required 50 \
      --json output.json \
      --html output.html
```

**功能**：
- 去除低质量碱基
- 修剪 FastQ 头尾
- 去除适配体序列
- 过滤短序列 (<50 bp)
- 生成质控统计

**输出**：
- `*_R1.fq.gz`, `*_R2.fq.gz`: 清洁 FastQ
- `*.json`, `*.html`: 质控报告

---

#### 步骤 1.2：FastQC 质量评估
```bash
fastqc -t 4 clean_R1.fq.gz clean_R2.fq.gz
```

**功能**：
- 评估序列质量分布
- 检查 GC 含量
- 识别重复序列
- 检测污染迹象

**输出**：
- `*_fastqc.zip`: 详细质控数据

---

#### 步骤 1.3：MultiQC 汇总报告
```bash
multiqc . --outdir qc --interactive
```

**功能**：
- 汇总所有样本 QC 数据
- 生成交互式 HTML 报告
- 快速识别异常样本

**输出**：
- `multiqc_report.html`: 全局 QC 报告

---

### 2️⃣ 序列比对（Mapping）模块

#### 步骤 2.1：创建参考基因组索引
```bash
# BWA-MEM2 索引
bwa-mem2 index hg19.fa

# GATK 序列字典
gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict

# Samtools 索引
samtools faidx hg19.fa
```

**输出**：
- `hg19.fa.0123`, `hg19.fa.amb`, 等（BWA 索引文件）
- `hg19.dict`: GATK 字典
- `hg19.fa.fai`: Samtools 索引

---

#### 步骤 2.2：BWA-MEM2 比对
```bash
bwa-mem2 mem -t 16 -R "@RG\tID:sample\tSM:sample" \
         hg19.fa \
         clean_R1.fq.gz clean_R2.fq.gz | \
samtools view -bS -o sample.bam
```

**参数说明**：
- `-t 16`: 使用 16 个线程
- `-R`: 添加读群组信息
- `@RG`: ID（读群组ID）、SM（样本名）

**输出**：
- `*.bam`: 未排序的 BAM 文件

---

#### 步骤 2.3：BAM 排序和索引
```bash
# 坐标排序
samtools sort -@ 4 -o sample.sorted.bam sample.bam

# 生成索引
samtools index sample.sorted.bam
```

**输出**：
- `*.sorted.bam`: 排序后的 BAM
- `*.sorted.bam.bai`: BAM 索引

---

#### 步骤 2.4：fgbio UMI 处理
```bash
# 提取 UMI
fgbio ExtractUmisFromBam \
  -i sample.sorted.bam \
  -o sample.umi.bam \
  -r 1M+P \
  -s RXrx

# 标记 PCR 重复
fgbio GroupReadsByUmi \
  -i sample.umi.bam \
  -o sample.grouped.bam \
  -s adjacency \
  --edits 1

# 去重和错误纠正
fgbio CallMolecularConsensusReads \
  -i sample.grouped.bam \
  -o sample.consensus.bam \
  -M 1 \
  -e 0.05
```

**功能**：
- 提取分子标签（UMI）
- 按 UMI 分组
- 去除 PCR 重复
- 生成分子共识序列
- 提高SNV检测准确性

**输出**：
- `*.consensus.bam`: 去重后的 BAM

---

### 3️⃣ SNV 检测（Single Nucleotide Variants）模块

#### 步骤 3.1：创建候选位点列表
```bash
# 从所有样本获取候选位点
samtools mpileup -f hg19.fa \
  tumor.bam normal.bam | \
  bcftools call -m > candidates.vcf
```

---

#### 步骤 3.2：MuTect2 SNV 检测
```bash
gatk Mutect2 \
  -R hg19.fa \
  -I tumor.bam \
  -I normal.bam \
  -normal NC \
  -L loci.bed \
  -germline-resource gnomAD.vcf.gz \
  -pon panel_of_normals.vcf.gz \
  -O tumor.snv.vcf.gz
```

**参数说明**：
- `-I`: 输入 BAM 文件
- `-normal`: 正常样本名称
- `-L`: 目标区域 BED 文件
- `-germline-resource`: 种系变异资源
- `-pon`: 正常样本 panel（用于识别伪阳性）

**输出**：
- `*.snv.vcf.gz`: 原始 SNV VCF

---

#### 步骤 3.3：SNV 过滤
```python
# 过滤标准
过滤规则:
  1. VAF (等位基因频率) ≥ 0.005 (0.5%)
  2. TLOD (肿瘤对照比值) ≥ 30
  3. 测序深度 ≥ 50
  4. 平均测序质量 ≥ 30
  5. 支持变异的 reads ≥ 10
  6. gnomAD 人群频率 ≤ 0.001

bcftools filter -i \
  "AF[*]>=0.005 & TLOD>=30 & DP>=50 & MMQ>=30 & AD[1]>=10" \
  tumor.snv.vcf.gz -o tumor.snv.filtered.vcf.gz
```

**输出**：
- `*.snv.filtered.vcf.gz`: 过滤后的 SNV

---

### 4️⃣ CNV 检测（Copy Number Variations）模块

#### 步骤 4.1：准备配体样本池
```bash
# 汇总多个正常样本创建参考
cnvkit.py batch normal.bam \
  -r hg19.fa \
  -t loci.bed \
  -o normal_reference
```

**功能**：
- 从正常样本创建 CNV 背景参考
- 减少假阳性检测

**输出**：
- `normal_reference.cnn`: 参考配置文件

---

#### 步骤 4.2：CNVKit 分析
```bash
cnvkit.py batch tumor.bam \
  -r normal_reference.cnn \
  -t loci.bed \
  -o tumor.cnv.cns \
  --diagram \
  --scatter
```

**流程**：
1. **Coverage**: 计算覆盖度
2. **Reference**: 创建参考
3. **Fix**: 标准化覆盖度
4. **Segment**: 分割连续区域
5. **Call**: 调用拷贝数

**输出**：
- `*.cnv.cns`: 分段 CNV 结果
- `*.diagram.pdf`: CNV 可视化图

---

#### 步骤 4.3：CNV 过滤
```python
# 过滤标准
过滤规则:
  1. log2 绝对值 ≥ 0.3
  2. 连续探针数 ≥ 5
  3. 置信区间 ≤ 1.0

bcftools filter -i \
  "abs(LOG2)>=0.3 & PROBES>=5 & CI<=1.0" \
  tumor.cnv.cns -o tumor.cnv.filtered.cns
```

**输出**：
- `*.cnv.filtered.cns`: 过滤后的 CNV

---

### 5️⃣ SV 检测（Structural Variations）模块

#### 步骤 5.1：Delly SV 检测
```bash
# 调用 SV
delly call -x hg19.excl \
  -g hg19.fa \
  -q 30 \
  tumor.bam normal.bam \
  -o tumor_normal.vcf

# 过滤
delly filter -f germline \
  -q 30 \
  -r 0.25 \
  tumor_normal.vcf \
  -o tumor_normal.delly.vcf
```

**检测的 SV 类型**：
- `DEL`: 删除
- `DUP`: 重复
- `INV`: 倒位
- `INS`: 插入
- `BND`: 断点

**输出**：
- `*.delly.vcf`: Delly SV 结果

---

#### 步骤 5.2：Manta SV 检测
```bash
# 配置
configManta.py --normalBam normal.bam \
               --tumorBam tumor.bam \
               --referenceFasta hg19.fa \
               --callRegions loci.bed.gz \
               --runDir manta_run

# 运行
manta_run/runWorkflow.py -m local -j 8
```

**功能**：
- 检测大型结构变异
- 整合配对末端和分割读数

**输出**：
- `*.manta.vcf.gz`: Manta SV 结果

---

#### 步骤 5.3：SV 过滤和合并
```bash
# 合并 Delly 和 Manta 结果
bcftools concat tumor.delly.vcf tumor.manta.vcf.gz | \
bcftools sort -o tumor.sv.vcf.gz

# 过滤
bcftools filter -i \
  "MAPQ>=30 & SVLEN>=50" \
  tumor.sv.vcf.gz -o tumor.sv.filtered.vcf.gz
```

**过滤标准**：
- MAPQ ≥ 30（比对质量）
- SVLEN ≥ 50（结构变异大小）

**输出**：
- `*.sv.filtered.vcf.gz`: 过滤后的 SV

---

### 6️⃣ 变异注释（Annotation）模块

#### 步骤 6.1：ANNOVAR 准备
```bash
# 格式转换
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  tumor.snv.filtered.vcf.gz > tumor.avinput
```

---

#### 步骤 6.2：ANNOVAR 注释
```bash
annotate_variation.pl \
  -buildver hg19 \
  -out tumor.snv \
  tumor.avinput \
  /annovar/humandb/

# 功能注释
table_annovar.pl \
  -buildver hg19 \
  -protocol refGene,avsnp147,clinvar,cosmic70,gnomad_genome \
  -operation g,f,f,f,f \
  -nastring . \
  -vcfinput \
  tumor.snv.filtered.vcf.gz \
  /annovar/humandb/ \
  -outfile tumor.snv
```

**注释数据库**：

| 数据库 | 功能 |
|--------|------|
| refGene | 基因注释 |
| avsnp147 | SNP 数据库 |
| clinvar | 临床意义 |
| cosmic70 | 癌症相关变异 |
| gnomad_genome | 种群频率 |
| dbscsnv11 | 剪接位点预测 |
| revel | 有害性评估 |

**输出**：
- `*.hg19_multianno.vcf`: 注释 VCF
- `*.hg19_multianno.txt`: 注释表格

---

#### 步骤 6.3：靶基因过滤
```bash
# 与靶基因列表交集
bcftools view -i "GENE in target_genes" \
  tumor.snv.hg19_multianno.vcf \
  -o tumor.snv.target.vcf
```

**目标基因列表**：
- 癌症相关基因（TP53, BRCA1/2, KRAS 等）
- 可操作性基因（用于治疗靶向）
- Panel 设计的特定基因

---

### 7️⃣ 报告生成（Reporting）模块

#### 步骤 7.1：单样本报告
```bash
python3 scripts/generate_report.py \
  --sample tumor_sample_id \
  --snv tumor.snv.filtered.vcf.gz \
  --cnv tumor.cnv.filtered.cns \
  --sv tumor.sv.filtered.vcf.gz \
  --coverage coverage_stats.txt \
  --umi-stats umi_stats.json \
  --fastp-json fastp.json \
  --output tumor_analysis_report.html
```

**报告内容**：
- 样本基本信息
- 测序质量统计
- 覆盖度分析
- UMI 去重统计
- SNV 摘要表格
- CNV 可视化
- SV 列表
- 注释结果

**输出**：
- `*_analysis_report.html`: 交互式 HTML 报告

---

#### 步骤 7.2：队列汇总报告
```bash
python3 scripts/generate_cohort_summary.py \
  --reports *_analysis_report.html \
  --output cohort_summary.html
```

**汇总内容**：
- 所有样本分析链接
- 全局统计数据
- 变异热点分析
- 样本间比较

**输出**：
- `cohort_summary.html`: 队列汇总报告

---

## 📈 关键质量指标

### QC 指标

| 指标 | 推荐值 | 说明 |
|------|--------|------|
| Q30 碱基比例 | >85% | 高质量碱基比例 |
| GC 含量 | 40-60% | 异常 GC 提示污染 |
| 重复率 | <20% | PCR 重复检测 |
| 适配体比例 | <5% | 适配体残留 |

### 比对指标

| 指标 | 推荐值 | 说明 |
|------|--------|------|
| 比对率 | >95% | 映射到参考基因组 |
| 去重后深度 | >1000x | 适用于液体活检 |
| 均匀性 | >80% | 覆盖分布均匀度 |
| MAPQ 中位数 | >30 | 比对质量 |

### UMI 指标

| 指标 | 说明 |
|------|------|
| UMI 去重率 | PCR 重复比例 |
| 分子家族大小 | UMI 分组统计 |
| 错误纠正率 | 纠正的错误比例 |

### 变异指标

| 指标 | SNV | CNV | SV |
|------|-----|-----|-----|
| 检测数量 | 10-500 | 1-50 | 1-100 |
| 平均 VAF | 0.5-5% | - | - |
| 假阳性率 | <1% | <5% | <10% |

---

## 🔄 参数优化建议

### 根据应用调整参数

#### 高灵敏度模式（液体活检）
```yaml
过滤参数:
  snv:
    min_vaf: 0.002        # 降低 VAF 阈值
    min_tlod: 20          # 降低 TLOD
    min_depth: 100        # 提高深度要求
    min_ad: 5             # 降低支持读数
```

#### 高特异性模式（临床报告）
```yaml
过滤参数:
  snv:
    min_vaf: 0.01         # 提高 VAF 阈值
    min_tlod: 50          # 提高 TLOD
    min_depth: 50         # 标准深度
    min_ad: 20            # 提高支持读数
```

---

## 📊 常见输出示例

### SNV VCF 格式示例
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr17   7577121 .       G       A       .       PASS    SOMATIC;TLOD=45.23;AF=0.025;MMQ=60;AD=10,400
chr13   32889611 .      C       T       .       PASS    SOMATIC;TLOD=38.15;AF=0.012;MMQ=60;AD=8,650
```

### CNV 结果格式示例
```
sample  chrom   start   end     log2    cn      call
tumor   chr1    100000  200000  0.58    3       gain
tumor   chr7    500000  600000  -1.02   1       loss
```

### 报告数据示例
```json
{
  "sample": "tumor_202611",
  "qc_metrics": {
    "total_reads": 50000000,
    "mapped_reads": 49500000,
    "duplicate_rate": 0.15,
    "mean_coverage": 2500
  },
  "variants": {
    "snv_count": 245,
    "cnv_count": 8,
    "sv_count": 3
  }
}
```

---

## 🎯 工作流优化

### 运行时间预估

| 步骤 | 样本数 | 预估时间 |
|------|--------|----------|
| QC | 1 | 30 min |
| 比对 | 1 | 2-3 hours |
| SNV 检测 | 1 | 1-2 hours |
| CNV 检测 | 1 | 30 min |
| SV 检测 | 1 | 1 hour |
| 注释 | 1 | 30 min |
| **总计** | **1** | **~6-8 hours** |

### 并行化策略
- 多样本并行运行：加速 `×(样本数)`
- 多线程：BWA (16 threads), MuTect2 (8 threads)
- 集群提交：SLURM/SGE 分布式计算

---

## 🔗 相关文档

- [Snakemake 文档](https://snakemake.readthedocs.io/)
- [GATK 最佳实践](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711)
- [CNVKit 文档](https://cnvkit.readthedocs.io/)
- [Delly 文档](https://github.com/dellytools/delly)
- [ANNOVAR 文档](http://annovar.openbioinformatics.org/en/latest/)

---

**版本**: 1.0  
**更新日期**: 2026-05-28
