#!/usr/bin/env snakemake
# -*- coding: utf-8 -*-
"""
ctDNA Panel Analysis Pipeline with UMI
=======================================
完整的 ctDNA 分析流程：
fastp质控 → BWA比对 → fgbio UMI处理 → SNV/CNV/SV检测
"""

import os
import re
from pathlib import Path

# ===== 加载配置 =====
configfile: "config.yaml"

import snakemake
snakemake.ancient = True


# ===== 全局变量定义 =====
TUMOR_SAMPLES = config["samples"]["tumor"]
NORMAL_SAMPLE = config["samples"]["normal"][0]
ALL_SAMPLES = TUMOR_SAMPLES + [NORMAL_SAMPLE]


DATA_DIR = config["paths"]["data_dir"]
OUTPUT_DIR = config["paths"]["output_dir"]

QC_DIR = os.path.join(OUTPUT_DIR, "qc")
CLEAN_DIR =  os.path.join(OUTPUT_DIR, "fastq")
BAM_DIR = os.path.join(OUTPUT_DIR, "bam")
VAR_DIR = os.path.join(OUTPUT_DIR, "variants")
FINAL_DIR = os.path.join(OUTPUT_DIR, "final")
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")
TMP_DIR = os.path.join(OUTPUT_DIR, "tmp")
REPORT_DIR = os.path.join(OUTPUT_DIR, "reports")


# 参考文件
REF_FASTA = config["reference"]["fasta"]
TARGET_BED = config["reference"]["target_bed"]


#import ipdb;ipdb.set_trace()

# ===== 通配符约束 =====
wildcard_constraints:
    sample="|".join(ALL_SAMPLES),
    tumor="|".join(TUMOR_SAMPLES)

# ===== 定义最终目标 =====
rule all:
    """最终输出目标"""
    input:
        # BWA-MEM2 索引文件（确保索引存在）
        #f"{config['reference']['bwa_mem2_index']}.0123", 
        # 质控报告
        expand(f"{QC_DIR}/multiqc_report.html"),
            
        # 最终变异结果
        expand(f"{VAR_DIR}/snv/{{tumor}}.snv.filtered.vcf.gz", tumor=TUMOR_SAMPLES),
        expand(f"{VAR_DIR}/cnv/{{tumor}}.cnv.filtered.cns", tumor=TUMOR_SAMPLES),

        # ── panelcn.MOPS CNV（扩增子 Panel 推荐，取消注释以启用）──
        # 使用前请确认：
        #   3. config.yaml 中 tools.rscript 路径正确
        #expand(f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.gene.tsv",    tumor=TUMOR_SAMPLES),
        #expand(f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.filtered.cns", tumor=TUMOR_SAMPLES),
        #expand(f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.report.txt",   tumor=TUMOR_SAMPLES),
        
        # ── SV 检测（扩增子 Panel 推荐使用 SvABA + 靶向融合检测）──
        #    将注释中的任意一行取消注释即可启用对应工具
        #    推荐优先级：靶向融合 > SvABA > GRIDSS2 > Manta > Delly

        # 靶向融合检测（专为扩增子 Panel 设计，推荐首选）
        #expand(f"{VAR_DIR}/sv/{{tumor}}.fusions.tsv", tumor=TUMOR_SAMPLES),
        #expand(f"{VAR_DIR}/sv/{{tumor}}.fusions.summary.txt", tumor=TUMOR_SAMPLES),

        # SvABA（靶向测序专用，推荐次选）
        #expand(f"{VAR_DIR}/sv/{{tumor}}.svaba.sv.vcf.gz", tumor=TUMOR_SAMPLES),

        # GRIDSS2（assembly-based，备选）
        #expand(f"{VAR_DIR}/sv/{{tumor}}.gridss.vcf.gz", tumor=TUMOR_SAMPLES),

        # Manta（已降阈值优化，备选）
        #expand(f"{VAR_DIR}/sv/{{tumor}}.manta.vcf.gz", tumor=TUMOR_SAMPLES),

        # Delly（已修复路径冲突，备选）
        #expand(f"{VAR_DIR}/sv/{{tumor}}.delly.done", tumor=TUMOR_SAMPLES),

        # SV 汇总报告（运行前取消上面至少一个工具的注释）
        #expand(f"{VAR_DIR}/sv/{{tumor}}.sv.summary.txt", tumor=TUMOR_SAMPLES),

        # 注释和分类结果
        expand(f"{VAR_DIR}/annovar/{{tumor}}.snv.hg19_multianno.vcf", tumor=TUMOR_SAMPLES),
        
        # 最终过滤结果（ANNOVAR 原始列筛选）
        expand(f"{FINAL_DIR}/{{tumor}}.snv.filtered.txt", tumor=TUMOR_SAMPLES),
        expand(f"{FINAL_DIR}/{{tumor}}.snv.filter_summary.txt", tumor=TUMOR_SAMPLES),

        # ── 标准上报格式（3' rule + HGVS + ClinVar 风格）───────────
        expand(f"{FINAL_DIR}/{{tumor}}.snv.standard_report.txt", tumor=TUMOR_SAMPLES),

        # ── 最终 Excel 临床报告（标准格式 + Tier 分级）────────────
        expand(f"{FINAL_DIR}/{{tumor}}.snv.final_report.xlsx", tumor=TUMOR_SAMPLES),



# ===================================================================
# rules
# ===================================================================

include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/snv.smk"
include: "rules/cnv.smk"
include: "rules/sv.smk"
include: "rules/annovar.smk"


# ===================================================================
# 主程序入口
# ===================================================================
if __name__ == "__main__":
    # 创建必要的输出目录
    dirs_to_create = [
        QC_DIR, CLEAN_DIR, RAW_BAM_DIR, UMI_DIR, FINAL_BAM_DIR,
        VAR_DIR, FINAL_DIR, LOG_DIR, TMP_DIR, REPORT_DIR
    ]
    for dir_path in dirs_to_create:
        os.makedirs(dir_path, exist_ok=True)
    
    # 创建变异子目录
    for subdir in ["snv", "cnv", "sv"]:
        os.makedirs(os.path.join(VAR_DIR, subdir), exist_ok=True)
    
    print("=" * 60)
    print("ctDNA Panel Analysis Pipeline")
    print("=" * 60)
    print(f"Tumor samples: {', '.join(TUMOR_SAMPLES)}")
    print(f"Normal sample: {NORMAL_SAMPLE}")
    print(f"Output directory: {config['paths']['output_dir']}")
    print("=" * 60)