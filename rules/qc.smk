# ===================================================================
# Part 1: 质量控制与 UMI 处理
# ===================================================================

rule fastqc_raw:
    """原始数据质量评估（过滤前）"""
    input:
        r1=f"{config['paths']['fastq_dir']}/{{sample}}/{{sample}}_R1.fq.gz",
        r2=f"{config['paths']['fastq_dir']}/{{sample}}/{{sample}}_R2.fq.gz"
    output:
        html1=f"{QC_DIR}/{{sample}}_R1_fastqc.html",
        html2=f"{QC_DIR}/{{sample}}_R2_fastqc.html"
    log:
        f"{LOG_DIR}/fastqc_raw_{{sample}}.log"
    threads: 2
    params:
        outdir=QC_DIR
    shell:
        """
        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            --quiet \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

# ===================================================================
# UMI 提取（先提取，再过滤）
# 说明：先提取 UMI 可以：
#   1. 避免 UMI 序列在过滤中被截断或修改
#   2. 保留完整的 UMI 信息用于下游去重
#   3. 将 UMI 添加到 read name 中，后续步骤自动保留
# ===================================================================

rule extract_umis:
    """
    提取双端 UMI 并合并为单个标签（在过滤之前）
    - 从原始 reads 中提取 UMI
    - 将 UMI 添加到 read name 中
    - 输出带 UMI 标签的 FASTQ 文件（仍包含接头和低质量碱基）
    """
    input:
        r1=f"{config['paths']['fastq_dir']}/{{sample}}/{{sample}}_R1.fq.gz",
        r2=f"{config['paths']['fastq_dir']}/{{sample}}/{{sample}}_R2.fq.gz"
    output:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.umi_raw.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.umi_raw.fastq.gz"
    log:
        f"{LOG_DIR}/umi_extract_{{sample}}.log"
    params:
        # UMI 模式配置（根据实际文库修改）
        umi_pattern = lambda wildcards: r'(?P<umi_1>.{6})'
    threads: config["resources"]["default_threads"]
    shell:
        """
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params.umi_pattern}' \
            --bc-pattern2='{params.umi_pattern}' \
            --stdin {input.r1} \
            --read2-in {input.r2} \
            --read2-out {output.r2} \
            --stdout {output.r1} \
            --ignore-read-pair-suffixes \
            --log {log} 2>&1
        
        echo "UMI extraction completed for {wildcards.sample}" >> {log}
        echo "Pattern R1: {params.umi_pattern}" >> {log}
        echo "Pattern R2: {params.umi_pattern}" >> {log}
        """

rule fastp_trim_umi:
    input:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.umi_raw.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.umi_raw.fastq.gz"
    output:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.umi.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.umi.fastq.gz",
        json=f"{QC_DIR}/{{sample}}_fastp.json",
        html=f"{QC_DIR}/{{sample}}_fastp.html"
    log:
        f"{LOG_DIR}/fastp_umi_{{sample}}.log"
    threads: config["resources"]["default_threads"]
    params:
        adapter_r1 = config.get("adapters", {}).get("r1", "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"),
        adapter_r2 = config.get("adapters", {}).get("r2", "AAGTCGGATCGTAGCCATGTCGTTCTGAGCCAAGGAGTTG")
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --adapter_sequence {params.adapter_r1} \
            --adapter_sequence_r2 {params.adapter_r2} \
            --trim_poly_x \
            --cut_front \
            --cut_front_window_size 1 \
            --cut_front_mean_quality 20 \
            --cut_tail \
            --cut_tail_window_size 4 \
            --cut_tail_mean_quality 20 \
            --overlap_len_require 15 \
            --overlap_diff_limit 3 \
            --correction \
            --length_required 50 \
            --average_qual 25 \
            --n_base_limit 5 \
            --unpaired1 {output.r1}.unpaired.fq \
            --unpaired2 {output.r2}.unpaired.fq \
            2>> {log}
        
        echo "fastp filtering completed for {wildcards.sample}" >> {log}
        python3 -c "
import json
with open('{output.json}') as f:
    d = json.load(f)
    s = d.get('summary', {{}})
    before = s.get('before_filtering', {{}})
    after = s.get('after_filtering', {{}})
    fr = d.get('filtering_result', {{}})
    print(f'Before filtering: {{before.get(\"total_reads\", \"?\")}} reads')
    print(f'After  filtering: {{after.get(\"total_reads\", \"?\")}} reads')
    print(f'Passed: {{fr.get(\"passed_filter_reads\", \"?\")}}')

" >> {log} 2>&1 || true
        """



rule fastqc_clean:
    """过滤后数据质量评估（带 UMI 的 clean reads）"""
    input:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.umi.clean.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.umi.clean.fastq.gz"
    output:
        html1=f"{QC_DIR}/{{sample}}_R1.umi.clean_fastqc.html",
        html2=f"{QC_DIR}/{{sample}}_R2.umi.clean_fastqc.html",
    log:
        f"{LOG_DIR}/fastqc_clean_{{sample}}.log"
    threads: 2
    params:
        outdir=QC_DIR
    shell:
        """
        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            --quiet \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """
