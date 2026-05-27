# ===================================================================
# Part 1: 质量控制
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
    #conda:
    #    "envs/environment.yaml"
    shell:
        """
        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            --quiet \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

rule fastp_trim:
    """
    使用 fastp 进行接头去除和质量过滤
    - 自动检测并去除接头序列
    - 去除低质量碱基
    - 过滤过短的 reads
    - 输出过滤后的 FASTQ 和 HTML 报告
    """
    input:
        r1=f"{config['paths']['fastq_dir']}/{{sample}}/{{sample}}_R1.fq.gz",
        r2=f"{config['paths']['fastq_dir']}/{{sample}}/{{sample}}_R2.fq.gz"
    output:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.clean.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.clean.fastq.gz",
        json=f"{QC_DIR}/{{sample}}_fastp.json",
        html=f"{QC_DIR}/{{sample}}_fastp.html"
    log:
        f"{LOG_DIR}/fastp_{{sample}}.log"
    threads: config["resources"]["default_threads"]
    #conda:
    #    "envs/environment.yaml"
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
            --length_required 50 \
            --average_qual 20 \
            --n_base_limit 5 \
            2> {log}
        
        # 记录过滤统计到日志
        echo "fastp filtering completed for {wildcards.sample}" >> {log}
        grep "reads passed filter" {output.json} >> {log} || true
        """

rule fastqc_clean:
    """过滤后数据质量评估"""
    input:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.clean.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.clean.fastq.gz"
    output:
        html1=f"{QC_DIR}/{{sample}}_R1.clean_fastqc.html",
        html2=f"{QC_DIR}/{{sample}}_R2.clean_fastqc.html",
    log:
        f"{LOG_DIR}/fastqc_clean_{{sample}}.log"
    threads: 2
    params:
        outdir=QC_DIR
    #conda:
    #    "envs/environment.yaml"
    shell:
        """
        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            --quiet \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

# UMI 正则表达式模式（6bp UMI，中间有固定碱基 T 和 A）
# 序列结构: UMI1(6bp) - T - Insert - A - UMI2(6bp)

rule extract_umis:
    """
    提取双端 UMI 并合并为单个 12bp 标签
    - 使用正则表达式定位 6bp UMI
    - 合并 R1 和 R2 的 UMI 为单个标签
    - 将合并后的 UMI 添加到 read name 中
    """
    input:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.clean.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.clean.fastq.gz"
    output:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.umi.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.umi.fastq.gz"
    log:
        f"{LOG_DIR}/umi_extract_{{sample}}.log"
    params:
        pattern=lambda wildcards: r'(?P<umi_1>.{6})T.*'
    threads: config["resources"]["default_threads"]
    shell:
        """
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params.pattern}' \
            --bc-pattern2='{params.pattern}' \
            --stdin {input.r1} \
            --read2-in {input.r2} \
            --read2-out {output.r2} \
            --stdout {output.r1} \
            --ignore-read-pair-suffixes \
            --log {log} 2>&1
        """