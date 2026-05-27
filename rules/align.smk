# ===================================================================
# Part 2: BWA 原始比对
# ===================================================================

rule bwa_mem2_index:
    """
    为 BWA-MEM2 构建参考基因组索引
    """
    input:
        fasta=REF_FASTA,
        fai=f"{REF_FASTA}.fai"
    output:
        idx=f"{config['reference']['bwa_mem2_index']}.0123"
    log:
        "logs/bwa_mem2_index.log"  
    params:
        prefix=config["reference"]["bwa_mem2_index"]
    threads: 6
    cache: True 
    resources:
        slurm_partition="q_fat",
        mem_mb=200000,
        mem_gb=200,         # 也可以使用 GB 单位
        disk_mb=100000     # 磁盘空间需求（可选）

    shell:
        """        
        bwa-mem2 index \
            -p {params.prefix} \
            {input.fasta} \
            > {log} 2>&1
        
        echo "BWA-MEM2 index built: {params.prefix}" >> {log}
        """


# ===== BWA-MEM2 比对，直接输出 BAM =====
rule bwa_mem2_mem:
    """
    BWA-MEM2 比对，直接输出 BAM（未排序）
    """
    input:
        r1=f"{CLEAN_DIR}/{{sample}}_R1.umi.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_R2.umi.fastq.gz", 
        idx=f"{config['reference']['bwa_mem2_index']}.0123"
    output:
        bam=f"{BAM_DIR}/{{sample}}.raw.unsorted.bam"
    log:
        f"{LOG_DIR}/bwa_mem2_{{sample}}.log"
    threads: 6
    resources:
        slurm_partition="q_fat,q_fat_l",
        mem_mb=120000,
        mem_gb=120,         # 也可以使用 GB 单位
        disk_mb=100000     # 磁盘空间需求（可选）
    params:
        index_prefix=config["reference"]["bwa_mem2_index"],
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}_ctDNA\tPL:MGISEQ"
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            -M \
            -R "{params.rg}" \
            {params.index_prefix} \
            {input.r1} \
            {input.r2} \
            2> {log} \
        | samtools view -bS \
            -o {output.bam}
        """


# ===== 排序 =====
rule samtools_sort:
    """
    BAM 排序
    """
    input:
        bam=f"{BAM_DIR}/{{sample}}.raw.unsorted.bam"
    output:
        bam=f"{BAM_DIR}/{{sample}}.raw.bam",
        bai=f"{BAM_DIR}/{{sample}}.raw.bam.bai",
        flagstat=f"{QC_DIR}/{{sample}}_raw_flagstat.txt"
    log:
        f"{LOG_DIR}/sort_{{sample}}.log"
    threads: 6
    resources:
        slurm_partition="q_fat,q_fat_l",
        mem_mb=120000,
        mem_gb=120,         # 也可以使用 GB 单位
        disk_mb=100000     # 磁盘空间需求（可选）
    params:
        tmp=TMP_DIR
    shell:
        """
        samtools sort \
            -@ {threads} \
            -O BAM \
            -o {output.bam} \
            {input.bam} \
            2> {log}
        
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}

        echo "BWA-MEM2 alignment completed at $(date)" >> {log}
        """

# ===================================================================
# Part 4: UMI Deduplication
# ===================================================================

rule umi_dedup:
    """
    基于 UMI 标签进行去重
    - 使用 directional 方法（适合链特异性 ctDNA 文库）
    - 编辑距离阈值设为 1（允许 1bp 测序错误）
    """
    input:
        bam=f"{BAM_DIR}/{{sample}}.raw.bam"
    output:
        dedup_bam=f"{BAM_DIR}/{{sample}}.dedup.bam",
        metrics=f"{QC_DIR}/{{sample}}_dedup_metrics_per_umi.tsv"
    log:
        f"{LOG_DIR}/umi_dedup_{{sample}}.log"
    params:
        method="directional",      # directional, unique, or percentile
        edit_distance=1,            # UMI 编辑距离阈值
        metrics=f"{QC_DIR}/{{sample}}_dedup_metrics"
    threads: 6
    resources:
        slurm_partition="q_fat,q_fat_l",
        mem_mb=120000,
        mem_gb=120,         # 也可以使用 GB 单位
        disk_mb=100000     # 磁盘空间需求（可选）
    run:
        # 执行去重
        shell(f"""
        umi_tools dedup \
            --method {params.method} \
            --edit-distance-threshold {params.edit_distance} \
            --output-stats={params.metrics} \
            --log {log} \
            -I {input.bam} \
            -S {output.dedup_bam} 2>> {log}
        """)
        
        # 索引去重后的 BAM
        shell(f"samtools index {output.dedup_bam} 2>> {log}")
        
        # 计算并记录去重统计
        shell(f"""
        echo "=== UMI Deduplication Statistics for {{wildcards.sample}} ===" >> {log}
        echo "Total reads before dedup: $(samtools view -c {input.bam})" >> {log}
        echo "Reads after dedup: $(samtools view -c {output.dedup_bam})" >> {log}
        echo "Deduplication rate: $(echo "scale=2; (1 - $(samtools view -c {output.dedup_bam})/$(samtools view -c {input.bam})) * 100" | bc)%" >> {log}
        """)


# ===================================================================
# Part 5: Post-Deduplication QC
# ===================================================================

rule flagstat_dedup:
    """对去重后的 BAM 文件进行 flagstat 统计"""
    input:
        bam=f"{BAM_DIR}/{{sample}}.dedup.bam"
    output:
        flagstat=f"{QC_DIR}/{{sample}}_dedup_flagstat.txt"
    log:
        f"{LOG_DIR}/flagstat_{{sample}}.log"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        """

rule multiqc_report:
    """汇总所有 QC 报告"""
    input:
        # FastQC 报告
        expand(f"{QC_DIR}/{{sample}}_R1_fastqc.html", sample=ALL_SAMPLES),
        expand(f"{QC_DIR}/{{sample}}_R2_fastqc.html", sample=ALL_SAMPLES),
        expand(f"{QC_DIR}/{{sample}}_R1.clean_fastqc.html", sample=ALL_SAMPLES),
        expand(f"{QC_DIR}/{{sample}}_R2.clean_fastqc.html", sample=ALL_SAMPLES),
        # fastp 报告
        expand(f"{QC_DIR}/{{sample}}_fastp.json", sample=ALL_SAMPLES),
        # 去重统计
        expand(f"{QC_DIR}/{{sample}}_dedup_metrics_per_umi.tsv", sample=ALL_SAMPLES)
    output:
        multiqc_html=f"{QC_DIR}/multiqc_report.html"
    log:
        f"{LOG_DIR}/multiqc.log"
    shell:
        """
        multiqc -f \
            --outdir {QC_DIR} \
            {QC_DIR} {BAM_DIR} \
            > {log} 2>&1
        """
