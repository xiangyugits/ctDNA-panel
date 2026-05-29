# ===================================================================
#  SNV/Indel 检测（Mutect2 配对模式）
# ===================================================================

rule call_snv_mutect2:
    """
    使用 Mutect2 检测体细胞 SNV 和 Indel
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        ref=REF_FASTA,
        bed=TARGET_BED,
        germline=config["reference"]["germline_resource"],
        germline_idx=config["reference"]["germline_resource_idx"]
    output:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.mutect2.vcf.gz",
        tbi=f"{VAR_DIR}/snv/{{tumor}}.mutect2.vcf.gz.tbi",
        stats=f"{VAR_DIR}/snv/{{tumor}}.mutect2.stats.tar.gz"
    log:
        f"{LOG_DIR}/mutect2_{{tumor}}.log"
    threads: config["resources"]["mutect_threads"]
    resources:
        slurm_partition="q_fat,q_fat_l",
        mem_mb=120000
    params:
        normal_name=NORMAL_SAMPLE,
        tmp_dir=TMP_DIR
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_mb}m -Djava.io.tmpdir={params.tmp_dir}" \
            Mutect2 \
            --reference {input.ref} \
            --input {input.tumor_bam} \
            --input {input.normal_bam} \
            --normal-sample {params.normal_name} \
            --germline-resource {input.germline} \
            --intervals {input.bed} \
            --f1r2-tar-gz {output.stats} \
            --output {output.vcf} \
            > {log} 2>&1
        """

rule filter_snv_mutect2:
    """过滤 Mutect2 结果"""
    input:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.mutect2.vcf.gz",
        stats=f"{VAR_DIR}/snv/{{tumor}}.mutect2.stats.tar.gz",
        ref=REF_FASTA
    output:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.mutect2.filtered.vcf.gz",
        tbi=f"{VAR_DIR}/snv/{{tumor}}.mutect2.filtered.vcf.gz.tbi"
    log:
        f"{LOG_DIR}/filter_snv_{{tumor}}.log"
    threads: 4
    params:
        contamination=config["variant_calling"]["snv"]["contamination_estimate"]
    resources:
        mem_mb=8000
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_mb}m" \
            FilterMutectCalls \
            --reference {input.ref} \
            --variant {input.vcf} \
            --contamination-estimate {params.contamination} \
            --output {output.vcf} \
            > {log} 2>&1
        """

rule ctdna_comprehensive_filter:
    """
    对 Mutect2 结果进行 ctDNA 综合过滤
    """
    input:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.mutect2.filtered.vcf.gz",
        tbi=f"{VAR_DIR}/snv/{{tumor}}.mutect2.filtered.vcf.gz.tbi"
    output:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.snv.filtered.vcf.gz",
        tbi=f"{VAR_DIR}/snv/{{tumor}}.snv.filtered.vcf.gz.tbi"
    log:
        f"{LOG_DIR}/ctdna_comprehensive_filter_{{tumor}}.log"
    params:
        min_af=config["filtering"]["snv"]['min_vaf'],
        min_tlod=config["filtering"]["snv"]['min_tlod'],
        min_depth=config["filtering"]["snv"]['min_depth'],
        min_mq=config["filtering"]["snv"]['min_mq'],
        min_ad=config["filtering"]["snv"]['min_ad']
    threads: 1
    resources:
        mem_mb=2000
    shell:
        """
        echo "=== ctDNA Filter: {wildcards.tumor} ===" > {log}
        echo "Total: $(bcftools view -H {input.vcf} 2>> {log} | wc -l)" >> {log}
        echo "PASS:  $(bcftools view -f PASS -H {input.vcf} 2>> {log} | wc -l)" >> {log}
        echo "Params: AF>={params.min_af} DP>={params.min_depth} TLOD>={params.min_tlod} MMQ>={params.min_mq} AD>={params.min_ad}" >> {log}
        
        # 修正：所有 FORMAT 字段使用 [样本:子字段] 格式
        bcftools view -f PASS {input.vcf} 2>> {log} \
        | bcftools filter \
            -i "FORMAT/AF[0:0] >= {params.min_af} && \
                FORMAT/DP[0:0] >= {params.min_depth} && \
                INFO/TLOD >= {params.min_tlod} && \
                INFO/MMQ >= {params.min_mq} && \
                FORMAT/AD[0:1] >= {params.min_ad}" \
            -Oz -o {output.vcf} \
            2>> {log}
        
        # 检查输出
        if [ -f {output.vcf} ]; then
            bcftools index -t {output.vcf} 2>> {log}
            echo "Output: $(bcftools view -H {output.vcf} 2>/dev/null | wc -l) variants" >> {log}
        else
            echo "WARNING: No output" >> {log}
        fi
        """