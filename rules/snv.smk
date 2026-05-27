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
        mem_mb=config["resources"]["mutect_memory_mb"]
    params:
        normal_name=NORMAL_SAMPLE,
        tmp_dir=TMP_DIR
    #conda:
    #    "envs/environment.yaml"
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
            --af-of-alleles-not-in-resource 0.0000025 \
            --minimum-allele-fraction 0.01 \
            --max-reads-per-alignment-start 0 \
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
    #conda:
    #    "envs/environment.yaml"
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
    基于实际可用的VCF字段
    """
    input:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.mutect2.filtered.vcf.gz",
        tbi=f"{VAR_DIR}/snv/{{tumor}}.mutect2.filtered.vcf.gz.tbi"
    output:
        vcf=f"{FINAL_DIR}/{{tumor}}.snv.filtered.vcf.gz",
        tbi=f"{FINAL_DIR}/{{tumor}}.snv.filtered.vcf.gz.tbi"
    log:
        f"{LOG_DIR}/ctdna_comprehensive_filter_{{tumor}}.log"
    params:
        # 核心参数（基于实际VCF字段）
        min_af=config["filtering"]["snv"].get("min_vaf", 0.005),
        min_tlod=config["filtering"]["snv"].get("min_tlod", 60),
        min_depth=config["filtering"]["snv"].get("min_depth", 10),
        min_mq=config["filtering"]["snv"].get("min_mq", 30),  # 对应 INFO/MMQ
        
        # 可选：AD过滤（支持突变reads数）
        min_ad=config["filtering"]["snv"].get("min_ad", 0),  # 新增，对应 FORMAT/AD[1]
        
        # 以下参数在标准Mutect2中不可用，建议设为0禁用
        #min_umi_families=config["filtering"]["snv"].get("min_umi_families", 0),
        #min_supporting_bases=config["filtering"]["snv"].get("min_supporting_bases", 0),
        #max_gnomad_af=config["filtering"]["snv"].get("max_gnomad_af", 1)  # 默认1表示不过滤
    threads: 2
    resources:
        mem_mb=4000
    shell:
        """
        # 构建基础过滤表达式
        FILTER_EXPR="FORMAT/AF < {params.min_af} || \
                     FORMAT/DP < {params.min_depth} || \
                     INFO/TLOD < {params.min_tlod} || \
                     INFO/MMQ < {params.min_mq}"
        
        # 添加AD过滤（如果需要）
        if [ {params.min_ad} -gt 0 ]; then
            FILTER_EXPR="$FILTER_EXPR || FORMAT/AD[:1] < {params.min_ad}"
        fi
        
        # 注意：STRAND_BIAS 已经由 GATK FilterMutectCalls 处理
        # 只需保留 PASS 位点即可，不需要额外过滤
        
        # 应用过滤
        bcftools view -f 'PASS' {input.vcf} | \
        bcftools filter \
            -e "$FILTER_EXPR" \
            -s "ctDNA_filter" \
            -Oz -o {output.vcf} 2> {log}
        
        bcftools index -t {output.vcf}
        """