
# ===================================================================
# CNV 检测（CNVkit 配对模式）
# ===================================================================

rule call_cnv_cnvkit:
    """
    使用 CNVkit 检测拷贝数变异
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        ref=REF_FASTA,
        bed="/home/zhangli_lab/zhouxiangyu/DATA/projects/CIBR-ZHANGLI/ctDNA/resource/target.bed"
    output:
        cnr=f"{VAR_DIR}/cnv/{{tumor}}.cnr",
        cns=f"{VAR_DIR}/cnv/{{tumor}}.cns"
    log:
        f"{LOG_DIR}/cnvkit_{{tumor}}.log"
    threads: 4
    resources:
        mem_mb=16000
    #conda:
    #    "envs/environment.yaml"
    shell:
        """
        ~/DATA/miniconda3/envs/cnvkit/bin/cnvkit.py batch \
            {input.tumor_bam} \
            --normal {input.normal_bam} \
            --targets {input.bed} \
            --fasta {input.ref} \
            --output-dir {VAR_DIR}/cnv/{wildcards.tumor}_cnvkit \
            --output-reference {VAR_DIR}/cnv/{wildcards.tumor}_reference.cnn \
            --diagram \
            --scatter \
            --processes {threads} \
            > {log} 2>&1
        
        # 复制标准命名的输出文件
        cp {VAR_DIR}/cnv/{wildcards.tumor}_cnvkit/{wildcards.tumor}.dedup.cnr {output.cnr}
        cp {VAR_DIR}/cnv/{wildcards.tumor}_cnvkit/{wildcards.tumor}.dedup.cns {output.cns}
        """



# ===================================================================
# 过滤
# ===================================================================


rule final_filter_cnv:
    """
    对 CNV 进行过滤
    """
    input:
        cns=f"{VAR_DIR}/cnv/{{tumor}}.cns",
        script="/home/zhangli_lab/zhouxiangyu/DATA/workflow/ctDNA-panel/scripts/filter_cnv.py"
    output:
        filtered=f"{FINAL_DIR}/{{tumor}}.cnv.filtered.cns"
    log:
        f"{LOG_DIR}/final_filter_cnv_{{tumor}}.log"
    params:
        min_log2_abs=config["filtering"]["cnv"]["min_log2_abs"],
        min_probes=config["filtering"]["cnv"]["min_probes"],
        max_ci=config["filtering"]["cnv"]["max_ci"]
    shell:
        """
        python3 {input.script} \
            --input {input.cns} \
            --output {output.filtered} \
            --min-log2-abs {params.min_log2_abs} \
            --min-probes {params.min_probes} \
            --max-ci {params.max_ci} \
            > {log} 2>&1
        """


