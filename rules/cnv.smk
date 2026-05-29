SCRIPT_DIR = os.path.join(workflow.basedir, "scripts")

# ===================================================================
# CNV 检测（CNVkit 配对模式）
# 注意：CNVkit 适合杂交捕获（Capture-based）Panel；
#       扩增子（Amplicon-based）Panel 建议使用下方 panelcn.MOPS 规则。
# ===================================================================

rule call_cnv_cnvkit:
    """
    使用 CNVkit 检测拷贝数变异
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        ref=REF_FASTA,
        bed=config["reference"].get("target_bed_plain",config["reference"]["target_bed"])
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
        script=os.path.join(SCRIPT_DIR, "filter_cnv.py")
    output:
        filtered=f"{VAR_DIR}/cnv/{{tumor}}.cnv.filtered.cns"
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


# ===================================================================
# CNV 检测（panelcn.MOPS —— 扩增子 Panel 专用，推荐替代 CNVkit）
#
# 适用场景：
#   - Amplicon-based ctDNA Panel（无 off-target reads）
#   - 单配对正常样本或多 PoN 正常样本均支持
#
# 依赖环境（需预先在服务器安装）：
#   conda create -n panelcnmops r-base=4.3
#   conda activate panelcnmops
#   Rscript -e "BiocManager::install('panelcn.MOPS')"
#   conda install -n panelcnmops r-optparse
#
# 激活方式：
#   在 Snakefile rule all 中取消以下两行的注释：
#     expand(f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.gene.tsv", tumor=TUMOR_SAMPLES),
#     expand(f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.cns",      tumor=TUMOR_SAMPLES),
# ===================================================================

# ------------------------------------------------------------------
# Step A：运行 panelcn.MOPS R 脚本（per-amplicon CNV 检测）
# ------------------------------------------------------------------
rule call_cnv_panelcnmops:
    """
    使用 panelcn.MOPS 对扩增子 Panel 进行 CNV 检测（推荐用于 Amplicon-based Panel）
    """
    input:
        tumor_bam  = f"{BAM_DIR}/{{tumor}}.dedup.bam",
        tumor_bai  = f"{BAM_DIR}/{{tumor}}.dedup.bam.bai",
        normal_bam = f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        normal_bai = f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam.bai",
        rscript    = os.path.join(SCRIPT_DIR, "run_panelcnmops.R"),
        bed        = TARGET_BED
    output:
        tsv_raw    = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}/{{tumor}}.panelcnmops.tsv",
        tsv_all    = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}/{{tumor}}.panelcnmops.all.tsv",
        rds        = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}/{{tumor}}.panelcnmops.rds"
    log:
        f"{LOG_DIR}/panelcnmops_call_{{tumor}}.log"
    params:
        outdir     = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}",
        sample     = "{tumor}",
        genome     = config["reference"]["genome"],
        min_rc     = config["filtering"]["cnv"].get("panelcnmops_min_rc", 1),
        alpha      = config["filtering"]["cnv"].get("panelcnmops_alpha", 0.05),
        threads    = 4,
        rbin       = config.get("tools", {}).get("rscript",
                         "~/DATA/miniconda3/envs/panelcnmops/bin/Rscript")
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        mkdir -p {params.outdir}
        {params.rbin} {input.rscript} \
            --tumor-bam    {input.tumor_bam} \
            --normal-bam   {input.normal_bam} \
            --target-bed   {input.bed} \
            --sample-name  {params.sample} \
            --output-dir   {params.outdir} \
            --genome       {params.genome} \
            --min-rc       {params.min_rc} \
            --alpha        {params.alpha} \
            --threads      {params.threads} \
            --log          {log} \
            2>> {log}
        """


# ------------------------------------------------------------------
# Step B：解析 panelcn.MOPS 结果，生成基因级报告 + .cns 兼容格式
# ------------------------------------------------------------------
rule parse_cnv_panelcnmops:
    """
    将 panelcn.MOPS per-amplicon 结果聚合为基因级 CNV 报告，
    并输出 CNVkit 兼容的 .cns 格式（可继续用 filter_cnv.py 过滤）
    """
    input:
        tsv    = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}/{{tumor}}.panelcnmops.tsv",
        script = os.path.join(SCRIPT_DIR, "parse_panelcnmops.py")
    output:
        gene_tsv = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.gene.tsv",
        cns      = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.cns",
        report   = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.report.txt"
    log:
        f"{LOG_DIR}/panelcnmops_parse_{{tumor}}.log"
    params:
        sample       = "{tumor}",
        outdir       = f"{VAR_DIR}/cnv/panelcnmops",
        min_log2_abs = config["filtering"]["cnv"].get("panelcnmops_min_log2_abs",
                           config["filtering"]["cnv"]["min_log2_abs"]),
        min_frac     = config["filtering"]["cnv"].get("panelcnmops_min_cnv_fraction", 0.3),
        min_pvalue   = config["filtering"]["cnv"].get("panelcnmops_min_pvalue", 0.05)
    shell:
        """
        python3 {input.script} \
            --input            {input.tsv} \
            --sample-name      {params.sample} \
            --output-dir       {params.outdir} \
            --min-log2-abs     {params.min_log2_abs} \
            --min-cnv-fraction {params.min_frac} \
            --min-pvalue       {params.min_pvalue} \
            > {log} 2>&1
        """


# ------------------------------------------------------------------
# Step C：对 panelcn.MOPS .cns 结果进行二次过滤（复用 filter_cnv.py）
# ------------------------------------------------------------------
rule filter_cnv_panelcnmops:
    """
    对 panelcn.MOPS 的 .cns 结果做最终过滤，输出最终报告级 CNV
    """
    input:
        cns    = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.cns",
        script = os.path.join(SCRIPT_DIR, "filter_cnv.py")
    output:
        filtered = f"{VAR_DIR}/cnv/panelcnmops/{{tumor}}.panelcnmops.filtered.cns"
    log:
        f"{LOG_DIR}/panelcnmops_filter_{{tumor}}.log"
    params:
        min_log2_abs = config["filtering"]["cnv"].get("panelcnmops_min_log2_abs",
                           config["filtering"]["cnv"]["min_log2_abs"]),
        min_probes   = config["filtering"]["cnv"].get("panelcnmops_min_amplicons",
                           config["filtering"]["cnv"]["min_probes"]),
        max_ci       = config["filtering"]["cnv"].get("panelcnmops_max_ci",
                           config["filtering"]["cnv"]["max_ci"])
    shell:
        """
        python3 {input.script} \
            --input        {input.cns} \
            --output       {output.filtered} \
            --min-log2-abs {params.min_log2_abs} \
            --min-probes   {params.min_probes} \
            --max-ci       {params.max_ci} \
            > {log} 2>&1
        """


