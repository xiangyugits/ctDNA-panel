# ===================================================================
# Part 8: ANNOVAR 变异注释
# ===================================================================

rule convert_vcf_to_annovar:
    """
    将 VCF 转换为 ANNOVAR 输入格式
    """
    input:
        vcf=f"{FINAL_DIR}/{{tumor}}.snv.filtered.vcf.gz"
    output:
        avinput=f"{FINAL_DIR}/{{tumor}}.snv.annovar_input"
    log:
        f"{LOG_DIR}/convert2annovar_{{tumor}}.log"
    params:
        annovar_dir=config["annotation"]["annovar"]["install_dir"]
    shell:
        """
        perl {params.annovar_dir}/convert2annovar.pl \
            -format vcf4 \
            -allsample \
            -withfreq \
            {input.vcf} \
            > {output.avinput} \
            2> {log}
        """

rule annotate_with_annovar:
    """
    使用 ANNOVAR 进行变异注释
    
    数据库：
    - refGene: 基因注释
    - avsnp147: dbSNP 147
    - clinvar: 临床意义 (20231231)
    - cosmic70: COSMIC v70 癌症突变
    - dbscsnv11: 剪接位点预测
    - revel: REVEL 致病性评分
    - gnomad_genome: gnomAD 全基因组频率
    """
    input:
        avinput=f"{FINAL_DIR}/{{tumor}}.snv.annovar_input"
    output:
        multianno=f"{FINAL_DIR}/{{tumor}}.snv.hg19_multianno.txt",
        vcf=f"{FINAL_DIR}/{{tumor}}.snv.annotated.vcf",
    log:
        f"{LOG_DIR}/annovar_{{tumor}}.log"
    params:
        annovar_dir=config["annotation"]["annovar"]["install_dir"],
        db_dir=config["annotation"]["annovar"]["db_dir"],
        buildver="hg19",
        # 协议：数据库名称
        protocols="refGene,avsnp147,clinvar,cosmic70,dbscsnv11,revel,gnomad_genome",
        # 操作：g=gene-based, f=filter-based
        operations="g,f,f,f,f,f,f"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        """
        perl {params.annovar_dir}/table_annovar.pl \
            {input.avinput} \
            {params.db_dir}/ \
            -buildver {params.buildver} \
            -out {FINAL_DIR}/{wildcards.tumor}.snv \
            -remove \
            -protocol {params.protocols} \
            -operation {params.operations} \
            -nastring . \
            -vcfinput \
            -thread {threads} \
            > {log} 2>&1
        """

rule filter_annovar_with_bcftools:
    """
    使用 bcftools 过滤 ANNOVAR 注释后的 VCF
    """
    input:
        vcf=f"{FINAL_DIR}/{{tumor}}.snv.annotated.vcf"
    output:
        filtered=f"{FINAL_DIR}/{{tumor}}.snv.annotated.filtered.vcf.gz",
        tbi=f"{FINAL_DIR}/{{tumor}}.snv.annotated.filtered.vcf.gz.tbi",
        summary=f"{FINAL_DIR}/{{tumor}}.snv.filtering_summary.txt"
    log:
        f"{LOG_DIR}/filter_annovar_{{tumor}}.log"
    params:
        # 过滤参数
        max_gnomad_af=config["filtering"]["snv"]["max_gnomad_af"],
        max_1000g_af=0.01,
        cancer_genes=config["annotation"]["post_annotation"]["gene_lists"]["cancer_genes"]
    shell:
        """
        # 1. 过滤人群频率
        # ANNOVAR 注释后的 INFO 字段包含 gnomAD_AF 等
        bcftools filter \
            -i 'INFO/gnomAD_AF <= {params.max_gnomad_af} || INFO/gnomAD_AF == "."' \
            {input.vcf} \
            2>> {log} \
        | bcftools filter \
            -i 'INFO/1000g2015aug_all <= {params.max_1000g_af} || INFO/1000g2015aug_all == "."' \
            2>> {log} \
        | bcftools filter \
            -i 'INFO/clinvar_20231231_clnsig !~ "Benign" || INFO/clinvar_20231231_clnsig == "."' \
            2>> {log} \
        | bcftools view \
            -i 'INFO/Func_refGene ~ "exonic" || INFO/Func_refGene ~ "splicing"' \
            2>> {log} \
        -Oz -o {output.filtered}
        
        # 2. 建立索引
        tabix -p vcf {output.filtered}
        
        # 3. 统计过滤结果
        echo "=== ANNOVAR bcftools Filtering Summary ===" > {output.summary}
        echo "Total variants in input:" >> {output.summary}
        bcftools stats {input.vcf} | grep "^SN" >> {output.summary}
        echo "" >> {output.summary}
        echo "Total variants after filtering:" >> {output.summary}
        bcftools stats {output.filtered} | grep "^SN" >> {output.summary}
        """

rule classify_variants_tier:
    """
    根据临床意义对变异进行分级
    """
    input:
        filtered=f"{FINAL_DIR}/{{tumor}}.snv.annotated.filtered.txt",
        script="/home/zhangli_lab/zhouxiangyu/DATA/workflow/ctDNA-panel/scripts/classify_variants.py"
    output:
        all_tiers=f"{FINAL_DIR}/{{tumor}}.snv.all_tiers.xlsx"
    log:
        f"{LOG_DIR}/classify_{{tumor}}.log"
    params:
        tier1_terms="Pathogenic,Likely_pathogenic,drug_response,risk_factor",
        tier2_terms="Uncertain_significance,conflicting"
    shell:
        """
        python3 {input.script} \
            --input {input.filtered} \
            --output-prefix {FINAL_DIR}/{wildcards.tumor}.snv \
            --tier1-terms {params.tier1_terms} \
            --tier2-terms {params.tier2_terms} \
            > {log} 2>&1
        """