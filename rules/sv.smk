FILTERING = config["filtering"]
SV_FILTER = FILTERING["sv"]
SV_DIR = f"{VAR_DIR}/sv"

# ============================================
# Rule 1: Delly 结构变异检测 (v0.8.x 语法)
# ============================================
rule delly_call:
    """
    使用 Delly v0.8.x 检测结构变异（原始BCF输出）
    语法: delly call -g ref -t SV_TYPE -o output.bcf -q MAPQ -s SUPPORT tumor.bam normal.bam
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        tumor_bai=f"{BAM_DIR}/{{tumor}}.dedup.bam.bai",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        normal_bai=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam.bai",
        ref=REF_FASTA
    output:
        bcf=f"{SV_DIR}/{{tumor}}.delly.raw.bcf"
    log:
        f"{LOG_DIR}/delly_call_{{tumor}}.log"
    threads: config["resources"]["delly_call_threads"]
    resources:
        mem_mb=config["resources"]["delly_call_memory"]
    params:
        # 使用逗号分隔，无空格
        sv_types=SV_FILTER.get("delly_sv_types", "DEL,DUP,INV,INS,TRA,BND"),
        min_mapq=SV_FILTER.get("delly_min_mapq", 10),
        min_support=SV_FILTER.get("delly_min_support", 3)
    shell:
        """
        module load delly
        
        # 打印过滤参数到日志
        echo "========================================" >> {log}
        echo "Delly Call Parameters:" >> {log}
        echo "========================================" >> {log}
        echo "Reference genome: {input.ref}" >> {log}
        echo "Tumor BAM: {input.tumor_bam}" >> {log}
        echo "Normal BAM: {input.normal_bam}" >> {log}
        echo "SV types: {params.sv_types}" >> {log}
        echo "Minimum MAPQ: {params.min_mapq}" >> {log}
        echo "Minimum support reads: {params.min_support}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "Memory: {resources.mem_mb} MB" >> {log}
        echo "========================================" >> {log}
        echo "" >> {log}
        
        # 执行Delly call
        delly call \
            -g {input.ref} \
            -t {params.sv_types} \
            -o {output.bcf} \
            -q {params.min_mapq} \
            -s {params.min_support} \
            {input.tumor_bam} {input.normal_bam} \
            2>> {log}
        
        # 打印结果统计
        echo "" >> {log}
        echo "========================================" >> {log}
        echo "Delly Call Results:" >> {log}
        echo "========================================" >> {log}
        if [ -f {output.bcf} ]; then
            echo "Output BCF file: {output.bcf}" >> {log}
            echo "BCF file size: $(du -h {output.bcf} | cut -f1)" >> {log}
            # 统计变异数量（如果bcftools可用）
            if command -v bcftools &> /dev/null; then
                SV_COUNT=$(bcftools view -H {output.bcf} 2>/dev/null | wc -l)
                echo "Number of raw SV calls: $SV_COUNT" >> {log}
            fi
        else
            echo "WARNING: Output BCF file not created!" >> {log}
        fi
        echo "========================================" >> {log}
        """


# ============================================
# Rule 2: Delly 体细胞变异过滤 (v0.8.x 语法 - 修正版)
# ============================================
rule delly_filter_somatic:
    """
    过滤体细胞结构变异 (v0.8.x)
    语法: delly filter -f somatic -o output.bcf -s tumor,normal input.bcf
    注意: -s 参数只能使用一次，样本用逗号分隔
    """
    input:
        bcf=rules.delly_call.output.bcf
    output:
        bcf=f"{SV_DIR}/{{tumor}}.delly.raw.bcf"
    log:
        f"{LOG_DIR}/delly_filter_{{tumor}}.log"
    threads: config["resources"]["delly_filter_threads"]
    resources:
        mem_mb=config["resources"]["delly_filter_memory"]
    params:
        # 使用逗号分隔的样本列表（不能有空格）
        samples=lambda wildcards: f"{wildcards.tumor},{NORMAL_SAMPLE}"
    shell:
        """
        module load delly
        
        delly filter \
            -f somatic \
            -o {output.bcf} \
            -s {params.samples} \
            {input.bcf} \
            2> {log}
        """


# ============================================
# Rule 3: BCF 转 VCF（全基因组）
# ============================================
rule delly_bcf_to_vcf:
    """
    将 BCF 格式转换为 VCF 格式（全基因组）
    """
    input:
        bcf=rules.delly_filter_somatic.output.bcf
    output:
        vcf=f"{SV_DIR}/{{tumor}}.delly.all.vcf.gz",
        tbi=f"{SV_DIR}/{{tumor}}.delly.all.vcf.gz.tbi"
    log:
        f"{LOG_DIR}/delly_bcf2vcf_{{tumor}}.log"
    threads: config["resources"]['default_threads']
    resources:
        mem_mb=config["resources"]['default_memory_mb']
    shell:
        """
        module load bcftools
        
        bcftools view -O z -o {output.vcf} {input.bcf}
        bcftools index {output.vcf}
        2> {log}
        """


# ============================================
# Rule 4: 过滤到靶向区域
# ============================================
rule delly_filter_targeted:
    """
    将全基因组 VCF 过滤到靶向区域
    """
    input:
        vcf=rules.delly_bcf_to_vcf.output.vcf,
        vcf_tbi=rules.delly_bcf_to_vcf.output.tbi,
        bed=TARGET_BED
    output:
        vcf=f"{SV_DIR}/{{tumor}}.delly.targeted.vcf.gz",
        tbi=f"{SV_DIR}/{{tumor}}.delly.targeted.vcf.gz.tbi",
        stats=f"{SV_DIR}/{{tumor}}.delly.targeted.stats.txt"
    log:
        f"{LOG_DIR}/delly_filter_targeted_{{tumor}}.log"
    threads: config["resources"]['default_threads']
    resources:
        mem_mb=config["resources"]['default_memory_mb']
    run:
        import subprocess
        import os
        
        # 检查BED文件是否存在且非空
        bed_file = input.bed
        if os.path.exists(bed_file) and os.path.getsize(bed_file) > 0:
            # 过滤到靶向区域
            cmd = f"""
            bcftools view -R {bed_file} {input.vcf} -O z -o {output.vcf}
            bcftools index {output.vcf}
            
            # 生成统计信息
            TOTAL=$(bcftools view -H {input.vcf} 2>/dev/null | wc -l)
            TARGETED=$(bcftools view -H {output.vcf} 2>/dev/null | wc -l)
            
            cat > {output.stats} << 'EOF'
Delly SV Statistics for {wildcards.tumor}
========================================
Total SVs (whole genome): $TOTAL
SVs in targeted regions: $TARGETED
Filtering BED file: {bed_file}
EOF
            """
            subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        else:
            # 如果没有BED文件，创建符号链接
            if not os.path.exists(output.vcf):
                os.symlink(input.vcf, output.vcf)
            if not os.path.exists(output.tbi):
                os.symlink(input.vcf_tbi, output.tbi)
            
            with open(output.stats, 'w') as f:
                f.write(f"No BED file provided. Using whole genome results for {wildcards.tumor}\n")


# ============================================
# Rule 5: 生成统计报告
# ============================================
rule delly_summarize:
    """
    生成 Delly 结果统计报告
    """
    input:
        vcf=rules.delly_filter_targeted.output.vcf,
        stats=rules.delly_filter_targeted.output.stats
    output:
        report=f"{SV_DIR}/{{tumor}}.delly.summary.html"
    log:
        f"{LOG_DIR}/delly_summary_{{tumor}}.log"
    threads: 1
    resources:
        mem_mb=2000
    run:
        import subprocess
        
        # 统计变异类型分布
        cmd = f"""
        # 统计变异类型
        if [ -s {input.vcf} ]; then
            bcftools query -f '%INFO/SVTYPE\\n' {input.vcf} 2>/dev/null | sort | uniq -c > {wildcards.tumor}.svtypes.tmp
        else
            touch {wildcards.tumor}.svtypes.tmp
        fi
        
        cat > {output.report} << 'EOF'
<html>
<head>
<title>Delly SV Report - {wildcards.tumor}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 20px; }}
h1 {{ color: #333; }}
h2 {{ color: #666; }}
pre {{ background-color: #f4f4f4; padding: 10px; border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>Delly Structural Variant Report</h1>
<h2>Sample: {wildcards.tumor}</h2>
<h2>Reference Genome: {GENOME}</h2>

<h3>Statistics</h3>
<pre>
$(cat {input.stats})
</pre>

<h3>SV Type Distribution</h3>
<pre>
$(cat {wildcards.tumor}.svtypes.tmp | awk '{{print $2": "$1}}')
</pre>

<h3>Configuration</h3>
<pre>
Min MAPQ: {SV_FILTER.get('delly_min_mapq', 10)}
Min Support: {SV_FILTER.get('delly_min_support', 3)}
SV Types: {SV_FILTER.get('delly_sv_types', 'DEL,DUP,INV,INS,TRA,BND')}
</pre>
</body>
</html>
EOF
        rm -f {wildcards.tumor}.svtypes.tmp
        """
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


# ============================================
# Rule 6: 清理临时文件
# ============================================
rule delly_cleanup:
    """
    清理 Delly 分析产生的临时文件
    """
    input:
        bcf_raw=rules.delly_call.output.bcf,
        bcf_somatic=rules.delly_filter_somatic.output.bcf
    output:
        done=f"{SV_DIR}/{{tumor}}.delly.cleanup.done"
    log:
        f"{LOG_DIR}/delly_cleanup_{{tumor}}.log"
    run:
        import os
        
        # 删除中间BCF文件（保留VCF）
        if os.path.exists(input.bcf_raw):
            os.remove(input.bcf_raw)
        if os.path.exists(input.bcf_somatic):
            os.remove(input.bcf_somatic)
        
        # 标记完成
        with open(output.done, 'w') as f:
            f.write(f"Cleanup completed for {wildcards.tumor}\n")


# ============================================
# Rule 7: 主规则（串联所有步骤）
# ============================================
rule call_sv_delly:
    """
    完整的 Delly 结构变异分析流程（全基因组 + 靶向过滤）
    """
    input:
        all_vcf=rules.delly_bcf_to_vcf.output.vcf,
        all_vcf_tbi=rules.delly_bcf_to_vcf.output.tbi,
        targeted_vcf=rules.delly_filter_targeted.output.vcf,
        targeted_vcf_tbi=rules.delly_filter_targeted.output.tbi,
        stats=rules.delly_filter_targeted.output.stats,
        report=rules.delly_summarize.output.report
    output:
        done=f"{SV_DIR}/{{tumor}}.delly.done"
    run:
        with open(output.done, 'w') as f:
            f.write(f"Delly analysis completed for {wildcards.tumor}\n")
            f.write(f"  - All VCF: {input.all_vcf}\n")
            f.write(f"  - Targeted VCF: {input.targeted_vcf}\n")
            f.write(f"  - Report: {input.report}\n")
            f.write(f"  - Stats: {input.stats}\n")

# ============================================
# Rule 7: manta 结构变异检测（备选方案）
# ============================================

rule call_sv_manta:
    """
    使用 Manta 检测结构变异（针对靶向测序优化）
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        tumor_bai=f"{BAM_DIR}/{{tumor}}.dedup.bam.bai",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        normal_bai=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam.bai",
        ref=REF_FASTA,
        bed=TARGET_BED
    output:
        vcf=f"{VAR_DIR}/sv/{{tumor}}.manta.vcf.gz",
        tbi=f"{VAR_DIR}/sv/{{tumor}}.manta.vcf.gz.tbi"
    log:
        f"{LOG_DIR}/manta_{{tumor}}.log"
    threads: 4
    resources:
        mem_mb=16000
    params:
        run_dir=f"{VAR_DIR}/sv/{{tumor}}_manta"
    shell:
        """
        module load manta
        
        # 关键修改1：不使用 --callRegions 和 --exome
        configManta.py \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam} \
            --referenceFasta {input.ref} \
            --runDir {params.run_dir} \
            > {log} 2>&1
        
        # 关键修改2：修改配置文件，降低阈值
        CONFIG_FILE="{params.run_dir}/configManta.py.ini"
        if [ -f "$CONFIG_FILE" ]; then
            # 降低高置信度reads阈值
            sed -i 's/minHqPairThreshold = 100/minHqPairThreshold = 10/g' $CONFIG_FILE
            sed -i 's/minHqMapq = 20/minHqMapq = 5/g' $CONFIG_FILE
            # 禁用统计检查（如果支持）
            echo "isSkipAlignmentStats = 1" >> $CONFIG_FILE
            echo "minCandidateRegionSize = 20" >> $CONFIG_FILE
            echo "maxCandidateRegionSize = 10000000" >> $CONFIG_FILE
        fi
        
        # 执行 Manta
        python2 {params.run_dir}/runWorkflow.py \
            -m local \
            -j {threads} \
            >> {log} 2>&1
        
        # 检查结果文件是否存在
        if [ -f {params.run_dir}/results/variants/somaticSV.vcf.gz ]; then
            cp {params.run_dir}/results/variants/somaticSV.vcf.gz {output.vcf}
            cp {params.run_dir}/results/variants/somaticSV.vcf.gz.tbi {output.tbi}
        else
            # 如果没有体细胞变异，创建空文件
            echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip > {output.vcf}
            tabix -p vcf {output.vcf}
        fi
        """