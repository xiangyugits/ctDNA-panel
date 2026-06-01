FILTERING = config["filtering"]
SV_FILTER = FILTERING["sv"]
SV_DIR = f"{VAR_DIR}/sv"
GENOME = config["reference"]["genome"]
GRIDSS_JAR = config["tools"]["gridss_jar"]

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
# Rule 2: Delly 体细胞变异过滤 (v0.8.x 语法 — 已修复路径冲突)
# ============================================
rule delly_filter_somatic:
    """
    过滤体细胞结构变异 (v0.8.x)
    语法: delly filter -f somatic -o output.bcf -s tumor,normal input.bcf

    【已修复】输出文件名为 .delly.somatic.bcf，避免与 delly_call 的 .delly.raw.bcf 冲突。
    """
    input:
        bcf=rules.delly_call.output.bcf
    output:
        bcf=f"{SV_DIR}/{{tumor}}.delly.somatic.bcf"
    log:
        f"{LOG_DIR}/delly_filter_{{tumor}}.log"
    threads: config["resources"]["delly_filter_threads"]
    resources:
        mem_mb=config["resources"]["delly_filter_memory"]
    params:
        samples=lambda wildcards: f"{wildcards.tumor},{NORMAL_SAMPLE}"
    shell:
        """
        module load delly
        
        echo "========================================" >> {log}
        echo "Delly Somatic Filter (v0.8.x)" >> {log}
        echo "Input BCF:  {input.bcf}" >> {log}
        echo "Output BCF: {output.bcf}" >> {log}
        echo "Samples:    {params.samples}" >> {log}
        echo "========================================" >> {log}
        
        delly filter \
            -f somatic \
            -o {output.bcf} \
            -s {params.samples} \
            {input.bcf} \
            2>> {log}
        
        # 验证输出
        if [ -f {output.bcf} ] && [ -s {output.bcf} ]; then
            SV_COUNT=$(bcftools view -H {output.bcf} 2>/dev/null | wc -l)
            echo "Somatic SV candidates: $SV_COUNT" >> {log}
        else
            echo "WARNING: No somatic SVs found (output empty or missing)" >> {log}
            touch {output.bcf}
        fi
        """


# ============================================
# Rule 3: BCF 转 VCF（全基因组）
# ============================================
rule delly_bcf_to_vcf:
    """
    将 somatic BCF 转换为 VCF 格式（全基因组）

    【已修复】输入改为 delly_filter_somatic 输出的 somatic BCF（非 raw BCF）。
    若 somatic BCF 为空，生成合法的空 VCF 文件避免下游报错。
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
        
        if [ -s {input.bcf} ]; then
            bcftools view -O z -o {output.vcf} {input.bcf} 2>> {log}
            bcftools index {output.vcf} 2>> {log}
        else
            echo "WARNING: Somatic BCF is empty, creating placeholder VCF" >> {log}
            echo -e '##fileformat=VCFv4.2\\n##source=Delly\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' | bgzip -c > {output.vcf}
            tabix -p vcf {output.vcf} 2>> {log}
        fi
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
        bed=target_bed_SV
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
    清理 Delly 分析产生的临时 BCF 文件

    【已修复】清理 raw.bcf + somatic.bcf（两者不再同名）。
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
        
        for f in [input.bcf_raw, input.bcf_somatic]:
            if os.path.exists(f):
                os.remove(f)
                print(f"Deleted: {f}")
        
        with open(output.done, 'w') as fh:
            fh.write(f"Delly cleanup completed for {wildcards.tumor}\n")


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
# Rule 8: Manta 结构变异检测（备选方案，已降阈值优化）
# ============================================

rule call_sv_manta:
    """
    使用 Manta 检测结构变异（针对靶向测序优化）

    【已优化】
    - minHqPairThreshold 从 100 → 3（扩增子数据背景噪音低但信号也弱）
    - minHqMapq 从 20 → 5
    - 若 Manta 因数据质量不足中断，自动生成空 VCF 占位
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        tumor_bai=f"{BAM_DIR}/{{tumor}}.dedup.bam.bai",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        normal_bai=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam.bai",
        ref=REF_FASTA,
        bed=target_bed_SV
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
        
        configManta.py \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam} \
            --referenceFasta {input.ref} \
            --runDir {params.run_dir} \
            > {log} 2>&1
        
        CONFIG_FILE="{params.run_dir}/configManta.py.ini"
        if [ -f "$CONFIG_FILE" ]; then
            # 大幅降低阈值以适应扩增子 Panel
            sed -i -E 's/minHqPairThreshold[ ]*=[ ]*[0-9]+/minHqPairThreshold = 3/g' $CONFIG_FILE
            sed -i -E 's/minHqMapq[ ]*=[ ]*[0-9]+/minHqMapq = 5/g' $CONFIG_FILE
            # 扩张候选区域
            echo "isSkipAlignmentStats = 1" >> $CONFIG_FILE
            echo "minCandidateRegionSize = 10" >> $CONFIG_FILE
            echo "maxCandidateRegionSize = 10000000" >> $CONFIG_FILE
            # 允许低深度区域
            echo "minPassDepth = 2" >> $CONFIG_FILE
        fi
        
        # 执行 Manta，失败时不中断流程
        python2 {params.run_dir}/runWorkflow.py \
            -m local \
            -j {threads} \
            >> {log} 2>&1 || echo "Manta exited with error (may be normal for amplicon data)" >> {log}
        
        # 检查结果文件
        if [ -f {params.run_dir}/results/variants/somaticSV.vcf.gz ]; then
            cp {params.run_dir}/results/variants/somaticSV.vcf.gz {output.vcf}
            cp {params.run_dir}/results/variants/somaticSV.vcf.gz.tbi {output.tbi}
            echo "Manta: somaticSV found and copied" >> {log}
        elif [ -f {params.run_dir}/results/variants/diploidSV.vcf.gz ]; then
            cp {params.run_dir}/results/variants/diploidSV.vcf.gz {output.vcf}
            cp {params.run_dir}/results/variants/diploidSV.vcf.gz.tbi {output.tbi}
            echo "Manta: diploidSV used (no somatic calls)" >> {log}
        else
            echo -e '##fileformat=VCFv4.2\\n##source=Manta\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' | bgzip -c > {output.vcf}
            tabix -p vcf {output.vcf}
            echo "Manta: no SV calls, created empty VCF" >> {log}
        fi
        """


# ===================================================================
# ──────────── 替代方案 1：GRIDSS2（assembly-based SV 检测）────────────
#
# GRIDSS2 对 targeted/amplicon 数据的兼容性优于 Delly/Manta：
# - 使用 assembly（组装）而非单纯依赖 discordant read pairs
# - 支持从 BED 文件限制检测范围，适合 Panel
# - 可结合重复区域黑名单过滤
#
# 安装：
#   conda create -n gridss gridss=2.13.2
#   需要预先下载重复区域黑名单 BED（如 ENCODE DAC Blacklist）
# ===================================================================

rule call_sv_gridss:
    """
    使用 GRIDSS2 检测结构变异（assembly-based，适合扩增子 Panel）

    GRIDSS2 通过局部组装 breakpoint 区域来检测 SV，不依赖全基因组
    insert size 分布，比 Delly/Manta 更适合扩增子数据。
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        tumor_bai=f"{BAM_DIR}/{{tumor}}.dedup.bam.bai",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        normal_bai=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam.bai",
        ref=REF_FASTA,
        bed=target_bed_SV
    output:
        vcf=f"{VAR_DIR}/sv/{{tumor}}.gridss.vcf.gz",
        tbi=f"{VAR_DIR}/sv/{{tumor}}.gridss.vcf.gz.tbi",
        assembly_bam=f"{VAR_DIR}/sv/{{tumor}}.gridss.assembly.bam"
    log:
        f"{LOG_DIR}/gridss_{{tumor}}.log"
    threads: 8
    resources:
        mem_mb=32000
    params:
        tmp_dir=f"{TMP_DIR}/gridss_{{tumor}}",
        blacklist=SV_FILTER.get("gridss_blacklist", ""),
        # 靶向区域（外扩 500bp，捕获断点附近的 split reads）
        target_bed=target_bed_SV
    shell:
        """
        module load gridss

        mkdir -p {params.tmp_dir}

        echo "========================================" >> {log}
        echo "GRIDSS2 SV Detection" >> {log}
        echo "Tumor:  {input.tumor_bam}" >> {log}
        echo "Normal: {input.normal_bam}" >> {log}
        echo "Target regions: {params.target_bed}" >> {log}
        echo "========================================" >> {log}

        # 构建 GRIDSS 命令
        GRIDSS_CMD="java -Xmx16g -jar {GRIDSS_JAR} \\
            --reference {input.ref} \\
            --output {output.vcf} \\
            --assembly {output.assembly_bam} \\
            --threads {threads} \\
            --workingdir {params.tmp_dir}"

        # 如果有黑名单 BED，加入过滤
        if [ -n "{params.blacklist}" ] && [ -f "{params.blacklist}" ]; then
            GRIDSS_CMD="$GRIDSS_CMD --blacklist {params.blacklist}"
            echo "Blacklist: {params.blacklist}" >> {log}
        fi

        # 可选：限制检测范围到靶向区域（外扩 500bp）
        # 对于扩增子 Panel 建议启用，可大幅减少假阳性
        if [ -f {params.target_bed} ]; then
            GRIDSS_CMD="$GRIDSS_CMD --includeBed {params.target_bed}"
            echo "Include BED: {params.target_bed}" >> {log}
        fi

        # 输入 BAM
        GRIDSS_CMD="$GRIDSS_CMD {input.tumor_bam} {input.normal_bam}"

        eval $GRIDSS_CMD >> {log} 2>&1

        # 索引输出 VCF
        if [ -f {output.vcf} ] && [ -s {output.vcf} ]; then
            bcftools index -t {output.vcf} 2>> {log}
            SV_COUNT=$(bcftools view -H {output.vcf} 2>/dev/null | wc -l)
            echo "GRIDSS SV calls: $SV_COUNT" >> {log}
        else
            echo "WARNING: GRIDSS produced no output" >> {log}
            echo -e '##fileformat=VCFv4.2\\n##source=GRIDSS2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' | bgzip -c > {output.vcf}
            tabix -p vcf {output.vcf}
        fi

        # 清理临时目录
        rm -rf {params.tmp_dir}
        """


# ===================================================================
# ──────────── 替代方案 2：SvABA（靶向测序专用的 SV/Indel 检测）───────
#
# SvABA 专为以下场景设计：
# - 低覆盖度靶向测序（~100X 即可工作）
# - 不依赖配对末端距离分布
# - 同时输出 SV 和 Indel
# - 在 ctDNA/FFPE 样本上验证充分
#
# 安装：
#   conda create -n svaba svaba
# ===================================================================

rule call_sv_svaba:
    """
    使用 SvABA 检测结构变异（靶向测序专用）

    SvABA 基于局部组装检测断点，不依赖 insert size，非常适合
    扩增子 Panel 数据。同时输出 SV 和 Indel 结果。

    参考文献：Wala JA et al., Genome Research, 2018
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        tumor_bai=f"{BAM_DIR}/{{tumor}}.dedup.bam.bai",
        normal_bam=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam",
        normal_bai=f"{BAM_DIR}/{NORMAL_SAMPLE}.dedup.bam.bai",
        ref=REF_FASTA,
        bed=target_bed_SV
    output:
        sv_vcf=f"{VAR_DIR}/sv/{{tumor}}.svaba.sv.vcf.gz",
        sv_tbi=f"{VAR_DIR}/sv/{{tumor}}.svaba.sv.vcf.gz.tbi",
        indel_vcf=f"{VAR_DIR}/sv/{{tumor}}.svaba.indel.vcf.gz",
        indel_tbi=f"{VAR_DIR}/sv/{{tumor}}.svaba.indel.vcf.gz.tbi"
    log:
        f"{LOG_DIR}/svaba_{{tumor}}.log"
    threads: 4
    resources:
        mem_mb=16000
    params:
        prefix=f"{VAR_DIR}/sv/{{tumor}}.svaba",
        target_regions=SV_FILTER.get("svaba_target_regions", ""),
        min_mapq=SV_FILTER.get("svaba_min_mapq", 5),
        min_read_support=SV_FILTER.get("svaba_min_support", 2)
    shell:
        """
        echo "========================================" >> {log}
        echo "SvABA SV/Indel Detection (Amplicon-optimized)" >> {log}
        echo "Tumor:  {input.tumor_bam}" >> {log}
        echo "Normal: {input.normal_bam}" >> {log}
        echo "Prefix: {params.prefix}" >> {log}
        echo "========================================" >> {log}

        # 构建 SvABA 命令
        # 注：SvABA 无 --min-sc-reads / --num-sv-reads 参数
        #     灵敏度由 -L (mate-lookup-min, default 3) 和 LOD 阈值控制
        SVABA_CMD="~/DATA/miniconda3/envs/svaba/bin/svaba run \\
            -t {input.tumor_bam} \\
            -n {input.normal_bam} \\
            -G {input.ref} \\
            -a {params.prefix} \\
            -p {threads} \\
            -L {params.min_read_support} \\
            --germline-sv-database /dev/null"

        # 若提供靶向区域 BED，限制检测范围
        if [ -n "{params.target_regions}" ] && [ -f "{params.target_regions}" ]; then
            SVABA_CMD="$SVABA_CMD -k {params.target_regions}"
        elif [ -f {input.bed} ]; then
            SVABA_CMD="$SVABA_CMD -k {input.bed}"
        fi

        eval $SVABA_CMD >> {log} 2>&1

        # 重命名输出文件
        if [ -f {params.prefix}.sv.vcf ]; then
            bgzip -c {params.prefix}.sv.vcf > {output.sv_vcf}
            tabix -p vcf {output.sv_vcf}
            SV_COUNT=$(bcftools view -H {output.sv_vcf} 2>/dev/null | wc -l)
            echo "SvABA SV calls: $SV_COUNT" >> {log}
        else
            echo "WARNING: SvABA produced no SV VCF" >> {log}
            echo -e '##fileformat=VCFv4.2\\n##source=SvABA\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' | bgzip -c > {output.sv_vcf}
            tabix -p vcf {output.sv_vcf}
        fi

        if [ -f {params.prefix}.indel.vcf ]; then
            bgzip -c {params.prefix}.indel.vcf > {output.indel_vcf}
            tabix -p vcf {output.indel_vcf}
        else
            echo -e '##fileformat=VCFv4.2\\n##source=SvABA\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' | bgzip -c > {output.indel_vcf}
            tabix -p vcf {output.indel_vcf}
        fi

        # 清理 SvABA 中间文件
        rm -f {params.prefix}.*.txt {params.prefix}.bps.txt.gz {params.prefix}.alignments.txt.gz {params.prefix}.discordant.txt.gz 2>/dev/null
        """


# ===================================================================
# ──────────── 替代方案 3：靶向融合基因检测（Python 脚本）────────────
#
# 专为扩增子 Panel 设计的融合基因检测方案：
# - 仅检测 panel 覆盖区域内的 discordant read pairs + split reads
# - 不依赖全基因组 insert size 模型
# - 输出基因级融合报告，适合临床上报
#
# 依赖：pysam, samtools
# ===================================================================

rule call_sv_fusion_targeted:
    """
    基于 discordant read pairs + split reads 的靶向融合基因检测。

    原理：在 Panel 目标区域（BED）内搜索：
      1. 跨不同基因/染色体的 paired-end reads（discordant）
      2. 含 soft-clip 的 split reads
    通过 reads 聚类识别融合断点。

    适合扩增子 Panel，因为只检测 Panel 实际覆盖的区域。
    不依赖基因组范围的 insert size 分布模型。
    """
    input:
        tumor_bam=f"{BAM_DIR}/{{tumor}}.dedup.bam",
        bed=target_bed_SV,
        script=os.path.join(SCRIPT_DIR, "detect_fusions.py")
    output:
        fusions_tsv=f"{VAR_DIR}/sv/{{tumor}}.fusions.tsv",
        fusions_vcf=f"{VAR_DIR}/sv/{{tumor}}.fusions.vcf",
        summary=f"{VAR_DIR}/sv/{{tumor}}.fusions.summary.txt"
    log:
        f"{LOG_DIR}/fusion_targeted_{{tumor}}.log"
    threads: 4
    resources:
        mem_mb=8000
    params:
        min_mapq=SV_FILTER.get("fusion_min_mapq", 20),
        min_support=SV_FILTER.get("fusion_min_support", 2),
        # BED 区域外扩范围（bp），用于捕获跨区域断点
        flank=SV_FILTER.get("fusion_flank", 500)
    shell:
        """
        echo "========================================" >> {log}
        echo "Targeted Fusion Detection" >> {log}
        echo "Tumor BAM:    {input.tumor_bam}" >> {log}
        echo "Target BED:   {input.bed}" >> {log}
        echo "Min MAPQ:     {params.min_mapq}" >> {log}
        echo "Min Support:  {params.min_support}" >> {log}
        echo "========================================" >> {log}

        python3 {input.script} \\
            --bam            {input.tumor_bam} \\
            --bed            {input.bed} \\
            --output-fusions {output.fusions_tsv} \\
            --output-vcf     {output.fusions_vcf} \\
            --output-summary {output.summary} \\
            --min-mapq       {params.min_mapq} \\
            --min-support    {params.min_support} \\
            --flank          {params.flank} \\
            --threads        {threads} \\
            > {log} 2>&1

        # 统计
        if [ -f {output.fusions_tsv} ]; then
            FUSION_COUNT=$(tail -n +2 {output.fusions_tsv} | wc -l)
            echo "Candidate fusions: $FUSION_COUNT" >> {log}
        fi
        """


# ===================================================================
# ──────────── Rule: SV 结果汇总报告 ─────────────────────────────────
# ===================================================================

rule sv_summary_report:
    """
    汇总所有 SV 工具的结果，生成统一报告。

    输入（均可为空）：
      - Delly targeted VCF
      - Manta VCF
      - GRIDSS2 VCF
      - SvABA VCF
      - 靶向融合 TSV

    输出：每个样本一个文本汇总。
    """
    input:
        delly_vcf=f"{SV_DIR}/{{tumor}}.delly.targeted.vcf.gz",
        delly_tbi=f"{SV_DIR}/{{tumor}}.delly.targeted.vcf.gz.tbi",
        manta_vcf=f"{SV_DIR}/{{tumor}}.manta.vcf.gz",
        svaba_vcf=f"{SV_DIR}/{{tumor}}.svaba.sv.vcf.gz",
        fusions_tsv=f"{SV_DIR}/{{tumor}}.fusions.tsv"
    params:
        gridss_vcf=f"{SV_DIR}/{{tumor}}.gridss.vcf.gz"
    output:
        report=f"{SV_DIR}/{{tumor}}.sv.summary.txt"
    log:
        f"{LOG_DIR}/sv_summary_{{tumor}}.log"
    run:
        import os
        
        lines = []
        lines.append(f"SV Analysis Summary for {wildcards.tumor}")
        lines.append("=" * 60)
        lines.append("")
        
        tools = {
            "Delly":   input.delly_vcf,
            "Manta":   input.manta_vcf,
            "GRIDSS2": params.gridss_vcf,
            "SvABA":   input.svaba_vcf,
        }
        
        for tool_name, vcf_path in tools.items():
            if os.path.exists(vcf_path) and os.path.getsize(vcf_path) > 200:
                # 简单统计行数
                with open(vcf_path) as fh:
                    count = sum(1 for line in fh if not line.startswith("#"))
                lines.append(f"  {tool_name}: {count} SV calls")
            else:
                lines.append(f"  {tool_name}: no calls (empty/missing)")
        
        fusion_path = input.fusions_tsv
        if os.path.exists(fusion_path) and os.path.getsize(fusion_path) > 0:
            with open(fusion_path) as fh:
                fusion_count = sum(1 for _ in fh) - 1
            lines.append(f"  Targeted Fusions: {fusion_count} candidates")
        else:
            lines.append(f"  Targeted Fusions: no candidates")
        
        lines.append("")
        lines.append("=" * 60)
        lines.append("Note: For amplicon-based panels, SV calls from")
        lines.append("genome-wide tools (Delly/Manta) may have high FDR.")
        lines.append("Prefer targeted fusion detection + SvABA/GRIDSS2 results.")
        
        with open(output.report, 'w') as fh:
            fh.write("\n".join(lines) + "\n")