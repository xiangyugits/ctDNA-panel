# ===================================================================
# Part 9: 生成分析报告
# ===================================================================

rule generate_sample_report:
    """
    生成单个样本的分析报告
    """
    input:
        snv=f"{FINAL_DIR}/{{tumor}}.snv.filtered.vcf.gz",
        cnv=f"{FINAL_DIR}/{{tumor}}.cnv.filtered.cns",
        sv=f"{FINAL_DIR}/{{tumor}}.sv.filtered.vcf.gz",
        coverage=f"{QC_DIR}/{{tumor}}_coverage.txt",
        umi_stats=f"{UMI_DIR}/{{tumor}}_umi_stats.json",
        fastp_json=f"{QC_DIR}/{{tumor}}_fastp.json",
        script="scripts/generate_report.py"
    output:
        report=f"{REPORT_DIR}/{{tumor}}_analysis_report.html"
    log:
        f"{LOG_DIR}/report_{{tumor}}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        python3 {input.script} \
            --sample {wildcards.tumor} \
            --snv {input.snv} \
            --cnv {input.cnv} \
            --sv {input.sv} \
            --coverage {input.coverage} \
            --umi-stats {input.umi_stats} \
            --fastp-json {input.fastp_json} \
            --output {output.report} \
            > {log} 2>&1
        """

rule generate_cohort_summary:
    """生成样本队列汇总报告"""
    input:
        reports=expand(f"{REPORT_DIR}/{{tumor}}_analysis_report.html", tumor=TUMOR_SAMPLES)
    output:
        summary=f"{REPORT_DIR}/cohort_summary.html"
    log:
        f"{LOG_DIR}/cohort_summary.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        echo "<html><body>" > {output.summary}
        echo "<h1>ctDNA Analysis Cohort Summary</h1>" >> {output.summary}
        echo "<h2>Samples Analyzed</h2>" >> {output.summary}
        echo "<ul>" >> {output.summary}
        for tumor in {' '.join(TUMOR_SAMPLES)}; do
            echo "<li><a href='${{tumor}}_analysis_report.html'>${{tumor}}</a></li>" >> {output.summary}
        done
        echo "</ul>" >> {output.summary}
        echo "<h2>Analysis Date: $(date)</h2>" >> {output.summary}
        echo "</body></html>" >> {output.summary}
        """