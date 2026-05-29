SCRIPT_DIR = os.path.join(workflow.basedir, "scripts")

# ===================================================================
# Part 8: ANNOVAR 变异注释
# ===================================================================

rule convert_vcf_to_annovar:
    """
    将 VCF 转换为 ANNOVAR 输入格式
    """
    input:
        vcf=f"{VAR_DIR}/snv/{{tumor}}.snv.filtered.vcf.gz",
    output:
        avinput=f"{VAR_DIR}/annovar/{{tumor}}.snv.annovar_input"
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
        avinput=f"{VAR_DIR}/annovar/{{tumor}}.snv.annovar_input"
    output:
        multianno=f"{VAR_DIR}/annovar/{{tumor}}.snv.hg19_multianno.txt",
        vcf=f"{VAR_DIR}/annovar/{{tumor}}.snv.hg19_multianno.vcf",
    log:
        f"{LOG_DIR}/annovar_{{tumor}}.log"
    params:
        annovar_dir=config["annotation"]["annovar"]["install_dir"],
        db_dir=config["annotation"]["annovar"]["db_dir"],
        out_prefix=f"{VAR_DIR}/annovar/{{tumor}}.snv",
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
            -out {params.out_prefix} \
            -protocol {params.protocols} \
            -operation {params.operations} \
            -nastring . \
            -vcfinput \
            -thread {threads} \
            > {log} 2>&1
        """
rule filter_by_genelist:
    """
    根据基因列表筛选 ANNOVAR 注释结果，并格式化输出
    输出格式：染色体、位置(3'rule)、基因(HGVS)、转录本(ClinVar)、
             cHGVS、pHGVS、类型(SNV/Insertion/Deletion/Complex)、VAF(两位小数)
    """
    input:
        multianno=f"{VAR_DIR}/annovar/{{tumor}}.snv.hg19_multianno.txt",
        gene_list=config["annotation"]["annovar"]["target_genes"]
    output:
        filtered=f"{FINAL_DIR}/{{tumor}}.snv.genelist_filtered.txt"
    log:
        f"{LOG_DIR}/filter_genelist_{{tumor}}.log"
    params:
        script=os.path.join(SCRIPT_DIR, "filter_by_genelist.py")
    shell:
        """
        python3 {params.script} \
            --multianno {input.multianno} \
            --gene-list {input.gene_list} \
            --output {output.filtered} \
            > {log} 2>&1
        """

rule filter_annovar_variants:
    """
    基于 ANNOVAR 注释结果过滤变异位点
    
    功能：
    - 保留功能性变异（exonic/splicing）
    - 过滤高人群频率变异（gnomAD EAS < 0.001）
    - REVEL 评分过滤（错义突变需要 ≥ 0.5）
    - 计算优先级分数（High/Medium/Low）
    - 只保留必要列
    """
    input:
        ann_txt=f"{VAR_DIR}/annovar/{{tumor}}.snv.hg19_multianno.txt"
    output:
        filtered=f"{FINAL_DIR}/{{tumor}}.snv.filtered.txt",
        summary=f"{FINAL_DIR}/{{tumor}}.snv.filter_summary.txt"
    log:
        f"{LOG_DIR}/filter_annovar_{{tumor}}.log"
    params:
        SampleName="{tumor}",
        # 保留的列
        keep_columns = [
            'Sample', 'Chr', 'Start', 'End', 'Ref', 'Alt',
            'Func.refGene', 'Gene', 'ExonicFunc.refGene', 'AAChange.refGene',
            'avsnp147', 'CLNSIG', 'cosmic70', 'REVEL', 'gnomAD_genome_EAS'
        ],
        # 有害功能效应
        deleterious_effects = [
            "nonsynonymous SNV", "stopgain", "stoploss",
            "frameshift insertion", "frameshift deletion",
            "nonframeshift insertion", "nonframeshift deletion", "startloss"
        ],
        # ClinVar 致病性关键词
        pathogenic_keywords = ['Pathogenic', 'Likely_pathogenic', 'risk_factor'],
        # 过滤阈值
        gnomad_threshold = config["filtering"]["snv"]["max_gnomad_af"],
        revel_threshold = 0.5
    run:
        import pandas as pd
        import numpy as np
        import warnings
        warnings.filterwarnings('ignore')
        
        # ============================================================
        # 辅助函数
        # ============================================================
        def get_float(val, default=np.nan):
            """安全转换为浮点数"""
            if pd.isna(val) or val == '.' or val == '':
                return default
            try:
                return float(val)
            except:
                return default
        
        def evaluate_variant(row):
            """评估变异，返回 -1(拒绝) 或 优先级分数(0-6)"""
            
            # 1. 功能区域检查
            if row['Func.refGene'] not in ['exonic', 'splicing', 'exonic;splicing']:
                return -1
            
            # 2. 功能效应检查
            if row['Func.refGene'] == 'exonic':
                if row['ExonicFunc.refGene'] not in params.deleterious_effects:
                    return -1
            
            # 3. 人群频率检查
            freq = get_float(row['gnomAD_genome_EAS'])
            if pd.notna(freq) and freq >= params.gnomad_threshold:
                return -1
            
            # 4. REVEL 检查（仅错义突变）
            if row['ExonicFunc.refGene'] == 'nonsynonymous SNV':
                revel = get_float(row['REVEL'])
                if pd.notna(revel) and revel < params.revel_threshold:
                    return -1
            
            # 5. 计算优先级分数
            priority = 0
            
            # ClinVar 致病性
            clnsig = str(row.get('CLNSIG', ''))
            if 'Pathogenic' in clnsig and 'Likely' not in clnsig:
                priority += 3
            elif any(kw in clnsig for kw in params.pathogenic_keywords):
                priority += 2
            
            # COSMIC 记录
            if pd.notna(row.get('cosmic70')) and row['cosmic70'] not in ['.', '']:
                priority += 2
            
            # REVEL 高分
            revel = get_float(row.get('REVEL'))
            if pd.notna(revel):
                priority += 2 if revel >= 0.75 else 1 if revel >= 0.5 else 0
            
            # 严重功能效应加分
            if row['ExonicFunc.refGene'] in ['stopgain', 'stoploss', 
                                              'frameshift insertion', 'frameshift deletion']:
                priority += 2
            
            # Novel 变异奖励（gnomAD 缺失）
            if pd.isna(get_float(row['gnomAD_genome_EAS'])):
                priority += 1
            
            return priority
        
        # ============================================================
        # 主处理逻辑
        # ============================================================
        
        # 读取数据
        df = pd.read_csv(input.ann_txt, sep='\t', low_memory=False)
        total_records = len(df)
        df.loc[:, 'Sample'] = params.SampleName  # 添加样本名称列
        
        # 保留需要的列
        available_cols = [col for col in params.keep_columns if col in df.columns]
        df = df[available_cols].copy()
        
        # 应用评估
        df['eval'] = df.apply(evaluate_variant, axis=1)
        
        # 分离通过/拒绝
        passed = df[df['eval'] != -1].copy()
        rejected = df[df['eval'] == -1].copy()
        
        # 添加优先级标签
        def get_priority(score):
            if score >= 3:
                return 'High'
            elif score >= 1:
                return 'Medium'
            else:
                return 'Low'
        
        passed['Priority'] = passed['eval'].apply(get_priority)
        
        # 按优先级排序
        passed = passed.sort_values('eval', ascending=False)
        rejected_sorted = rejected.sort_values('eval', ascending=False)
        
        # 合并（拒绝的排在后面）
        final_df = pd.concat([passed, rejected_sorted], ignore_index=True)
        
        # 删除临时列
        #final_df = final_df.drop(columns=['eval'])
        
        # 保存结果
        final_df.to_csv(output.filtered, sep='\t', index=False)
        
        # ============================================================
        # 生成统计报告
        # ============================================================
        with open(output.summary, 'w') as f:
            f.write("="*60 + "\n")
            f.write(f"变异过滤统计报告 - {wildcards.tumor}\n")
            f.write("="*60 + "\n\n")
            
            f.write(f"总变异数: {total_records}\n")
            f.write(f"通过过滤: {len(passed)} ({len(passed)/total_records*100:.1f}%)\n")
            f.write(f"拒绝过滤: {len(rejected)} ({len(rejected)/total_records*100:.1f}%)\n\n")
            
            f.write("优先级分布:\n")
            if len(passed) > 0:
                high_cnt = (passed['Priority'] == 'High').sum()
                medium_cnt = (passed['Priority'] == 'Medium').sum()
                low_cnt = (passed['Priority'] == 'Low').sum()
                f.write(f"  High: {high_cnt} ({high_cnt/len(passed)*100:.1f}%)\n")
                f.write(f"  Medium: {medium_cnt} ({medium_cnt/len(passed)*100:.1f}%)\n")
                f.write(f"  Low: {low_cnt} ({low_cnt/len(passed)*100:.1f}%)\n")
            else:
                f.write("  无通过变异的变异\n")
            
            f.write(f"\n过滤参数:\n")
            f.write(f"  GNOMAD_THRESHOLD: {params.gnomad_threshold}\n")
            f.write(f"  REVEL_THRESHOLD: {params.revel_threshold}\n")
            f.write(f"  保留列数: {len(available_cols)}\n")
        
        # 记录日志
        with open(log[0], 'w') as f:
            f.write(f"过滤完成时间: {pd.Timestamp.now()}\n")
            f.write(f"输入文件: {input.ann_txt}\n")
            f.write(f"输出文件: {output.filtered}\n")
            f.write(f"通过变异数: {len(passed)}\n")

# ============================================================
# 合并所有样本的过滤结果
# ============================================================
rule merge_filtered_variants:
    """
    合并所有样本的过滤结果
    """
    input:
        expand(f"{FINAL_DIR}/filtered/{{tumor}}.snv.filtered.txt", tumor=config['samples'])
    output:
        merged=f"{FINAL_DIR}/filtered/all_samples.filtered.xlsx"
    run:
        import pandas as pd
        
        all_dfs = []
        for f in input:
            df = pd.read_csv(f, sep='\t')
            all_dfs.append(df)
        
        merged_df = pd.concat(all_dfs, ignore_index=True)
        
        # 保存为 Excel（便于查看）
        with pd.ExcelWriter(output.merged) as writer:
            merged_df.to_excel(writer, sheet_name='All_Variants', index=False)
            
            # 添加优先级统计 sheet
            priority_stats = merged_df['Priority'].value_counts().to_frame()
            priority_stats.to_excel(writer, sheet_name='Priority_Stats')
