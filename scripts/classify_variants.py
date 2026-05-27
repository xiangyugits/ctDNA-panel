#!/usr/bin/env python3
"""
classify_variants.py - 根据临床意义对 ANNOVAR 注释的变异进行分级
分级标准：基于 ClinVar、COSMIC、功能预测等证据
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


class VariantClassifier:
    """变异临床分级器"""
    
    def __init__(self, tier1_terms=None, tier2_terms=None):
        self.tier1_terms = tier1_terms or [
            'Pathogenic', 'Likely_pathogenic',
            'drug_response', 'risk_factor'
        ]
        self.tier2_terms = tier2_terms or [
            'Uncertain_significance', 'conflicting'
        ]
        
        self.tiers = {
            'Tier1_actionable': [],    # 强临床意义，可操作
            'Tier2_potential': [],     # 潜在临床意义
            'Tier3_unknown': [],       # 意义不明确
            'Tier4_benign': []         # 良性或可能良性
        }
        
        self.stats = defaultdict(int)
    
    def _safe_float(self, value, default=0.0):
        """安全转换为浮点数"""
        try:
            if pd.isna(value) or value == '.' or value == '':
                return default
            return float(value)
        except (ValueError, TypeError):
            return default
    
    def _safe_str(self, value, default=''):
        """安全转换为字符串"""
        try:
            if pd.isna(value) or value == '.':
                return default
            return str(value)
        except (ValueError, TypeError):
            return default
    
    def _get_clinvar_significance(self, row):
        """获取 ClinVar 临床意义"""
        clinvar_cols = [col for col in row.index if 'clinvar' in col.lower()]
        for col in clinvar_cols:
            val = self._safe_str(row[col])
            if val:
                return val
        return ''
    
    def _get_cosmic_info(self, row):
        """获取 COSMIC 信息"""
        cosmic_cols = [col for col in row.index if 'cosmic' in col.lower()]
        for col in cosmic_cols:
            val = self._safe_str(row[col])
            if val and val != '.':
                return val
        return ''
    
    def _get_revel_score(self, row):
        """获取 REVEL 评分"""
        revel_cols = [col for col in row.index if 'revel' in col.lower()]
        for col in revel_cols:
            return self._safe_float(row[col])
        return 0.0
    
    def _get_splice_prediction(self, row):
        """获取剪接位点预测"""
        splice_cols = [col for col in row.index if 'dbscsnv' in col.lower()]
        predictions = []
        for col in splice_cols:
            val = self._safe_str(row[col])
            if val and val != '.':
                predictions.append(val)
        return ','.join(predictions) if predictions else ''
    
    def _is_pathogenic_clinvar(self, clinvar_sig):
        """判断 ClinVar 是否为致病"""
        if not clinvar_sig:
            return False
        clinvar_lower = clinvar_sig.lower()
        
        # 明确的致病变异
        pathogenic_terms = [
            'pathogenic', 'likely_pathogenic',
            'drug_response', 'risk_factor'
        ]
        if any(term in clinvar_lower for term in pathogenic_terms):
            if 'benign' not in clinvar_lower:
                return True
        return False
    
    def _is_benign_clinvar(self, clinvar_sig):
        """判断 ClinVar 是否为良性"""
        if not clinvar_sig:
            return False
        clinvar_lower = clinvar_sig.lower()
        
        benign_terms = ['benign', 'likely_benign']
        if any(term in clinvar_lower for term in benign_terms):
            if 'pathogenic' not in clinvar_lower:
                return True
        return False
    
    def _is_hotspot_mutation(self, cosmic_info, row):
        """判断是否为热点突变"""
        if cosmic_info:
            return True
        
        # 检查是否有 COSMIC ID
        for col in row.index:
            if 'cosmic' in col.lower() and 'id' in col.lower():
                val = self._safe_str(row[col])
                if val and val != '.':
                    return True
        return False
    
    def _is_damaging_prediction(self, row):
        """判断功能预测是否为有害"""
        damaging_scores = []
        
        # REVEL 评分 (>0.5 认为有害)
        revel = self._get_revel_score(row)
        if revel > 0.5:
            damaging_scores.append(f'REVEL={revel:.3f}')
        
        # SIFT 预测
        sift_cols = [col for col in row.index if 'sift' in col.lower() and 'pred' in col.lower()]
        for col in sift_cols:
            val = self._safe_str(row[col])
            if val.upper() in ['D', 'DAMAGING']:
                damaging_scores.append(f'SIFT={val}')
        
        # PolyPhen2 预测
        polyphen_cols = [col for col in row.index if 'polyphen' in col.lower() and 'pred' in col.lower()]
        for col in polyphen_cols:
            val = self._safe_str(row[col])
            if val.upper() in ['D', 'PROBABLY_DAMAGING', 'POSSIBLY_DAMAGING']:
                damaging_scores.append(f'PolyPhen2={val}')
        
        return damaging_scores
    
    def _is_function_important(self, row):
        """判断功能影响是否重要"""
        func_cols = [col for col in row.index if 'func' in col.lower() and 'refgene' in col.lower()]
        exonic_func_cols = [col for col in row.index if 'exonicfunc' in col.lower()]
        
        important_functions = [
            'stopgain', 'stoploss', 'frameshift',
            'nonsynonymous', 'nonframeshift',
            'splicing', 'splice'
        ]
        
        for col in func_cols + exonic_func_cols:
            val = self._safe_str(row[col]).lower()
            for func in important_functions:
                if func in val:
                    return True, val
        return False, ''
    
    def classify_variant(self, row):
        """
        对单个变异进行分级
        
        返回: (tier, reasons)
        """
        reasons = []
        
        # 收集证据
        clinvar_sig = self._get_clinvar_significance(row)
        cosmic_info = self._get_cosmic_info(row)
        damaging_preds = self._is_damaging_prediction(row)
        is_important, func_type = self._is_function_important(row)
        splice_pred = self._get_splice_prediction(row)
        
        # Tier 1: 强临床意义，可操作
        if self._is_pathogenic_clinvar(clinvar_sig):
            reasons.append(f'ClinVar: {clinvar_sig}')
            return 'Tier1_actionable', reasons
        
        if self._is_hotspot_mutation(cosmic_info, row):
            reasons.append(f'COSMIC hotspot: {cosmic_info}')
            if damaging_preds:
                reasons.append('Damaging: ' + '; '.join(damaging_preds))
            return 'Tier1_actionable', reasons
        
        # 功能丧失突变（移码、无义）+ 低人群频率
        if func_type in ['stopgain', 'frameshift', 'stoploss']:
            gnomad_cols = [col for col in row.index if 'gnomad' in col.lower() and 'af' in col.lower()]
            max_af = 0
            for col in gnomad_cols:
                max_af = max(max_af, self._safe_float(row[col]))
            
            if max_af < 0.01:
                reasons.append(f'LOF mutation ({func_type})')
                reasons.append(f'gnomAD AF={max_af:.4f}')
                return 'Tier1_actionable', reasons
        
        # 剪接位点突变
        if splice_pred and 'damaging' in splice_pred.lower():
            reasons.append(f'Splice site damaging: {splice_pred}')
            return 'Tier1_actionable', reasons
        
        # Tier 2: 潜在临床意义
        clinvar_lower = clinvar_sig.lower()
        if any(term.lower() in clinvar_lower for term in self.tier2_terms):
            reasons.append(f'ClinVar: {clinvar_sig}')
            return 'Tier2_potential', reasons
        
        if damaging_preds and is_important:
            reasons.append('Damaging: ' + '; '.join(damaging_preds))
            return 'Tier2_potential', reasons
        
        # Tier 4: 良性
        if self._is_benign_clinvar(clinvar_sig):
            reasons.append(f'ClinVar benign: {clinvar_sig}')
            return 'Tier4_benign', reasons
        
        # Tier 3: 意义不明确（默认）
        if is_important:
            reasons.append('Functional impact: ' + func_type)
        else:
            reasons.append('No clear evidence')
        
        return 'Tier3_unknown', reasons
    
    def classify_dataframe(self, df):
        """对整个数据框进行分类"""
        print(f"Classifying {len(df)} variants...")
        
        tier_names = []
        tier_reasons = []
        
        for idx, row in df.iterrows():
            tier, reasons = self.classify_variant(row)
            tier_names.append(tier)
            tier_reasons.append('; '.join(reasons))
            self.stats[tier] += 1
        
        df['TIER'] = tier_names
        df['TIER_REASONS'] = tier_reasons
        
        # 收集各分级
        for tier in ['Tier1_actionable', 'Tier2_potential', 'Tier3_unknown', 'Tier4_benign']:
            self.tiers[tier] = df[df['TIER'] == tier].to_dict('records')
        
        return df
    
    def write_tier_output(self, prefix):
        """输出分级结果"""
        print("\nWriting output files...")
        
        all_variants = []
        
        for tier in ['Tier1_actionable', 'Tier2_potential', 'Tier3_unknown', 'Tier4_benign']:
            variants = self.tiers[tier]
            if variants:
                tier_df = pd.DataFrame(variants)
                
                # 输出 TSV
                suffix = tier.replace('_', '_')
                output_file = f"{prefix}.{suffix}.txt"
                tier_df.to_csv(output_file, sep='\t', index=False)
                print(f"  {tier}: {len(variants)} variants -> {output_file}")
                
                # 收集所有变异
                for var in variants:
                    var['Final_Tier'] = tier
                    all_variants.append(var)
        
        # 输出 Excel（所有分级）
        if all_variants:
            all_df = pd.DataFrame(all_variants)
            excel_file = f"{prefix}.all_tiers.xlsx"
            
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # 按分级写入不同的 sheet
                for tier in ['Tier1_actionable', 'Tier2_potential', 
                            'Tier3_unknown', 'Tier4_benign']:
                    tier_df = all_df[all_df['Final_Tier'] == tier]
                    if not tier_df.empty:
                        # Excel sheet 名最长 31 字符
                        sheet_name = tier.replace('_', ' ')[:31]
                        tier_df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # 所有变异汇总
                all_df.to_excel(writer, sheet_name='All Variants', index=False)
            
            print(f"  All tiers: {len(all_variants)} variants -> {excel_file}")
    
    def generate_summary(self, prefix):
        """生成分级摘要"""
        summary_file = f"{prefix}.classification_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("Variant Clinical Classification Summary\n")
            f.write("=" * 70 + "\n\n")
            
            f.write("Tier Definitions:\n")
            f.write("-" * 70 + "\n")
            f.write("Tier 1 - Actionable: Strong clinical significance\n")
            f.write("  - Pathogenic/Likely pathogenic in ClinVar\n")
            f.write("  - COSMIC hotspot mutation\n")
            f.write("  - Loss-of-function mutation (frameshift, stop-gain)\n")
            f.write("  - Damaging splice site mutation\n\n")
            
            f.write("Tier 2 - Potential: Potential clinical significance\n")
            f.write("  - Uncertain significance in ClinVar\n")
            f.write("  - Damaging prediction + functional impact\n\n")
            
            f.write("Tier 3 - Unknown: Unknown significance\n")
            f.write("  - Functional impact but no clear evidence\n")
            f.write("  - No clinical database evidence\n\n")
            
            f.write("Tier 4 - Benign: Likely benign\n")
            f.write("  - Benign/Likely benign in ClinVar\n")
            f.write("  - High population frequency\n\n")
            
            f.write("=" * 70 + "\n")
            f.write("Classification Results:\n")
            f.write("-" * 70 + "\n")
            
            total = sum(self.stats.values())
            for tier in ['Tier1_actionable', 'Tier2_potential', 
                         'Tier3_unknown', 'Tier4_benign']:
                count = self.stats[tier]
                pct = count / max(total, 1) * 100
                f.write(f"  {tier:>20}: {count:>6} ({pct:>5.1f}%)\n")
            
            f.write(f"  {'─' * 35}\n")
            f.write(f"  {'Total':>20}: {total:>6}\n")
            f.write("=" * 70 + "\n\n")
            
            # Tier 1 详细列表
            if self.tiers['Tier1_actionable']:
                f.write("Tier 1 - Actionable Variants:\n")
                f.write("-" * 70 + "\n")
                for var in self.tiers['Tier1_actionable']:
                    gene = var.get('Gene.refGene', var.get('Gene', 'Unknown'))
                    func = var.get('ExonicFunc.refGene', var.get('Function', 'Unknown'))
                    clinvar = var.get('CLINSIG', var.get('clinvar', ''))
                    
                    # 尝试从各种可能的列获取位置信息
                    chrom = var.get('Chr', var.get('#Chr', var.get('CHROM', '?')))
                    start = var.get('Start', var.get('POS', '?'))
                    ref = var.get('Ref', var.get('REF', '?'))
                    alt = var.get('Alt', var.get('ALT', '?'))
                    
                    f.write(f"  {gene}: {chrom}:{start} {ref}>{alt}\n")
                    f.write(f"    Function: {func}\n")
                    if clinvar:
                        f.write(f"    ClinVar: {clinvar}\n")
                    reasons = var.get('TIER_REASONS', '')
                    if reasons:
                        f.write(f"    Reasons: {reasons}\n")
                    f.write("\n")
        
        print(f"  Summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Classify variants into clinical tiers based on ANNOVAR annotation"
    )
    
    parser.add_argument("--input", "-i", required=True,
                       help="Filtered ANNOVAR results file")
    parser.add_argument("--output-prefix", "-o", required=True,
                       help="Output file prefix")
    
    parser.add_argument("--tier1-terms",
                       default="Pathogenic,Likely_pathogenic,drug_response,risk_factor",
                       help="Comma-separated ClinVar terms for Tier 1")
    parser.add_argument("--tier2-terms",
                       default="Uncertain_significance,conflicting",
                       help="Comma-separated ClinVar terms for Tier 2")
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        return
    
    # 读取数据
    print(f"Loading variants from {args.input}")
    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    print(f"  Loaded {len(df)} variants")
    
    # 打印列名（调试用）
    print(f"  Available columns: {', '.join(df.columns[:20])}...")
    
    # 初始化分类器
    classifier = VariantClassifier(
        tier1_terms=args.tier1_terms.split(','),
        tier2_terms=args.tier2_terms.split(',')
    )
    
    # 执行分类
    df = classifier.classify_dataframe(df)
    
    # 输出结果
    classifier.write_tier_output(args.output_prefix)
    classifier.generate_summary(args.output_prefix)
    
    # 打印统计
    print("\n" + "=" * 50)
    print("Classification Complete")
    print("=" * 50)
    total = sum(classifier.stats.values())
    for tier in ['Tier1_actionable', 'Tier2_potential', 
                 'Tier3_unknown', 'Tier4_benign']:
        count = classifier.stats[tier]
        pct = count / max(total, 1) * 100
        print(f"  {tier}: {count} ({pct:.1f}%)")
    print(f"  Total: {total}")
    print("=" * 50)


if __name__ == "__main__":
    main()