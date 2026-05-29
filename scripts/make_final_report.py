#!/usr/bin/env python3
"""
make_final_report.py - 整合标准格式变异结果与 Tier 分级，生成最终临床报告

输入：
  1. filter_by_genelist.py 产出的标准格式 TSV（12 列）
  2. classify_variants.py 产出的分级结果目录

输出：
  Excel 工作簿，包含：
    - 标准格式 sheet（所有变异，含 Tier 列）
    - Tier1_Actionable sheet
    - Tier2_Potential sheet  
    - Tier3_Unknown sheet
    - Tier4_Benign sheet
    - 统计摘要 sheet
"""

import pandas as pd
import argparse
import sys
import re
from pathlib import Path
from typing import Optional, Tuple


# ============================================================
# 配置
# ============================================================

TIER_LABELS = {
    'Tier1_actionable': 'Tier 1 - 强临床意义（可操作）',
    'Tier2_potential':  'Tier 2 - 潜在临床意义',
    'Tier3_unknown':    'Tier 3 - 临床意义不明',
    'Tier4_benign':     'Tier 4 - 良性/可能良性'
}

TIER_COLORS = {
    'Tier1_actionable': '#FF6B6B',  # 红色
    'Tier2_potential':  '#FFD93D',  # 黄色
    'Tier3_unknown':    '#6BCB77',  # 绿色
    'Tier4_benign':     '#4D96FF'   # 蓝色
}


def load_classification_results(prefix_dir: str, sample_name: str) -> Optional[pd.DataFrame]:
    """
    从 classify_variants.py 输出中加载 Tier 分级结果。

    查找文件：{prefix_dir}/{sample_name}.classification.all_tiers.xlsx
    或：{prefix_dir}/all_tiers.xlsx
    """
    candidates = [
        Path(prefix_dir) / f"{sample_name}.classification.all_tiers.xlsx",
        Path(prefix_dir) / "all_tiers.xlsx",
        Path(prefix_dir) / f"{sample_name}.all_tiers.xlsx",
    ]

    # 同时搜索 *_all_tiers.xlsx
    for p in Path(prefix_dir).glob("*_all_tiers.xlsx"):
        candidates.append(p)

    for path in candidates:
        if path.exists():
            print(f"[INFO] Loading Tier classification from: {path}")
            try:
                df = pd.read_excel(path, sheet_name='All Variants')
                print(f"[INFO]   Loaded {len(df)} classified variants")
                return df
            except Exception as e:
                print(f"[WARN] Failed to read {path}: {e}")
                continue

    print("[WARN] No Tier classification file found.")
    return None


def load_tier_tsv_files(tier_dir: str) -> Optional[pd.DataFrame]:
    """
    从单个 TSV 分级文件加载 Tier 信息。
    文件格式：{tier_dir}/*.Tier1_actionable.txt 等
    """
    all_rows = []

    tier_files = {
        'Tier1_actionable': 'Tier1_actionable',
        'Tier2_potential': 'Tier2_potential',
        'Tier3_unknown': 'Tier3_unknown',
        'Tier4_benign': 'Tier4_benign'
    }

    for tier_name, suffix in tier_files.items():
        pattern = f"*.{suffix}.txt"
        matches = list(Path(tier_dir).glob(pattern))
        for m in matches:
            try:
                tdf = pd.read_csv(m, sep='\t', low_memory=False)
                tdf['TIER'] = tier_name
                all_rows.append(tdf)
                print(f"[INFO]   {tier_name}: {len(tdf)} variants from {m.name}")
            except Exception as e:
                print(f"[WARN] Failed to read {m}: {e}")

    if all_rows:
        return pd.concat(all_rows, ignore_index=True)
    return None


def build_variant_key(df: pd.DataFrame) -> pd.DataFrame:
    """
    构建用于匹配 Tier 的复合键。
    匹配逻辑：Chr + Start(原始) + Ref + Alt + Gene
    """
    df = df.copy()

    def make_key(row):
        chrom = str(row.get('Chr', '')).replace('chr', '').replace('Chr', '')
        pos = str(row.get('Start', ''))
        ref = str(row.get('Ref', ''))
        alt = str(row.get('Alt', ''))
        gene = str(row.get('Gene', ''))
        return f"{chrom}:{pos}:{ref}>{alt}:{gene}"

    df['_variant_key'] = df.apply(make_key, axis=1)
    return df


def merge_with_tiers(
    standard_df: pd.DataFrame,
    tier_df: Optional[pd.DataFrame]
) -> pd.DataFrame:
    """
    将 Tier 分级信息合并到标准格式数据框中。
    """
    if tier_df is None or tier_df.empty:
        standard_df['Tier'] = 'Unclassified'
        standard_df['Tier_Reason'] = ''
        return standard_df

    # 为 tier_df 构建匹配键
    tier_df = build_variant_key(tier_df)
    standard_df = build_variant_key(standard_df)

    # 尝试多级匹配
    # Level 1: 精确位置 + 基因
    tier_map = {}
    tier_reason_map = {}

    for _, row in tier_df.iterrows():
        key = row.get('_variant_key', '')
        tier = row.get('TIER', row.get('Final_Tier', 'Unclassified'))
        reason = row.get('TIER_REASONS', row.get('TIER_REASON', ''))

        if key and key not in tier_map:
            tier_map[key] = tier
            tier_reason_map[key] = reason

    standard_df['Tier'] = standard_df['_variant_key'].map(tier_map).fillna('Unclassified')
    standard_df['Tier_Reason'] = standard_df['_variant_key'].map(tier_reason_map).fillna('')

    # Level 2: 基因级备选匹配
    unmatched = standard_df[standard_df['Tier'] == 'Unclassified']
    if len(unmatched) > 0:
        gene_tier_map = {}
        for _, row in tier_df.iterrows():
            gene = str(row.get('Gene.refGene', row.get('Gene', ''))).split(',')[0].strip()
            tier = row.get('TIER', row.get('Final_Tier', ''))
            if gene and gene not in gene_tier_map:
                gene_tier_map[gene] = tier

        for idx in unmatched.index:
            gene = standard_df.loc[idx, 'Gene']
            if gene in gene_tier_map:
                standard_df.loc[idx, 'Tier'] = gene_tier_map[gene]
                standard_df.loc[idx, 'Tier_Reason'] = '(gene-level match)'

    # 清理辅助列
    standard_df = standard_df.drop(columns=['_variant_key'], errors='ignore')

    return standard_df


def write_final_excel(
    df: pd.DataFrame,
    output_file: str,
    sample_name: str
):
    """
    写入最终 Excel 报告，包含多个 sheet。
    """
    # 标准格式列
    standard_cols = [
        "样本编号", "Chr", "Start", "End", "Ref", "Alt",
        "Gene", "Type", "Transcript", "cHGVS", "pHGVS", "VAF",
        "Tier", "Tier_Reason"
    ]

    # 确保所有列存在
    for col in standard_cols:
        if col not in df.columns:
            df[col] = '.'

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # ===== Sheet 1: 标准格式（全部变异） =====
        all_df = df[standard_cols].copy()
        all_df.to_excel(writer, sheet_name='Standard_Report', index=False)

        # ===== Sheet 2-5: 按 Tier 分 sheet =====
        tier_order = [
            'Tier1_actionable',
            'Tier2_potential',
            'Tier3_unknown',
            'Tier4_benign'
        ]

        for tier in tier_order:
            tier_df = df[df['Tier'] == tier][standard_cols].copy()
            sheet_name = tier.replace('_', ' ')[:31]  # Excel sheet name limit
            if len(tier_df) > 0:
                tier_df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"[INFO]   {tier}: {len(tier_df)} variants")

        # Unclassified
        unclass_df = df[df['Tier'] == 'Unclassified'][standard_cols].copy()
        if len(unclass_df) > 0:
            unclass_df.to_excel(writer, sheet_name='Unclassified', index=False)
            print(f"[INFO]   Unclassified: {len(unclass_df)} variants")

        # ===== Sheet 最后: 统计摘要 =====
        summary_data = _build_summary(df, sample_name)
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)

    print(f"\n[INFO] Final report written: {output_file}")


def _build_summary(df: pd.DataFrame, sample_name: str) -> list:
    """构建统计摘要数据"""
    rows = []

    rows.append({'Category': 'Sample', 'Value': sample_name})
    rows.append({'Category': '', 'Value': ''})

    rows.append({'Category': '=== Variant Summary ===', 'Value': ''})
    rows.append({'Category': 'Total Variants', 'Value': len(df)})

    # 按类型
    for vtype, count in df['Type'].value_counts().items():
        rows.append({'Category': f'  {vtype}', 'Value': count})

    rows.append({'Category': '', 'Value': ''})
    rows.append({'Category': '=== Tier Classification ===', 'Value': ''})

    tier_order = ['Tier1_actionable', 'Tier2_potential', 'Tier3_unknown', 'Tier4_benign']
    for tier in tier_order:
        count = (df['Tier'] == tier).sum()
        pct = count / max(len(df), 1) * 100
        label = TIER_LABELS.get(tier, tier)
        rows.append({'Category': f'  {label}', 'Value': f'{count} ({pct:.1f}%)'})

    unclass = (df['Tier'] == 'Unclassified').sum()
    if unclass > 0:
        rows.append({'Category': '  Unclassified', 'Value': unclass})

    rows.append({'Category': '', 'Value': ''})
    rows.append({'Category': '=== Gene Distribution ===', 'Value': ''})

    for gene, count in df['Gene'].value_counts().head(20).items():
        rows.append({'Category': f'  {gene}', 'Value': count})

    # VAF 统计
    rows.append({'Category': '', 'Value': ''})
    rows.append({'Category': '=== VAF Statistics ===', 'Value': ''})

    vafs = []
    for v in df['VAF']:
        try:
            vafs.append(float(v))
        except (ValueError, TypeError):
            pass

    if vafs:
        rows.append({'Category': '  Min', 'Value': f'{min(vafs):.4f}'})
        rows.append({'Category': '  Max', 'Value': f'{max(vafs):.4f}'})
        rows.append({'Category': '  Median', 'Value': f'{sorted(vafs)[len(vafs)//2]:.4f}'})
        rows.append({'Category': '  Mean', 'Value': f'{sum(vafs)/len(vafs):.4f}'})

    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Generate final clinical report integrating standard format + Tier classification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 基本用法
  python make_final_report.py \\
      -i sample1.snv.standard_report.txt \\
      -o sample1.snv.final_report.xlsx \\
      --tier-dir ./variants/classification/

  # 无 Tier 分级（仅标准格式输出到 Excel）
  python make_final_report.py \\
      -i sample1.snv.standard_report.txt \\
      -o sample1.snv.final_report.xlsx
        """
    )

    parser.add_argument("--input", "-i", required=True,
                        help="Standard format TSV from filter_by_genelist.py")
    parser.add_argument("--output", "-o", required=True,
                        help="Output Excel report file")
    parser.add_argument("--tier-dir", "-t", default=None,
                        help="Directory containing Tier classification output files")
    parser.add_argument("--sample-name", "-s", default=None,
                        help="Sample name (auto-detected from input if not specified)")

    args = parser.parse_args()

    # 检查输入
    if not Path(args.input).exists():
        print(f"[ERROR] Input file not found: {args.input}")
        sys.exit(1)

    # 读取标准格式数据
    print(f"[INFO] Loading standard report: {args.input}")
    standard_df = pd.read_csv(args.input, sep='\t', low_memory=False)
    print(f"[INFO]   Loaded {len(standard_df)} variants in standard format")

    # 确定样本名
    sample_name = args.sample_name
    if sample_name is None:
        sample_name = standard_df['样本编号'].iloc[0] if len(standard_df) > 0 else 'unknown'

    # 加载 Tier 分级
    tier_df = None
    if args.tier_dir and Path(args.tier_dir).exists():
        print(f"[INFO] Searching for Tier classification in: {args.tier_dir}")

        # 方法 1: 加载 Excel
        tier_df = load_classification_results(args.tier_dir, sample_name)

        # 方法 2: 加载 TSV
        if tier_df is None:
            tier_df = load_tier_tsv_files(args.tier_dir)

    if tier_df is None:
        print("[WARN] No Tier classification loaded. Output will only contain standard format.")
    else:
        print(f"[INFO] Loaded {len(tier_df)} Tier-classified variants")

    # 合并
    final_df = merge_with_tiers(standard_df, tier_df)

    # 输出
    write_final_excel(final_df, args.output, sample_name)

    # 统计
    print(f"\n{'='*60}")
    print(f"Final Report Summary - {sample_name}")
    print(f"{'='*60}")
    print(f"  Total variants: {len(final_df)}")

    for tier in ['Tier1_actionable', 'Tier2_potential', 'Tier3_unknown', 'Tier4_benign', 'Unclassified']:
        count = (final_df['Tier'] == tier).sum()
        if count > 0:
            label = TIER_LABELS.get(tier, tier)
            print(f"  {label}: {count}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
