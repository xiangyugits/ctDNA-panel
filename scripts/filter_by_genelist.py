#!/usr/bin/env python3
"""
filter_by_genelist.py - 根据基因列表筛选 ANNOVAR 结果并格式化输出
"""

import pandas as pd
import argparse
from pathlib import Path


def filter_and_format(
    multianno_file: str,
    gene_list_file: str,
    output_file: str
):
    """
    根据基因列表筛选 ANNOVAR 注释结果并格式化输出

    输出列：
    - 样本编号
    - Chr: 纯数字
    - Start: 3' rule
    - End
    - Ref
    - Alt
    - Gene: HGVS 命名
    - Type: SNV/Insertion/Deletion/Complex
    - Transcript: ClinVar 写法
    - cHGVS: ClinVar 写法
    - pHGVS: ClinVar 写法
    - VAF (%): 小数点后两位
    """

    # ===== 读取基因列表 =====
    if not Path(gene_list_file).exists():
        raise FileNotFoundError(f"Gene list not found: {gene_list_file}")

    with open(gene_list_file) as f:
        target_genes = set(line.strip() for line in f if line.strip())

    print(f"Target genes: {len(target_genes)}")

    # ===== 读取 ANNOVAR 结果 =====
    if not Path(multianno_file).exists():
        raise FileNotFoundError(f"ANNOVAR file not found: {multianno_file}")

    df = pd.read_csv(multianno_file, sep='\t', low_memory=False)
    print(f"Total variants in ANNOVAR: {len(df)}")

    # ===== 查找基因列 =====
    gene_col = None
    for col in ['Gene.refGene', 'Gene_refGene', 'Gene.knownGene', 'GENE']:
        if col in df.columns:
            gene_col = col
            break

    if gene_col is None:
        print("ERROR: No gene column found in ANNOVAR output")
        print(f"Available columns: {df.columns.tolist()}")
        return

    # ===== 筛选目标基因 =====
    def match_gene(gene_str):
        if pd.isna(gene_str) or gene_str == '.':
            return False
        genes = set(str(gene_str).split(','))
        genes = set(g.strip() for g in genes)
        return bool(genes & target_genes)

    mask = df[gene_col].apply(match_gene)
    filtered = df[mask].copy()
    print(f"Matched variants: {len(filtered)}")

    if filtered.empty:
        print("Warning: No variants matched the gene list!")
        # 创建空输出
        with open(output_file, 'w') as f:
            f.write("样本编号\tChr\tStart\tEnd\tRef\tAlt\tGene\tType\tTranscript\tcHGVS\tpHGVS\tVAF (%)\n")
        return

    # ===== 查找各列 =====
    chr_col = None
    for col in ['Chr', '#Chr', 'CHROM']:
        if col in filtered.columns:
            chr_col = col
            break

    start_col = None
    for col in ['Start', 'POS']:
        if col in filtered.columns:
            start_col = col
            break

    end_col = None
    for col in ['End', 'END']:
        if col in filtered.columns:
            end_col = col
            break

    ref_col = None
    for col in ['Ref', 'REF']:
        if col in filtered.columns:
            ref_col = col
            break

    alt_col = None
    for col in ['Alt', 'ALT']:
        if col in filtered.columns:
            alt_col = col
            break

    func_col = None
    for col in ['ExonicFunc.refGene', 'ExonicFunc_refGene', 'Func.refGene', 'Func_refGene']:
        if col in filtered.columns:
            func_col = col
            break

    aachange_col = None
    for col in ['AAChange.refGene', 'AAChange_refGene']:
        if col in filtered.columns:
            aachange_col = col
            break

    clinvar_col = None
    for col in ['clinvar_20231231_clnsig', 'CLINSIG', 'clinvar']:
        if col in filtered.columns:
            clinvar_col = col
            break

    # ===== 查找 VAF 列 =====
    vaf_col = None
    for col in ['AF', 'VAF', 'ALT_FREQ', 'gnomAD_AF']:
        if col in filtered.columns:
            vaf_col = col
            break
    if vaf_col is None:
        for col in filtered.columns:
            if 'af' in col.lower() or 'vaf' in col.lower() or 'freq' in col.lower():
                vaf_col = col
                break

    # ===== 提取样本编号 =====
    # 从文件名或路径推断样本名
    sample_name = Path(multianno_file).stem.replace('.snv.hg19_multianno', '').replace('.hg19_multianno', '')

    # ===== 构建输出 =====
    result_rows = []

    for _, row in filtered.iterrows():
        # --- Chr: 纯数字 ---
        chrom = str(row[chr_col]) if chr_col else '.'
        chrom = chrom.replace('chr', '').replace('Chr', '').replace('CHR', '')

        # --- Start (3' rule) ---
        start = row[start_col] if start_col else '.'

        # --- End ---
        end = row[end_col] if end_col else '.'

        # --- Ref ---
        ref = str(row[ref_col]) if ref_col else '.'

        # --- Alt ---
        alt = str(row[alt_col]) if alt_col else '.'

        # --- Gene (HGVS 命名, 取第一个基因) ---
        gene_str = str(row[gene_col]) if pd.notna(row[gene_col]) else '.'
        gene = gene_str.split(',')[0].strip()

        # --- Transcript ---
        transcript = '.'
        if aachange_col and pd.notna(row[aachange_col]) and row[aachange_col] != '.':
            parts = str(row[aachange_col]).split(':')
            if len(parts) >= 1:
                transcript = parts[0].split(',')[0].strip()

        # --- cHGVS ---
        chgvs = '.'
        if aachange_col and pd.notna(row[aachange_col]) and row[aachange_col] != '.':
            aa_str = str(row[aachange_col])
            # 尝试匹配 c. 开头的部分
            for part in aa_str.split(':'):
                part = part.strip()
                if part.startswith('c.'):
                    chgvs = part
                    break
            # 如果没找到 c.，但 AAChange 包含 exon 信息
            if chgvs == '.':
                # 从多个转录本中找
                for segment in aa_str.split(','):
                    for part in segment.split(':'):
                        part = part.strip()
                        if part.startswith('c.'):
                            chgvs = part
                            break

        # --- pHGVS ---
        phgvs = '.'
        if aachange_col and pd.notna(row[aachange_col]) and row[aachange_col] != '.':
            aa_str = str(row[aachange_col])
            for part in aa_str.split(':'):
                part = part.strip()
                if part.startswith('p.'):
                    phgvs = part
                    break
            if phgvs == '.':
                for segment in aa_str.split(','):
                    for part in segment.split(':'):
                        part = part.strip()
                        if part.startswith('p.'):
                            phgvs = part
                            break

        # --- Type: SNV/Insertion/Deletion/Complex ---
        vtype = 'SNV'  # 默认
        if func_col and pd.notna(row[func_col]) and row[func_col] != '.':
            func_val = str(row[func_col]).lower()
            if 'frameshift' in func_val:
                if 'deletion' in func_val:
                    vtype = 'Deletion'
                elif 'insertion' in func_val:
                    vtype = 'Insertion'
                else:
                    vtype = 'Complex'
            elif 'nonframeshift' in func_val:
                if 'deletion' in func_val:
                    vtype = 'Deletion'
                elif 'insertion' in func_val:
                    vtype = 'Insertion'
                else:
                    vtype = 'Complex'
            elif 'stopgain' in func_val or 'stoploss' in func_val:
                vtype = 'SNV'
            elif 'nonsynonymous' in func_val:
                vtype = 'SNV'
            elif 'synonymous' in func_val:
                vtype = 'SNV'
            elif 'splicing' in func_val or 'splice' in func_val:
                vtype = 'SNV'
        elif ref_col and alt_col:
            ref_len = len(str(row[ref_col]))
            alt_len = len(str(row[alt_col]))
            if ref_len == 1 and alt_len == 1:
                vtype = 'SNV'
            elif alt_len > ref_len:
                vtype = 'Insertion'
            elif alt_len < ref_len:
                vtype = 'Deletion'
            else:
                vtype = 'Complex'

        # --- VAF (%) - 两位小数 ---
        vaf_str = '.'
        if vaf_col and pd.notna(row[vaf_col]) and row[vaf_col] != '.':
            try:
                vaf_val = float(row[vaf_col])
                if vaf_val > 1:
                    vaf_val = vaf_val / 100
                vaf_str = f"{vaf_val:.2f}"
            except (ValueError, TypeError):
                vaf_str = '.'

        result_rows.append([
            sample_name,
            chrom,
            start,
            end,
            ref,
            alt,
            gene,
            vtype,
            transcript,
            chgvs,
            phgvs,
            vaf_str
        ])

    # ===== 输出 =====
    columns = ["样本编号", "Chr", "Start", "End", "Ref", "Alt", "Gene",
               "Type", "Transcript", "cHGVS", "pHGVS", "VAF (%)"]
    result_df = pd.DataFrame(result_rows, columns=columns)
    result_df.to_csv(output_file, sep='\t', index=False)

    print(f"\nOutput: {output_file}")
    print(f"Total variants: {len(result_df)}")

    # ===== 统计 =====
    print(f"\nGene distribution:")
    for gene, count in result_df['Gene'].value_counts().items():
        print(f"  {gene}: {count}")

    print(f"\nType distribution:")
    for vtype, count in result_df['Type'].value_counts().items():
        print(f"  {vtype}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter ANNOVAR results by gene list and format output"
    )
    parser.add_argument("--multianno", "-m", required=True,
                       help="ANNOVAR hg19_multianno.txt file")
    parser.add_argument("--gene-list", "-g", required=True,
                       help="Target gene list file (one gene per line)")
    parser.add_argument("--output", "-o", required=True,
                       help="Output formatted file")

    args = parser.parse_args()

    filter_and_format(
        multianno_file=args.multianno,
        gene_list_file=args.gene_list,
        output_file=args.output
    )


if __name__ == "__main__":
    main()