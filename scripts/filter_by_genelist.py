#!/usr/bin/env python3
"""
filter_by_genelist.py - 根据基因列表筛选 ANNOVAR 结果并转换为标准上报格式

标准上报格式（12 列，TSV）：
    样本编号 | Chr | Start | End | Ref | Alt | Gene |
    Type | Transcript | cHGVS | pHGVS | VAF

报告规范：
  1. 3' rule：所有变异位点位置按靠近基因转录方向 3' 端表述
  2. Gene：遵循 HGVS 基因命名规则
  3. Transcript / cHGVS / pHGVS：参考 ClinVar 写法
  4. Type：SNV / Insertion / Deletion / Complex
  5. VAF：保留两位小数，无需填写 %
"""

import pandas as pd
import argparse
import re
from pathlib import Path
from typing import Optional, Tuple, Dict, List


# ============================================================
# 3' rule 位置计算
# ============================================================

def apply_3prime_rule(original_start: int, ref: str, alt: str) -> int:
    """
    根据 3' rule 计算标准化 Start 位置。

    规则：
    - SNV（ref_len==1, alt_len==1）：位置不变
    - Insertion（alt_len > ref_len）：位置 = original_start + ref_len
      （插入位的 3' 侧碱基）
    - Deletion（ref_len > alt_len）：位置 = original_start + ref_len - 1
      （被删除片段中最靠近 3' 端的碱基）
    - Complex（ref_len == alt_len 且 >1 或多碱基替换）：
      位置 = original_start + max(ref_len, alt_len) - 1

    Returns:
        3' rule 调整后的 Start 坐标
    """
    ref_len = len(ref) if ref and ref != '.' else 0
    alt_len = len(alt) if alt and alt != '.' else 0

    if ref_len == 0 or alt_len == 0:
        return original_start  # 异常情况，不调整

    if ref_len == 1 and alt_len == 1:
        # SNV：位置不变
        return original_start
    elif alt_len > ref_len:
        # Insertion：3' 侧 = 插入点下一位
        return original_start + ref_len
    elif ref_len > alt_len:
        # Deletion：最 3' 被删碱基
        return original_start + ref_len - 1
    else:
        # Complex：3' 侧位置
        return original_start + max(ref_len, alt_len) - 1


# ============================================================
# HGVS / ClinVar 风格转录本与命名提取
# ============================================================

def parse_aachange(aa_str: str) -> Dict[str, str]:
    """
    从 ANNOVAR AAChange.refGene 字段提取转录本、cHGVS、pHGVS。

    AAChange.refGene 格式示例：
        NM_000546.6(TP53):exon5:c.524G>A:p.R175H,
        NM_001126112.2(TP53):exon5:c.524G>A:p.R175H

    策略：
    1. 多个转录本按逗号分隔
    2. 每个转录本解析为 {transcript, exon, cHGVS, pHGVS}
    3. 优先选择第一个转录本（ANNOVAR 默认优先输出 canonical）
    4. 转录本名称遵循 ClinVar 偏好（NM_ 转录本优先）
    """
    result = {
        'transcript': '.',
        'cHGVS': '.',
        'pHGVS': '.',
        'all_transcripts': []
    }

    if not aa_str or pd.isna(aa_str) or str(aa_str).strip() == '.':
        return result

    aa_str = str(aa_str).strip()

    # 按逗号分隔多个转录本
    entries = [e.strip() for e in aa_str.split(',') if e.strip()]

    parsed_entries = []
    for entry in entries:
        # 格式：Transcript(Gene):exonN:c.XXX:p.XXX
        # 可能有不完整的情况，用正则强健解析
        parts = entry.split(':')

        transcript = '.'
        exon = '.'
        chgvs = '.'
        phgvs = '.'

        for part in parts:
            part = part.strip()

            # 转录本：NM_xxxx 或 NR_xxxx 或 XM_xxxx
            if re.match(r'^[NXY]M_\d+', part):
                # 提取纯 transcript（去掉括号中的基因名）
                m = re.match(r'^([NXY]M_\d+\.?\d*)', part)
                if m:
                    transcript = m.group(1)

            # exon
            elif part.lower().startswith('exon'):
                exon = part

            # cHGVS
            elif part.startswith('c.'):
                chgvs = part

            # pHGVS
            elif part.startswith('p.'):
                phgvs = part

        if transcript != '.' or chgvs != '.' or phgvs != '.':
            parsed_entries.append({
                'transcript': transcript,
                'exon': exon,
                'cHGVS': chgvs,
                'pHGVS': phgvs
            })

    result['all_transcripts'] = parsed_entries

    if parsed_entries:
        # 优先 NM_ 转录本（ClinVar 标准）
        nm_entries = [e for e in parsed_entries if e['transcript'].startswith('NM_')]
        selected = nm_entries[0] if nm_entries else parsed_entries[0]

        result['transcript'] = selected['transcript']
        result['cHGVS'] = selected['cHGVS']
        result['pHGVS'] = selected['pHGVS']

    return result


def extract_hgvs_gene_name(gene_str: str) -> str:
    """
    提取 HGVS 标准基因名称。

    规则：
    - ANNOVAR Gene.refGene 中的第一个基因名即为 HGVS 官方符号
    - 多基因注释用逗号分隔，取第一个
    - 特殊处理：如 "MIRxxx" 等非编码基因保留原名
    """
    if not gene_str or pd.isna(gene_str) or str(gene_str).strip() == '.':
        return '.'

    genes = [g.strip() for g in str(gene_str).split(',') if g.strip()]
    if not genes:
        return '.'

    gene = genes[0]

    # 标准化：去除常见后缀
    # 某些数据库会在基因名后加编号，如 TP53;TP53 取第一个
    if ';' in gene:
        gene = gene.split(';')[0].strip()

    return gene


# ============================================================
# 变异类型分类
# ============================================================

def classify_variant_type(
    ref: str,
    alt: str,
    func: Optional[str] = None,
    exonic_func: Optional[str] = None
) -> str:
    """
    根据 Ref/Alt 序列和功能注释判定变异类型。

    返回：SNV / Insertion / Deletion / Complex
    """
    ref = str(ref) if ref and ref != '.' else ''
    alt = str(alt) if alt and alt != '.' else ''

    if not ref or not alt:
        # 尝试从功能注释推断
        if func or exonic_func:
            combined = str(func or '') + ' ' + str(exonic_func or '')
            combined = combined.lower()
            if 'insertion' in combined:
                return 'Insertion'
            elif 'deletion' in combined:
                return 'Deletion'
            elif 'nonsynonymous' in combined or 'stop' in combined:
                return 'SNV'
        return 'SNV'  # 默认

    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1 and alt_len == 1:
        return 'SNV'
    elif alt_len > ref_len:
        # 检查是否纯插入（ref 是 alt 的前缀）
        if alt.startswith(ref):
            return 'Insertion'
        else:
            return 'Complex'
    elif ref_len > alt_len:
        return 'Deletion'
    else:
        # ref_len == alt_len > 1: MNV 或 Complex
        return 'Complex'


# ============================================================
# VAF 格式化
# ============================================================

def format_vaf(vaf_value) -> str:
    """
    格式化 VAF：保留两位小数，无需 % 号。

    输入：
      - 0.0~1.0 范围的比例值
      - >1 的百分数值（自动/100）
      - '.' 或缺失值返回 '.'

    Returns:
        str: 如 "0.35", "0.05", "."
    """
    if vaf_value is None or pd.isna(vaf_value) or vaf_value == '.' or vaf_value == '':
        return '.'

    try:
        v = float(vaf_value)
        if v > 1.0:
            v = v / 100.0
        if v < 0:
            return '.'
        return f"{v:.2f}"
    except (ValueError, TypeError):
        return '.'


# ============================================================
# 主处理函数
# ============================================================

def filter_and_format(
    multianno_file: str,
    gene_list_file: str,
    output_file: str,
    sample_name: Optional[str] = None,
    vaf_col_override: Optional[str] = None
):
    """
    主流程：读入 ANNOVAR 结果 → 基因列表筛选 → 标准上报格式输出
    """
    # ===== 读取基因列表 =====
    if not Path(gene_list_file).exists():
        raise FileNotFoundError(f"Gene list file not found: {gene_list_file}")

    with open(gene_list_file) as f:
        target_genes = set(line.strip() for line in f if line.strip() and not line.startswith('#'))

    print(f"[INFO] Target genes loaded: {len(target_genes)}")
    if len(target_genes) <= 10:
        print(f"       Genes: {', '.join(sorted(target_genes))}")

    # ===== 读取 ANNOVAR 结果 =====
    if not Path(multianno_file).exists():
        raise FileNotFoundError(f"ANNOVAR file not found: {multianno_file}")

    df = pd.read_csv(multianno_file, sep='\t', low_memory=False)
    print(f"[INFO] Total variants in ANNOVAR output: {len(df)}")

    # ===== 列名自动检测 =====
    col_map = _detect_columns(df)
    print(f"[INFO] Detected columns: {col_map}")

    # ===== 基因筛选 =====
    gene_col = col_map['gene']
    if gene_col not in df.columns:
        # 回退查找
        for candidate in ['Gene.refGene', 'Gene_refGene', 'Gene.knownGene', 'GENE', 'Gene']:
            if candidate in df.columns:
                gene_col = candidate
                col_map['gene'] = candidate
                break
        else:
            raise KeyError(f"No gene column found. Available: {list(df.columns)[:20]}")

    def match_gene(gene_str):
        if pd.isna(gene_str) or str(gene_str).strip() in ('.', ''):
            return False
        genes = set(g.strip() for g in str(gene_str).split(',') if g.strip())
        return bool(genes & target_genes)

    mask = df[gene_col].apply(match_gene)
    filtered = df[mask].copy()
    print(f"[INFO] Variants matched to gene list: {len(filtered)}")

    # ===== 确定样本名 =====
    if sample_name is None:
        sample_name = Path(multianno_file).stem
        sample_name = re.sub(r'\.snv\.hg19_multianno$', '', sample_name)
        sample_name = re.sub(r'\.hg19_multianno$', '', sample_name)

    # ===== 构建输出行 =====
    output_columns = [
        "样本编号", "Chr", "Start", "End", "Ref", "Alt",
        "Gene", "Type", "Transcript", "cHGVS", "pHGVS", "VAF"
    ]

    result_rows = []
    skipped_count = 0

    for _, row in filtered.iterrows():
        try:
            # --- Chr ---
            chrom_raw = str(row[col_map['chr']]) if col_map['chr'] in row.index else '.'
            chrom = chrom_raw.replace('chr', '').replace('Chr', '').replace('CHR', '')
            # 确保纯数字（去除非数字字符，但保留 XYMT）
            if chrom.upper() not in ('X', 'Y', 'M', 'MT') and not chrom.isdigit():
                chrom = re.sub(r'[^0-9XYMT]', '', chrom.upper())

            # --- Start (原始) ---
            try:
                original_start = int(row[col_map['start']])
            except (ValueError, TypeError, KeyError):
                original_start = 0

            # --- End (原始) ---
            try:
                original_end = int(row[col_map['end']]) if col_map['end'] in row.index else original_start
            except (ValueError, TypeError):
                original_end = original_start

            # --- Ref ---
            ref = str(row[col_map['ref']]) if col_map['ref'] in row.index else '.'
            if ref in ('.', 'nan', ''):
                ref = '.'

            # --- Alt ---
            alt = str(row[col_map['alt']]) if col_map['alt'] in row.index else '.'
            if alt in ('.', 'nan', ''):
                alt = '.'

            # --- 3' rule 调整 Start/End ---
            if ref != '.' and alt != '.':
                adjusted_start = apply_3prime_rule(original_start, ref, alt)

                ref_len = len(ref)
                alt_len = len(alt)

                if alt_len > ref_len:
                    # Insertion：End = 3' Start + 插入长度（实际基因组范围）
                    insert_len = alt_len - ref_len
                    adjusted_end = original_start + ref_len + insert_len - 1
                elif ref_len > alt_len:
                    # Deletion：End 同 3' Start
                    adjusted_end = adjusted_start
                else:
                    # SNV 或 Complex：End = Start
                    adjusted_end = adjusted_start
            else:
                adjusted_start = original_start
                adjusted_end = original_end

            # --- Gene (HGVS 命名) ---
            gene = extract_hgvs_gene_name(row[gene_col])

            # --- Type ---
            func_val = None
            exonic_func_val = None
            if col_map.get('func') and col_map['func'] in row.index:
                func_val = row[col_map['func']]
            if col_map.get('exonic_func') and col_map['exonic_func'] in row.index:
                exonic_func_val = row[col_map['exonic_func']]

            vtype = classify_variant_type(ref, alt, func_val, exonic_func_val)

            # --- Transcript / cHGVS / pHGVS (ClinVar 风格) ---
            aachange_raw = None
            if col_map.get('aachange') and col_map['aachange'] in row.index:
                aachange_raw = row[col_map['aachange']]

            aachange_parsed = parse_aachange(aachange_raw)
            transcript = aachange_parsed['transcript']
            chgvs = aachange_parsed['cHGVS']
            phgvs = aachange_parsed['pHGVS']

            # 如果 AAChange 解析失败，尝试从 VCF INFO 列恢复
            if transcript == '.':
                # 尝试从其他可能包含转录本信息的列获取
                for candidate_col in ['Transcript', 'TRANSCRIPT', 'NM']:
                    if candidate_col in row.index:
                        val = str(row[candidate_col])
                        if val not in ('.', 'nan', ''):
                            transcript = val
                            break

            # --- VAF ---
            vaf_raw = '.'
            if vaf_col_override and vaf_col_override in row.index:
                vaf_raw = row[vaf_col_override]
            elif col_map.get('vaf') and col_map['vaf'] in row.index:
                vaf_raw = row[col_map['vaf']]

            vaf_str = format_vaf(vaf_raw)

            result_rows.append([
                sample_name,
                chrom,
                adjusted_start,
                adjusted_end,
                ref,
                alt,
                gene,
                vtype,
                transcript,
                chgvs,
                phgvs,
                vaf_str
            ])
        except Exception as e:
            skipped_count += 1
            if skipped_count <= 3:
                print(f"[WARN] Skipped variant due to error: {e}")
            continue

    # ===== 输出 =====
    result_df = pd.DataFrame(result_rows, columns=output_columns)

    # 去重（同样本、同基因、同 cHGVS 只保留第一个）
    result_df = result_df.drop_duplicates(
        subset=["样本编号", "Gene", "cHGVS"],
        keep='first'
    )

    # 按基因、染色体、位置排序
    result_df['_chr_order'] = result_df['Chr'].apply(_chr_sort_key)
    result_df = result_df.sort_values(
        ['样本编号', '_chr_order', 'Start', 'Gene']
    ).drop(columns=['_chr_order'])

    result_df.to_csv(output_file, sep='\t', index=False, na_rep='.')

    print(f"\n[INFO] Standard report written: {output_file}")
    print(f"[INFO] Total variants: {len(result_df)}")
    if skipped_count > 0:
        print(f"[WARN] Skipped variants: {skipped_count}")

    # ===== 统计 =====
    if len(result_df) > 0:
        print(f"\n{'='*60}")
        print(f"Gene distribution:")
        for gene, count in result_df['Gene'].value_counts().head(20).items():
            print(f"  {gene}: {count}")

        print(f"\nType distribution:")
        for vtype, count in result_df['Type'].value_counts().items():
            print(f"  {vtype}: {count}")

        print(f"\nVAF range: {_vaf_range_summary(result_df)}")
        print(f"{'='*60}")

    return result_df


# ============================================================
# 辅助函数
# ============================================================

def _detect_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    """
    自动检测 ANNOVAR 输出中的关键列名。
    返回列名映射字典。
    """
    cols = {c: None for c in [
        'chr', 'start', 'end', 'ref', 'alt',
        'gene', 'func', 'exonic_func', 'aachange', 'vaf'
    ]}

    # 染色体
    for c in ['Chr', '#Chr', 'CHROM', 'CHR']:
        if c in df.columns:
            cols['chr'] = c
            break

    # 起始位置
    for c in ['Start', 'POS', 'POSITION']:
        if c in df.columns:
            cols['start'] = c
            break

    # 终止位置
    for c in ['End', 'END']:
        if c in df.columns:
            cols['end'] = c
            break

    # 参考碱基
    for c in ['Ref', 'REF', 'REFERENCE']:
        if c in df.columns:
            cols['ref'] = c
            break

    # 变异碱基
    for c in ['Alt', 'ALT', 'ALLELE']:
        if c in df.columns:
            cols['alt'] = c
            break

    # 基因
    for c in ['Gene.refGene', 'Gene_refGene', 'Gene.knownGene', 'GENE', 'Gene']:
        if c in df.columns:
            cols['gene'] = c
            break

    # 功能区
    for c in ['Func.refGene', 'Func_refGene', 'Func.knownGene']:
        if c in df.columns:
            cols['func'] = c
            break

    # 外显子功能
    for c in ['ExonicFunc.refGene', 'ExonicFunc_refGene', 'ExonicFunc.knownGene']:
        if c in df.columns:
            cols['exonic_func'] = c
            break

    # AA 变化
    for c in ['AAChange.refGene', 'AAChange_refGene', 'AAChange.knownGene']:
        if c in df.columns:
            cols['aachange'] = c
            break

    # VAF
    vaf_candidates = [
        'AF', 'VAF', 'ALT_FREQ', 'gnomAD_AF',
        'AF_TUMOR', 'TUMOR_AF', 'tumor_AF'
    ]
    for c in vaf_candidates:
        if c in df.columns:
            cols['vaf'] = c
            break
    if cols['vaf'] is None:
        for c in df.columns:
            if 'af' in c.lower() and ('tumor' in c.lower() or 'vaf' in c.lower()):
                cols['vaf'] = c
                break
    if cols['vaf'] is None:
        for c in df.columns:
            if 'vaf' in c.lower() or 'freq' in c.lower():
                cols['vaf'] = c
                break

    return cols


def _chr_sort_key(chrom: str) -> int:
    """染色体排序键：1-22, X=23, Y=24, M/MT=25, 其他=99"""
    chrom = str(chrom).upper().strip()
    if chrom.isdigit():
        return int(chrom)
    elif chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    elif chrom in ('M', 'MT'):
        return 25
    else:
        return 99


def _vaf_range_summary(df: pd.DataFrame) -> str:
    """VAF 列统计摘要"""
    vafs = []
    for v in df['VAF']:
        try:
            vafs.append(float(v))
        except (ValueError, TypeError):
            pass
    if not vafs:
        return "N/A"
    return f"min={min(vafs):.4f}, max={max(vafs):.4f}, median={sorted(vafs)[len(vafs)//2]:.4f}"


# ============================================================
# 命令行入口
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Filter ANNOVAR results by gene list and format to standard clinical report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 基本用法
  python filter_by_genelist.py \\
      -m sample1.snv.hg19_multianno.txt \\
      -g target.genes.txt \\
      -o sample1.snv.standard_report.txt

  # 指定样本名和 VAF 列
  python filter_by_genelist.py \\
      -m sample1.snv.hg19_multianno.txt \\
      -g target.genes.txt \\
      -o sample1.snv.standard_report.txt \\
      --sample-name T202611 \\
      --vaf-col AF_TUMOR
        """
    )
    parser.add_argument("--multianno", "-m", required=True,
                        help="ANNOVAR hg19_multianno.txt output file")
    parser.add_argument("--gene-list", "-g", required=True,
                        help="Target gene list (one HGVS gene symbol per line)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output standard report file (TSV format)")
    parser.add_argument("--sample-name", "-s", default=None,
                        help="Sample identifier (auto-detected from filename if not specified)")
    parser.add_argument("--vaf-col", default=None,
                        help="VAF column name override (auto-detected if not specified)")

    args = parser.parse_args()

    filter_and_format(
        multianno_file=args.multianno,
        gene_list_file=args.gene_list,
        output_file=args.output,
        sample_name=args.sample_name,
        vaf_col_override=args.vaf_col
    )


if __name__ == "__main__":
    main()
