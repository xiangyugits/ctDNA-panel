#!/usr/bin/env python3
"""
filter_cnv.py - 过滤 CNVkit 检测的拷贝数变异结果
用于 ctDNA Panel 分析，基于 log2 比值、探针数和置信区间过滤
"""

import argparse
import sys
from pathlib import Path

def filter_cnv(
    input_file: str,
    output_file: str,
    min_log2_abs: float = 0.3,
    min_probes: int = 5,
    max_ci: float = 1.0,
    min_depth: float = 50,
    gene_list: str = None,
    gene_column: str = "gene"
):
    """
    过滤 CNVkit 输出的 .cns 文件
    
    Parameters
    ----------
    input_file : str
        CNVkit 输出的 .cns 文件路径
    output_file : str
        过滤后的输出文件路径
    min_log2_abs : float
        log2 比值绝对值的最小阈值 (默认 0.3)
        0.3 约等于 1.23 倍拷贝数变化
    min_probes : int
        最少连续探针数 (默认 5)
    max_ci : float
        最大置信区间宽度 (默认 1.0)
        置信区间过宽说明拷贝数估计不可靠
    min_depth : float
        最小覆盖深度 (默认 50)
    gene_list : str
        可选，关注的基因列表文件（一行一个基因）
        如果提供，只保留包含这些基因的 CNV
    gene_column : str
        .cns 文件中基因列的名称 (默认 "gene")
    """
    
    # 读取基因列表
    target_genes = set()
    if gene_list:
        if Path(gene_list).exists():
            with open(gene_list) as f:
                target_genes = {line.strip() for line in f if line.strip()}
            print(f"Loaded {len(target_genes)} target genes from {gene_list}")
        else:
            print(f"Warning: gene list file not found: {gene_list}")
    
    # 统计
    stats = {
        "total": 0,
        "passed_log2": 0,
        "passed_probes": 0,
        "passed_ci": 0,
        "passed_depth": 0,
        "passed_gene": 0,
        "final": 0,
        "amplification": 0,
        "deletion": 0
    }
    
    with open(input_file) as f_in, open(output_file, 'w') as f_out:
        # 读取 header
        header = f_in.readline().strip()
        f_out.write(header + "\n")
        
        # 解析列名
        cols = header.split('\t')
        col_idx = {name: i for i, name in enumerate(cols)}
        
        # 检查必需的列
        required_cols = ["log2", "probes", "ci"]
        for col in required_cols:
            if col not in col_idx:
                # 尝试常见别名
                if col == "log2" and "log2" not in col_idx:
                    if "log2ratio" in col_idx:
                        col_idx["log2"] = col_idx["log2ratio"]
                    else:
                        raise ValueError(f"Cannot find 'log2' or 'log2ratio' column in {input_file}")
                elif col == "probes":
                    if "probes" not in col_idx and "n_probes" in col_idx:
                        col_idx["probes"] = col_idx["n_probes"]
                    else:
                        raise ValueError(f"Cannot find 'probes' column in {input_file}")
                elif col == "ci":
                    if "ci" not in col_idx:
                        print("Warning: 'ci' column not found, skipping CI filter")
                        col_idx["ci"] = None
        
        # 逐行过滤
        for line in f_in:
            if not line.strip():
                continue
            
            stats["total"] += 1
            fields = line.strip().split('\t')
            
            # 解析各字段
            try:
                log2_val = float(fields[col_idx["log2"]])
                probes = int(fields[col_idx["probes"]])
                ci_val = None
                if col_idx.get("ci") is not None:
                    try:
                        ci_val = float(fields[col_idx["ci"]])
                    except (ValueError, IndexError):
                        ci_val = None
            except (ValueError, IndexError) as e:
                print(f"Warning: skipping line {stats['total']}: {e}")
                continue
            
            # 1. log2 比值过滤
            if abs(log2_val) < min_log2_abs:
                continue
            stats["passed_log2"] += 1
            
            # 2. 探针数过滤
            if probes < min_probes:
                continue
            stats["passed_probes"] += 1
            
            # 3. 置信区间过滤
            if ci_val is not None and ci_val > max_ci:
                continue
            stats["passed_ci"] += 1
            
            # 4. 深度过滤（如果列存在）
            if "depth" in col_idx:
                try:
                    depth = float(fields[col_idx["depth"]])
                    if depth < min_depth:
                        continue
                except (ValueError, IndexError):
                    pass
            stats["passed_depth"] += 1
            
            # 5. 基因过滤
            if target_genes:
                gene_str = fields[col_idx.get(gene_column, -1)] if gene_column in col_idx else ""
                if not gene_str:
                    continue
                genes = set(gene_str.split(','))
                if not genes.intersection(target_genes):
                    continue
            stats["passed_gene"] += 1
            
            # 通过所有过滤
            f_out.write(line)
            stats["final"] += 1
            
            # 统计扩增/缺失
            if log2_val > 0:
                stats["amplification"] += 1
            else:
                stats["deletion"] += 1
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Filter CNVkit .cns output for ctDNA panel analysis"
    )
    
    parser.add_argument("--input", "-i", required=True,
                       help="CNVkit .cns input file")
    parser.add_argument("--output", "-o", required=True,
                       help="Filtered output file")
    
    # 过滤参数
    parser.add_argument("--min-log2-abs", type=float, default=0.3,
                       help="Minimum absolute log2 ratio (default: 0.3)")
    parser.add_argument("--min-probes", type=int, default=5,
                       help="Minimum number of probes in segment (default: 5)")
    parser.add_argument("--max-ci", type=float, default=1.0,
                       help="Maximum confidence interval width (default: 1.0)")
    parser.add_argument("--min-depth", type=float, default=50,
                       help="Minimum coverage depth (default: 50)")
    
    # 基因过滤
    parser.add_argument("--gene-list",
                       help="File with target genes (one per line)")
    parser.add_argument("--gene-column", default="gene",
                       help="Gene column name in .cns file (default: gene)")
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # 执行过滤
    print("=" * 60)
    print("CNV Filtering Parameters")
    print("=" * 60)
    print(f"  Input: {args.input}")
    print(f"  Output: {args.output}")
    print(f"  Min |log2|: {args.min_log2_abs}")
    print(f"  Min probes: {args.min_probes}")
    print(f"  Max CI: {args.max_ci}")
    print(f"  Min depth: {args.min_depth}")
    if args.gene_list:
        print(f"  Gene list: {args.gene_list}")
    print()
    
    stats = filter_cnv(
        input_file=args.input,
        output_file=args.output,
        min_log2_abs=args.min_log2_abs,
        min_probes=args.min_probes,
        max_ci=args.max_ci,
        min_depth=args.min_depth,
        gene_list=args.gene_list,
        gene_column=args.gene_column
    )
    
    # 打印统计
    print("=" * 60)
    print("CNV Filtering Summary")
    print("=" * 60)
    print(f"  Total segments:        {stats['total']:>8}")
    print(f"  Passed log2 filter:    {stats['passed_log2']:>8}")
    print(f"  Passed probes filter:  {stats['passed_probes']:>8}")
    print(f"  Passed CI filter:      {stats['passed_ci']:>8}")
    print(f"  Passed depth filter:   {stats['passed_depth']:>8}")
    if args.gene_list:
        print(f"  Passed gene filter:    {stats['passed_gene']:>8}")
    print(f"  ─────────────────────────────")
    print(f"  Final passed:          {stats['final']:>8}")
    print(f"    Amplifications:      {stats['amplification']:>8}")
    print(f"    Deletions:           {stats['deletion']:>8}")
    print(f"  Filter rate:           {(1 - stats['final']/max(stats['total'],1))*100:>7.1f}%")
    print("=" * 60)
    
    if stats['final'] == 0:
        print("\nWarning: No CNV segments passed filtering!")


if __name__ == "__main__":
    main()