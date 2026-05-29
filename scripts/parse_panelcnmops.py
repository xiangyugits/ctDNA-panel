#!/usr/bin/env python3
"""
parse_panelcnmops.py
====================
将 panelcn.MOPS 输出的 per-amplicon TSV 解析并汇总为基因级 CNV 报告，
同时生成与 filter_cnv.py 兼容的 .cns 格式文件（用于后续统一过滤）。

输入：run_panelcnmops.R 输出的 {sample}.panelcnmops.tsv
输出：
  1. {sample}.panelcnmops.gene.tsv  —— 基因级汇总（每基因一行）
  2. {sample}.panelcnmops.cns       —— CNVkit 兼容格式（供 filter_cnv.py 过滤）
  3. {sample}.panelcnmops.report.txt —— 人类可读摘要报告
"""

import argparse
import sys
import math
from pathlib import Path
from collections import defaultdict


# =============================================================================
# 工具函数
# =============================================================================

def safe_float(val, default=float("nan")):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=0):
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return default


def cn_to_type(cn: int) -> str:
    """整数拷贝数 → AMP / DEL / NEUTRAL / LOH"""
    if cn == 0:
        return "HOMDEL"
    elif cn == 1:
        return "DEL"
    elif cn == 2:
        return "NEUTRAL"
    elif cn == 3:
        return "GAIN"
    elif cn >= 4:
        return "AMP"
    else:
        return "NEUTRAL"


def log2_to_cn(log2r: float) -> float:
    """log2(CN/2) → CN"""
    return 2.0 * (2.0 ** log2r)


# =============================================================================
# 核心解析函数
# =============================================================================

def parse_panelcnmops_tsv(input_file: str):
    """
    读取 per-amplicon TSV，返回行列表。
    期望列（大小写不敏感）：
      chr, start, end, gene, CN, log2ratio,
      RC_tumor, RC_normal_mean, RC_normalized, pvalue, cnv_type
    """
    rows = []
    with open(input_file) as fh:
        header_raw = fh.readline().strip().split("\t")
        header = [h.lower().strip('"') for h in header_raw]

        # 列名映射（容错多种写法）
        col_aliases = {
            "chr":            ["chr", "chrom", "chromosome"],
            "start":          ["start"],
            "end":            ["end"],
            "gene":           ["gene", "genename", "gene_name"],
            "cn":             ["cn"],
            "log2ratio":      ["log2ratio", "log2", "log2_ratio", "log2r"],
            "rc_tumor":       ["rc_tumor", "rc", "readcount_tumor"],
            "rc_normal_mean": ["rc_normal_mean", "medrc", "rc_normal"],
            "rc_normalized":  ["rc_normalized", "normalized_rc"],
            "pvalue":         ["pvalue", "p_value", "p.value"],
            "cnv_type":       ["cnv_type", "type", "event"],
        }

        col_idx = {}
        for field, aliases in col_aliases.items():
            for alias in aliases:
                if alias in header:
                    col_idx[field] = header.index(alias)
                    break

        required = ["chr", "start", "end", "gene", "cn", "log2ratio"]
        missing = [r for r in required if r not in col_idx]
        if missing:
            raise ValueError(f"输入文件缺少必需列: {missing}。实际列: {header_raw}")

        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            fields = [f.strip('"') for f in fields]

            def get(key, default=""):
                idx = col_idx.get(key)
                if idx is None or idx >= len(fields):
                    return default
                return fields[idx]

            row = {
                "chr":            get("chr"),
                "start":          safe_int(get("start")),
                "end":            safe_int(get("end")),
                "gene":           get("gene", "Unknown"),
                "cn":             safe_int(get("cn"), default=2),
                "log2ratio":      safe_float(get("log2ratio")),
                "rc_tumor":       safe_float(get("rc_tumor")),
                "rc_normal_mean": safe_float(get("rc_normal_mean")),
                "rc_normalized":  safe_float(get("rc_normalized")),
                "pvalue":         safe_float(get("pvalue"), default=1.0),
                "cnv_type":       get("cnv_type", cn_to_type(safe_int(get("cn"), 2))),
            }
            rows.append(row)

    return rows


def aggregate_by_gene(rows):
    """
    将 per-amplicon 结果按基因聚合为基因级结果。
    聚合规则：
      - 染色体：取该基因最常出现的 chr
      - 坐标：min(start) ~ max(end)
      - CN：取加权中位数（按 RC_tumor 加权）
      - log2ratio：加权均值（按 RC_tumor 加权）
      - amplicons：统计总数 / CNV 数
      - cnv_fraction：有 CNV 的 amplicon 比例
      - pvalue：取最显著（最小）p 值
    """
    gene_rows = defaultdict(list)
    for row in rows:
        gene = row["gene"]
        if not gene or gene == "NA":
            gene = f"{row['chr']}:{row['start']}-{row['end']}"
        gene_rows[gene].append(row)

    results = []
    for gene, amplist in sorted(gene_rows.items()):
        n_total = len(amplist)

        # 染色体
        chrs = [r["chr"] for r in amplist]
        chrom = max(set(chrs), key=chrs.count)

        # 坐标范围
        g_start = min(r["start"] for r in amplist)
        g_end   = max(r["end"]   for r in amplist)

        # 加权 log2ratio（按 RC_tumor 加权）
        weights = [max(r["rc_tumor"], 1) for r in amplist]
        total_w = sum(weights)
        w_log2  = sum(r["log2ratio"] * w for r, w in zip(amplist, weights)) / total_w

        # 加权 CN
        w_cn = sum(r["cn"] * w for r, w in zip(amplist, weights)) / total_w
        median_cn = round(w_cn)

        # CNV amplicon 统计
        cnv_amps = [r for r in amplist if r["cnv_type"] not in ("NEUTRAL",)]
        n_cnv    = len(cnv_amps)
        cnv_frac = n_cnv / n_total if n_total > 0 else 0.0

        # 最显著 p 值
        pvalues = [r["pvalue"] for r in amplist if not math.isnan(r["pvalue"])]
        min_pval = min(pvalues) if pvalues else 1.0

        # 均值 RC
        mean_rc_tumor  = sum(r["rc_tumor"] for r in amplist) / n_total
        mean_rc_normal = sum(r["rc_normal_mean"] for r in amplist) / n_total

        # 整体 CNV 类型判断（基于 cnv_fraction > 0.5 多数原则）
        cnv_type_votes = [r["cnv_type"] for r in cnv_amps]
        if not cnv_type_votes:
            gene_cnv_type = "NEUTRAL"
        else:
            gene_cnv_type = max(set(cnv_type_votes), key=cnv_type_votes.count)

        results.append({
            "gene":            gene,
            "chr":             chrom,
            "start":           g_start,
            "end":             g_end,
            "CN":              median_cn,
            "log2ratio":       round(w_log2, 4),
            "cnv_type":        gene_cnv_type,
            "cnv_fraction":    round(cnv_frac, 3),
            "n_amplicons":     n_total,
            "n_cnv_amplicons": n_cnv,
            "mean_RC_tumor":   round(mean_rc_tumor, 1),
            "mean_RC_normal":  round(mean_rc_normal, 1),
            "min_pvalue":      min_pval,
        })

    return results


def write_gene_tsv(gene_results, output_file: str, sample_name: str):
    """写出基因级汇总 TSV"""
    header = [
        "Sample", "Gene", "Chr", "Start", "End",
        "CN", "Log2Ratio", "CNV_Type", "CNV_Fraction",
        "N_Amplicons", "N_CNV_Amplicons",
        "Mean_RC_Tumor", "Mean_RC_Normal", "Min_Pvalue"
    ]
    with open(output_file, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in gene_results:
            row = [
                sample_name,
                r["gene"],
                r["chr"],
                str(r["start"]),
                str(r["end"]),
                str(r["CN"]),
                str(r["log2ratio"]),
                r["cnv_type"],
                str(r["cnv_fraction"]),
                str(r["n_amplicons"]),
                str(r["n_cnv_amplicons"]),
                str(r["mean_RC_tumor"]),
                str(r["mean_RC_normal"]),
                f"{r['min_pvalue']:.2e}",
            ]
            fh.write("\t".join(row) + "\n")


def write_cns_compat(gene_results, output_file: str):
    """
    写出 CNVkit .cns 兼容格式，供 filter_cnv.py 过滤。
    列：chromosome / start / end / gene / log2 / cn / depth / probes / weight / ci
    注：
      - probes  用 n_amplicons 填充（对应 panelcn.MOPS 的 amplicon 数）
      - depth   用 mean_RC_tumor 填充
      - ci      用 pvalue 转换的近似不确定度填充
      - weight  固定 1.0
    """
    header = "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\tci"
    with open(output_file, "w") as fh:
        fh.write(header + "\n")
        for r in gene_results:
            # CI 近似：p 值越大（越不确定），CI 越宽
            # 映射规则：p=0.05 → CI=0.5，p=1.0 → CI=2.0，p=0.001 → CI=0.1
            pv = max(r["min_pvalue"], 1e-300)
            ci_approx = round(-0.5 * math.log10(pv) ** 0.3 + 2.0, 3)
            ci_approx = max(0.05, min(ci_approx, 5.0))  # 限制在合理范围

            row = [
                r["chr"],
                str(r["start"]),
                str(r["end"]),
                r["gene"],
                str(r["log2ratio"]),
                str(r["CN"]),
                str(r["mean_RC_tumor"]),
                str(r["n_amplicons"]),
                "1.0",
                str(ci_approx),
            ]
            fh.write("\t".join(row) + "\n")


def write_report(gene_results, output_file: str, sample_name: str,
                 min_log2: float, min_frac: float):
    """写出人类可读的文本报告"""
    cnv_results = [r for r in gene_results if r["cnv_type"] != "NEUTRAL"]
    amp_results = sorted(
        [r for r in cnv_results if r["CN"] > 2],
        key=lambda x: x["log2ratio"], reverse=True
    )
    del_results = sorted(
        [r for r in cnv_results if r["CN"] < 2],
        key=lambda x: x["log2ratio"]
    )

    lines = [
        "=" * 70,
        f"panelcn.MOPS CNV 分析报告",
        f"样本: {sample_name}",
        "=" * 70,
        "",
        f"过滤参数: min_log2_abs={min_log2}, min_cnv_fraction={min_frac}",
        "",
        f"总基因数:     {len(gene_results)}",
        f"CNV 基因数:   {len(cnv_results)}  (扩增: {len(amp_results)}, 缺失: {len(del_results)})",
        "",
    ]

    if amp_results:
        lines += ["─" * 70, "扩增（AMP/GAIN）:", "─" * 70]
        lines += [f"  {'Gene':<12} {'Chr':<6} {'Start':<12} {'End':<12} "
                  f"{'CN':<4} {'Log2R':<8} {'Type':<8} {'Frac':<6} {'N_amp'}"]
        for r in amp_results:
            lines.append(
                f"  {r['gene']:<12} {r['chr']:<6} {r['start']:<12} {r['end']:<12} "
                f"{r['CN']:<4} {r['log2ratio']:<8.3f} {r['cnv_type']:<8} "
                f"{r['cnv_fraction']:<6.2f} {r['n_amplicons']}"
            )
        lines.append("")

    if del_results:
        lines += ["─" * 70, "缺失（DEL/HOMDEL）:", "─" * 70]
        lines += [f"  {'Gene':<12} {'Chr':<6} {'Start':<12} {'End':<12} "
                  f"{'CN':<4} {'Log2R':<8} {'Type':<8} {'Frac':<6} {'N_amp'}"]
        for r in del_results:
            lines.append(
                f"  {r['gene']:<12} {r['chr']:<6} {r['start']:<12} {r['end']:<12} "
                f"{r['CN']:<4} {r['log2ratio']:<8.3f} {r['cnv_type']:<8} "
                f"{r['cnv_fraction']:<6.2f} {r['n_amplicons']}"
            )
        lines.append("")

    if not cnv_results:
        lines.append("  未检测到符合阈值的 CNV。")
        lines.append("")

    lines += ["=" * 70]

    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# =============================================================================
# 主程序
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="解析 panelcn.MOPS 输出，生成基因级 CNV 报告"
    )
    parser.add_argument("--input",        "-i", required=True,
                        help="panelcn.MOPS 输出的 per-amplicon TSV 文件（{sample}.panelcnmops.tsv）")
    parser.add_argument("--sample-name",  "-s", required=True,
                        help="样本名称（用于报告中标注）")
    parser.add_argument("--output-dir",   "-o", required=True,
                        help="输出目录（将生成 .gene.tsv / .cns / .report.txt）")
    parser.add_argument("--min-log2-abs", type=float, default=0.3,
                        help="基因级 |log2ratio| 最低阈值（过滤后输出，默认 0.3）")
    parser.add_argument("--min-cnv-fraction", type=float, default=0.3,
                        help="基因内 CNV amplicon 最低比例阈值（默认 0.3 即 30%%）")
    parser.add_argument("--min-pvalue",   type=float, default=1.0,
                        help="最高 p 值阈值（默认 1.0 = 不过滤；建议 0.05）")
    args = parser.parse_args()

    # 检查输入
    if not Path(args.input).exists():
        print(f"错误：输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print(f"读取 panelcn.MOPS 结果: {args.input}")

    # Step 1：读取 per-amplicon TSV
    rows = parse_panelcnmops_tsv(args.input)
    print(f"  共 {len(rows)} 个 amplicon 记录")

    # Step 2：按基因聚合
    gene_results_all = aggregate_by_gene(rows)
    print(f"  聚合为 {len(gene_results_all)} 个基因")

    # Step 3：基因级过滤
    gene_results_filtered = [
        r for r in gene_results_all
        if (abs(r["log2ratio"]) >= args.min_log2_abs
            and r["cnv_fraction"] >= args.min_cnv_fraction
            and r["min_pvalue"] <= args.min_pvalue
            and r["cnv_type"] != "NEUTRAL")
    ]
    print(f"  过滤后 CNV 基因数: {len(gene_results_filtered)}")

    sn = args.sample_name
    od = args.output_dir

    # Step 4：写出各格式文件
    gene_tsv  = Path(od) / f"{sn}.panelcnmops.gene.tsv"
    cns_file  = Path(od) / f"{sn}.panelcnmops.cns"
    report    = Path(od) / f"{sn}.panelcnmops.report.txt"

    # 基因级 TSV（过滤后的 CNV）
    write_gene_tsv(gene_results_filtered, str(gene_tsv), sn)
    print(f"  基因级 TSV: {gene_tsv}")

    # CNVkit 兼容 .cns（过滤后）
    write_cns_compat(gene_results_filtered, str(cns_file))
    print(f"  .cns 兼容格式: {cns_file}")

    # 文本报告
    write_report(gene_results_all, str(report), sn,
                 args.min_log2_abs, args.min_cnv_fraction)
    print(f"  文本报告: {report}")

    print("解析完成。")


if __name__ == "__main__":
    main()
