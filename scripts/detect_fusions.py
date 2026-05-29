#!/usr/bin/env python3
"""
Targeted Fusion Detection for Amplicon-based ctDNA Panels
==========================================================

原理：在 Panel 目标区域内搜索 discordant read pairs + split reads，
     通过聚类识别候选融合断点。

适用场景：
  - 扩增子 Panel（无 off-target reads，无全基因组 insert size 模型）
  - 仅检测 Panel 实际覆盖的基因区域内的融合事件

输出：
  1. .fusions.tsv  — 基因级融合报告（基因对、断点位置、支持reads数）
  2. .fusions.vcf  — 断点级 VCF（可用 IGV 可视化验证）
  3. .summary.txt  — 人类可读摘要

依赖：pip install pysam
"""

import argparse
import sys
import os
from collections import defaultdict
from itertools import combinations


def parse_bed(bed_path):
    """解析 BED 文件，返回 {chr: [(start, end, gene), ...]}"""
    regions = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3] if len(parts) >= 4 else "."
            regions[chrom].append((start, end, gene))
    return regions


def expand_regions(regions, flank):
    """将 BED 区域外扩 flank bp，用于捕获跨区域断点"""
    expanded = defaultdict(list)
    for chrom, intervals in regions.items():
        for start, end, gene in intervals:
            expanded[chrom].append((max(0, start - flank), end + flank, gene))
    return expanded


def region_contains(regions, chrom, pos):
    """检查位置是否在目标区域内（模糊匹配）"""
    if chrom not in regions:
        return False, None
    for start, end, gene in regions[chrom]:
        if start <= pos <= end:
            return True, gene
    return False, None


def find_gene_for_region(regions, chrom, pos):
    """查找 pos 对应的基因"""
    if chrom not in regions:
        return None
    for start, end, gene in regions[chrom]:
        if start <= pos <= end:
            return gene
    return None


def detect_fusions(bam_path, regions, expanded_regions,
                   min_mapq=20, min_support=2, flank=500):
    """
    核心检测逻辑：
    1. 遍历 BAM 文件所有 reads
    2. 找出 discordant read pairs（两条 read 映射到不同 chr 或不同基因）
    3. 找出含 soft-clip 的 split reads
    4. 聚类断点坐标
    5. 过滤低支持度的候选
    """
    try:
        import pysam
    except ImportError:
        print("ERROR: pysam is required. Install: pip install pysam", file=sys.stderr)
        sys.exit(1)

    bam = pysam.AlignmentFile(bam_path, "rb")

    # --- 收集 discordant pairs 和 split reads ---
    discordant_pairs = []    # [(chrA, posA, geneA, chrB, posB, geneB, read_id)]
    split_reads = []         # [(chr, pos, gene, clip_seq, clip_len, read_id)]

    read_count = 0
    for read in bam.fetch():
        read_count += 1
        if read_count % 1000000 == 0:
            print(f"  Processed {read_count:,} reads...", file=sys.stderr)

        # 质控
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        if read.is_qcfail:
            continue

        chrom = read.reference_name
        pos = read.reference_start
        mate_chrom = read.next_reference_name

        in_region, gene = region_contains(expanded_regions, chrom, pos)
        if not in_region:
            continue

        # --- Split reads (含 soft-clip) ---
        if read.cigartuples:
            first_op, first_len = read.cigartuples[0]
            last_op, last_len = read.cigartuples[-1]
            if first_op == 4 and first_len >= 20:  # 5' soft-clip >= 20bp
                split_reads.append((chrom, pos, gene,
                                    read.query_sequence[:first_len],
                                    first_len, read.query_name))
            if last_op == 4 and last_len >= 20:  # 3' soft-clip >= 20bp
                clip_start = len(read.query_sequence) - last_len
                split_reads.append((chrom, pos, gene,
                                    read.query_sequence[clip_start:],
                                    last_len, read.query_name))

        # --- Discordant read pairs ---
        if read.is_paired and read.is_proper_pair:
            continue  # 跳过正常的 concordant pairs

        if read.is_read1:  # 只处理 read1 避免重复计数
            mate_in_region, mate_gene = region_contains(
                expanded_regions, mate_chrom, read.next_reference_start)

            # 判断是否为"discordant"：不同 chr 或同 chr 但不同基因
            is_discordant = False
            if chrom != mate_chrom:
                is_discordant = True
            elif mate_in_region and gene != mate_gene and gene and mate_gene:
                # 同染色体不同基因
                distance = abs(pos - read.next_reference_start)
                if distance > 10000:  # >10kb 才考虑
                    is_discordant = True

            if is_discordant:
                discordant_pairs.append((
                    chrom, pos, gene,
                    mate_chrom, read.next_reference_start, mate_gene or ".",
                    read.query_name
                ))

    bam.close()
    print(f"  Total reads processed: {read_count:,}", file=sys.stderr)
    print(f"  Discordant pairs:      {len(discordant_pairs):,}", file=sys.stderr)
    print(f"  Split reads:           {len(split_reads):,}", file=sys.stderr)

    # --- 聚类断点 ---
    fusions = cluster_breakpoints(discordant_pairs, split_reads,
                                  regions, min_support)
    return fusions


def cluster_breakpoints(discordant_pairs, split_reads,
                        regions, min_support):
    """
    聚类断点坐标：
      - 按 (chrA, geneA, chrB, geneB) 分组
      - 每组内按 position 聚类（同一断点 +- 200bp）
      - 合并 discordant pairs + split reads 支持
    """
    # --- 按基因对分组 ---
    gene_pairs = defaultdict(lambda: {
        "breakpoints": defaultdict(list),  # (chrA, posA_cluster, chrB, posB_cluster) -> [reads]
        "split_reads": [],
        "discordant": 0
    })

    # Discordant pairs → 基因对
    for (chrA, posA, geneA, chrB, posB, geneB, read_id) in discordant_pairs:
        key = tuple(sorted([geneA, geneB]))
        gene_pairs[key]["discordant"] += 1
        # 聚类：500bp 窗口
        posA_bin = (posA // 500) * 500
        posB_bin = (posB // 500) * 500
        bp_key = (chrA, posA_bin, chrB, posB_bin)
        gene_pairs[key]["breakpoints"][bp_key].append(read_id)

    # Split reads → 基因对（按 clip 位置的基因判断）
    for (chr_s, pos_s, gene_s, clip_seq, clip_len, read_id) in split_reads:
        # 单端 split read，配对到最近的目标区域基因
        # 简化处理：记录到 gene_s 的自融合
        pass  # split reads 单独处理较复杂，先做 discordant 聚类

    # --- 过滤 + 构建融合列表 ---
    fusion_results = []
    for (geneA, geneB), data in sorted(gene_pairs.items()):
        # 合并所有断点的支持 reads
        total_unique_reads = set()
        for bp_key, read_ids in data["breakpoints"].items():
            total_unique_reads.update(read_ids)

        support_count = len(total_unique_reads)
        if support_count < min_support:
            continue

        # 选择最主要的断点
        best_bp = max(data["breakpoints"].items(),
                      key=lambda x: len(x[1]))
        (chrA, posA_cluster, chrB, posB_cluster), _ = best_bp

        fusion_results.append({
            "geneA": geneA,
            "geneB": geneB,
            "chrA": chrA,
            "posA": posA_cluster,
            "chrB": chrB,
            "posB": posB_cluster,
            "supporting_reads": support_count,
            "discordant_pairs": data["discordant"],
            "distinct_breakpoints": len(data["breakpoints"])
        })

    # 按支持 reads 数排序
    fusion_results.sort(key=lambda x: x["supporting_reads"], reverse=True)
    return fusion_results


def write_outputs(fusions, output_fusions, output_vcf, output_summary):
    """输出三种格式的结果文件"""

    # --- TSV 格式（基因级融合报告）---
    with open(output_fusions, 'w') as fh:
        fh.write("\t".join([
            "Gene_A", "Gene_B", "Chr_A", "Breakpoint_A",
            "Chr_B", "Breakpoint_B", "Supporting_Reads",
            "Discordant_Pairs", "Distinct_Breakpoints", "Fusion_Type"
        ]) + "\n")

        for f in fusions:
            fusion_type = "INTER_CHR" if f["chrA"] != f["chrB"] else "INTRA_CHR"
            fh.write("\t".join([
                f["geneA"], f["geneB"],
                f["chrA"], str(f["posA"]),
                f["chrB"], str(f["posB"]),
                str(f["supporting_reads"]),
                str(f["discordant_pairs"]),
                str(f["distinct_breakpoints"]),
                fusion_type
            ]) + "\n")

    # --- VCF 格式（每条融合一个记录）---
    with open(output_vcf, 'w') as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=detect_fusions.py\n")
        fh.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        fh.write("##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakend\">\n")
        fh.write("##INFO=<ID=GENE_A,Number=1,Type=String,Description=\"Gene A\">\n")
        fh.write("##INFO=<ID=GENE_B,Number=1,Type=String,Description=\"Gene B\">\n")
        fh.write("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Supporting reads\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for i, f in enumerate(fusions):
            fusion_id = f"FUSION_{i + 1}"
            mate_id = f"FUSION_{i + 1}_MATE"
            svtype = "BND"
            alt = f"N[{f['chrB']}:{f['posB']}["
            info = (f"SVTYPE={svtype};MATEID={mate_id};"
                    f"GENE_A={f['geneA']};GENE_B={f['geneB']};"
                    f"SUPPORT={f['supporting_reads']}")
            fh.write(f"{f['chrA']}\t{f['posA']}\t{fusion_id}\tN\t{alt}\t.\tPASS\t{info}\n")

            # Mate breakend
            alt_mate = f"N[{f['chrA']}:{f['posA']}["
            info_mate = (f"SVTYPE={svtype};MATEID={fusion_id};"
                         f"GENE_A={f['geneA']};GENE_B={f['geneB']};"
                         f"SUPPORT={f['supporting_reads']}")
            fh.write(f"{f['chrB']}\t{f['posB']}\t{mate_id}\tN\t{alt_mate}\t.\tPASS\t{info_mate}\n")

    # --- 文本摘要 ---
    with open(output_summary, 'w') as fh:
        fh.write("Targeted Fusion Detection Summary\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Total candidate fusions: {len(fusions)}\n\n")

        inter_chr = sum(1 for f in fusions if f["chrA"] != f["chrB"])
        intra_chr = sum(1 for f in fusions if f["chrA"] == f["chrB"])
        fh.write(f"  Inter-chromosomal: {inter_chr}\n")
        fh.write(f"  Intra-chromosomal: {intra_chr}\n\n")

        if fusions:
            fh.write("Top Candidates:\n")
            fh.write("-" * 60 + "\n")
            for f in fusions[:10]:
                ftype = "INTER" if f["chrA"] != f["chrB"] else "INTRA"
                fh.write(
                    f"  {f['geneA']}--{f['geneB']} [{ftype}] "
                    f"support={f['supporting_reads']} "
                    f"({f['chrA']}:{f['posA']} <-> {f['chrB']}:{f['posB']})\n"
                )


def main():
    parser = argparse.ArgumentParser(
        description="Targeted Fusion Detection for Amplicon ctDNA Panels"
    )
    parser.add_argument("--bam", required=True,
                        help="Tumor BAM file (deduplicated)")
    parser.add_argument("--bed", required=True,
                        help="Panel target regions BED file")
    parser.add_argument("--output-fusions", required=True,
                        help="Output fusion TSV file")
    parser.add_argument("--output-vcf", required=True,
                        help="Output fusion VCF file")
    parser.add_argument("--output-summary", required=True,
                        help="Output summary text file")
    parser.add_argument("--min-mapq", type=int, default=20,
                        help="Minimum mapping quality (default: 20)")
    parser.add_argument("--min-support", type=int, default=2,
                        help="Minimum supporting reads per fusion (default: 2)")
    parser.add_argument("--flank", type=int, default=500,
                        help="BED region flank extension in bp (default: 500)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads (reserved for future use)")
    args = parser.parse_args()

    # 验证输入
    for path, name in [(args.bam, "BAM"), (args.bed, "BED")]:
        if not os.path.exists(path):
            print(f"ERROR: {name} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    print(f"Loading target regions from: {args.bed}", file=sys.stderr)
    regions = parse_bed(args.bed)
    expanded_regions = expand_regions(regions, args.flank)

    total_intervals = sum(len(v) for v in regions.values())
    print(f"  Target regions: {total_intervals} intervals across "
          f"{len(regions)} chromosomes", file=sys.stderr)
    print(f"  Flank extension: {args.flank} bp", file=sys.stderr)
    print(f"  Min MAPQ: {args.min_mapq}", file=sys.stderr)
    print(f"  Min support: {args.min_support}", file=sys.stderr)
    print("", file=sys.stderr)

    print(f"Scanning BAM: {args.bam}", file=sys.stderr)
    fusions = detect_fusions(
        args.bam, regions, expanded_regions,
        min_mapq=args.min_mapq,
        min_support=args.min_support,
        flank=args.flank
    )

    print(f"\nWriting results...", file=sys.stderr)
    write_outputs(fusions, args.output_fusions,
                  args.output_vcf, args.output_summary)

    print(f"Done. {len(fusions)} candidate fusions detected.", file=sys.stderr)
    print(f"  Fusions TSV:    {args.output_fusions}", file=sys.stderr)
    print(f"  Fusions VCF:    {args.output_vcf}", file=sys.stderr)
    print(f"  Summary:        {args.output_summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
