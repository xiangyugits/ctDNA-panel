#!/usr/bin/env python3
"""
prepare_gene_annotation.py — 从 ANNOVAR refGene 数据库生成基因注释 BED

输入：
  1. 3 列 target BED（chr, start, end，无基因名）
  2. ANNOVAR humandb 目录中的 refGene.txt（或任意 refGene 格式文件）

输出：
  4 列基因注释 BED（chr, start, end, gene_name）
  可直接用作 panelcn.MOPS 的 --gene-annotation 参数

算法：
  - 从 refGene 提取每个基因的全长坐标（所有转录本的 min start → max end）
  - 对每个 amplicon 找到重叠的基因，重叠 bp 最大的优先
  - 若一个 amplicon 跨越多个基因，输出重叠最大者

用法：
  python3 prepare_gene_annotation.py \
      --target-bed  /path/to/loci.bed \
      --refgene     /path/to/annovar/humandb/hg19_refGene.txt \
      --output      /path/to/gene_annotation.bed

替代方案（若已安装 bedtools）：
  bedtools intersect \
      -a loci.bed \
      -b <(awk -F'\t' '{print $3"\t"$5"\t"$6"\t"$13}' hg19_refGene.txt | sort -k1,1 -k2,2n) \
      -wa -wb | \
      awk '{print $1"\t"$2"\t"$3"\t"$7}' | sort -u > gene_annotation.bed
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_refgene(path, genome="hg19"):
    """
    解析 ANNOVAR refGene.txt，返回 {gene_name: (chr, start, end)}
    
    refGene.txt 格式（ANNOVAR）：
      col 1:  bin
      col 2:  transcript name (NM_...)
      col 3:  chrom
      col 4:  strand
      col 5:  txStart  (0-based)
      col 6:  txEnd
      col 7:  cdsStart
      col 8:  cdsEnd
      col 9:  exonCount
      col 10: exonStarts
      col 11: exonEnds
      col 12: score
      col 13: name2 (gene symbol)  ← 目标列
      col 14: cdsStartStat
      col 15: cdsEndStat
      col 16: exonFrames
    
    返回每个基因的合并坐标（取所有转录本的 min start, max end）
    """
    gene_coords = {}  # gene_name → (chr, min_start, max_end)
    n_transcripts = 0
    
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 13:
                continue
            
            try:
                chrom = cols[2]
                tx_start = int(cols[4])  # 0-based
                tx_end = int(cols[5])
                gene_name = cols[12]  # name2
            except (ValueError, IndexError):
                continue
            
            n_transcripts += 1
            
            if gene_name not in gene_coords:
                gene_coords[gene_name] = [chrom, tx_start, tx_end]
            else:
                prev = gene_coords[gene_name]
                if chrom != prev[0]:
                    # 同名基因在多个染色体上（罕见），保留首次出现的
                    continue
                prev[1] = min(prev[1], tx_start)
                prev[2] = max(prev[2], tx_end)
    
    print(f"  解析 refGene: {n_transcripts} 个转录本 → {len(gene_coords)} 个基因", file=sys.stderr)
    return gene_coords


def normalise_chr(chrom):
    """统一 chr 前缀：确保都以 'chr' 开头"""
    chrom = chrom.strip()
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom


def overlap_length(a_start, a_end, b_start, b_end):
    """计算两个区间的重叠长度（0-based，半开区间）"""
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    return max(0, overlap_end - overlap_start)


def assign_genes(target_bed, gene_coords):
    """
    为每个 amplicon 分配基因名
    
    策略：
      1. 找到所有与该 amplicon 有重叠的基因
      2. 选重叠 bp 最多的基因
      3. 若无重叠基因，选距离最近的基因（仅在 amplicon 完全在基因间区时）
    """
    results = []
    n_multi = 0
    n_none = 0
    
    for i, (chrom, start, end) in enumerate(target_bed):
        chrom_norm = normalise_chr(chrom)
        best_gene = None
        best_overlap = 0
        all_overlaps = []
        
        for gene_name, (g_chr, g_start, g_end) in gene_coords.items():
            if normalise_chr(g_chr) != chrom_norm:
                continue
            ov = overlap_length(start, end, g_start, g_end)
            if ov > 0:
                all_overlaps.append((gene_name, ov))
                if ov > best_overlap:
                    best_overlap = ov
                    best_gene = gene_name
        
        if best_gene is None:
            # 无重叠：找距离最近的基因
            n_none += 1
            min_dist = float("inf")
            for gene_name, (g_chr, g_start, g_end) in gene_coords.items():
                if normalise_chr(g_chr) != chrom_norm:
                    continue
                if end <= g_start:
                    dist = g_start - end
                elif start >= g_end:
                    dist = start - g_end
                else:
                    dist = 0
                if dist < min_dist:
                    min_dist = dist
                    best_gene = gene_name
        
        if len(all_overlaps) > 1:
            n_multi += 1
        
        results.append((chrom, start, end, best_gene or "intergenic"))
    
    print(f"  基因分配: {len(results)} amplicons", file=sys.stderr)
    print(f"    无重叠（归入最近基因）: {n_none}", file=sys.stderr)
    print(f"    多基因重叠（已选最大重叠者）: {n_multi}", file=sys.stderr)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="从 ANNOVAR refGene 数据库生成 panel 基因注释 BED"
    )
    parser.add_argument("--target-bed", required=True,
                        help="3 列 target BED（chr, start, end）")
    parser.add_argument("--refgene", required=True,
                        help="ANNOVAR refGene.txt 文件路径（humandb/hg19_refGene.txt）")
    parser.add_argument("--output", required=True,
                        help="输出 4 列基因注释 BED（chr, start, end, gene_name）")
    parser.add_argument("--genome", default="hg19",
                        help="参考基因组版本（仅用于日志，默认 hg19）")
    args = parser.parse_args()
    
    # -----------------------------------------------------------------
    # Step 1: 读取 target BED
    # -----------------------------------------------------------------
    target_amplicons = []
    with open(args.target_bed, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                chrom = cols[0].strip()
                start = int(cols[1])
                end = int(cols[2])
                target_amplicons.append((chrom, start, end))
            except ValueError:
                continue
    
    print(f"Target BED: {len(target_amplicons)} amplicons", file=sys.stderr)
    
    # -----------------------------------------------------------------
    # Step 2: 解析 refGene
    # -----------------------------------------------------------------
    gene_coords = parse_refgene(args.refgene, args.genome)
    
    # -----------------------------------------------------------------
    # Step 3: 分配基因名
    # -----------------------------------------------------------------
    annotated = assign_genes(target_amplicons, gene_coords)
    
    # -----------------------------------------------------------------
    # Step 4: 写出
    # -----------------------------------------------------------------
    with open(args.output, "w") as f:
        for chrom, start, end, gene in annotated:
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
    
    # 统计基因覆盖
    gene_counts = defaultdict(int)
    for _, _, _, gene in annotated:
        gene_counts[gene] += 1
    
    print(f"\n输出: {args.output}", file=sys.stderr)
    print(f"  覆盖基因: {len(gene_counts)}", file=sys.stderr)
    print(f"  Amplicon/基因 Top5:", file=sys.stderr)
    for gene, count in sorted(gene_counts.items(), key=lambda x: -x[1])[:5]:
        print(f"    {gene}: {count} amplicons", file=sys.stderr)


if __name__ == "__main__":
    main()
