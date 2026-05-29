#!/usr/bin/env Rscript
# =============================================================================
# run_panelcnmops.R
# 使用 panelcn.MOPS 对扩增子 ctDNA Panel 进行 CNV 检测
#
# 适用场景：
#   - 扩增子测序（Amplicon-based Panel），无 off-target reads
#   - 单配对正常对照（tumor vs NC），或多正常样本 PoN 模式
#   - ctDNA / 液体活检低频变异
#
# 依赖：
#   panelcn.MOPS (Bioconductor), GenomicRanges, Rsamtools, BiocParallel
#
# 安装（如未安装）：
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("panelcn.MOPS")
# =============================================================================

suppressPackageStartupMessages({
  library(panelcn.MOPS)
  library(GenomicRanges)
  library(Rsamtools)
  library(BiocParallel)
  library(optparse)
})

# =============================================================================
# 命令行参数解析
# =============================================================================

option_list <- list(
  make_option("--tumor-bam",    type="character", dest="tumor_bam",
              help="肿瘤样本 BAM 文件路径（必须已建索引）[必需]"),
  make_option("--normal-bam",   type="character", dest="normal_bam",
              help="正常对照 BAM 文件路径，支持逗号分隔多个（PoN 模式）[必需]"),
  make_option("--target-bed",   type="character", dest="target_bed",
              help="Panel 目标区域 BED 文件路径（不要 gz 压缩）[必需]"),
  make_option("--sample-name",  type="character", dest="sample_name",
              help="肿瘤样本名称（用于输出文件命名）[必需]"),
  make_option("--output-dir",   type="character", dest="output_dir",
              help="输出目录路径 [必需]"),
  make_option("--genome",       type="character", dest="genome",    default="hg19",
              help="参考基因组版本，hg19 或 hg38 [默认: hg19]"),
  make_option("--min-rc",       type="integer",   dest="min_rc",    default=1,
              help="每个 amplicon 的最低 read count 阈值 [默认: 1]"),
  make_option("--alpha",        type="numeric",   dest="alpha",     default=0.05,
              help="CNV 检测显著性阈值 [默认: 0.05]"),
  make_option("--threads",      type="integer",   dest="threads",   default=4,
              help="并行线程数 [默认: 4]"),
  make_option("--log",          type="character", dest="log",       default="",
              help="日志输出文件（留空则输出到 stderr）")
)

opt <- parse_args(OptionParser(option_list=option_list,
                               description="panelcn.MOPS CNV detection for ctDNA amplicon panel"))

# 参数检查
for (req in c("tumor_bam", "normal_bam", "target_bed", "sample_name", "output_dir")) {
  if (is.null(opt[[req]])) {
    stop(sprintf("缺少必需参数 --%s，请运行 --help 查看帮助", gsub("_", "-", req)))
  }
}

# =============================================================================
# 日志配置
# =============================================================================

log_msg <- function(msg, level="INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", ts, level, msg)
  if (nchar(opt$log) > 0) {
    cat(line, "\n", file=opt$log, append=TRUE)
  } else {
    message(line)
  }
}

log_msg("panelcn.MOPS CNV 分析开始")
log_msg(sprintf("肿瘤样本:  %s", opt$tumor_bam))
log_msg(sprintf("正常对照:  %s", opt$normal_bam))
log_msg(sprintf("目标区域:  %s", opt$target_bed))
log_msg(sprintf("样本名称:  %s", opt$sample_name))
log_msg(sprintf("输出目录:  %s", opt$output_dir))

# =============================================================================
# 输出目录
# =============================================================================

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# 并行设置
# =============================================================================

if (.Platform$OS.type == "unix") {
  register(MulticoreParam(workers=opt$threads))
} else {
  register(SnowParam(workers=opt$threads))
}

# =============================================================================
# Step 1：加载并验证目标区域 BED 文件
# =============================================================================

log_msg("Step 1: 读取目标区域 BED 文件")

# panelcn.MOPS 要求 BED 为未压缩格式，列顺序: chr / start / end / name
bed_raw <- read.table(opt$target_bed, header=FALSE, sep="\t",
                      stringsAsFactors=FALSE, comment.char="#")

# 确保至少有 4 列（chr, start, end, name）
if (ncol(bed_raw) < 3) {
  stop("BED 文件至少需要 3 列（chr, start, end），建议有第 4 列 name")
}

# 构造标准 4 列 BED（若无 name 列则自动生成）
if (ncol(bed_raw) >= 4) {
  bed4 <- bed_raw[, 1:4]
} else {
  bed4 <- bed_raw[, 1:3]
  bed4$V4 <- paste0(bed_raw[,1], ":", bed_raw[,2], "-", bed_raw[,3])
}
colnames(bed4) <- c("chr", "start", "end", "name")

# 确保 chr 前缀统一（panelcn.MOPS 内部用 chr 前缀）
if (!grepl("^chr", bed4$chr[1])) {
  bed4$chr <- paste0("chr", bed4$chr)
}

# 写出临时标准化 BED（panelcn.MOPS 从文件读取）
tmp_bed <- file.path(opt$output_dir, paste0(opt$sample_name, "_panelcnmops_target.bed"))
write.table(bed4, tmp_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
log_msg(sprintf("  目标 amplicon 数量: %d", nrow(bed4)))

# =============================================================================
# Step 2：统计各 amplicon 的 read count
# =============================================================================

log_msg("Step 2: 统计各样本 amplicon read count")

# 解析正常对照 BAM 列表（支持逗号分隔多个）
normal_bams <- trimws(strsplit(opt$normal_bam, ",")[[1]])
normal_bams <- normal_bams[file.exists(normal_bams)]
if (length(normal_bams) == 0) {
  stop("未找到任何正常对照 BAM 文件，请检查路径")
}
log_msg(sprintf("  正常对照 BAM 数量: %d", length(normal_bams)))

if (!file.exists(opt$tumor_bam)) {
  stop(sprintf("肿瘤 BAM 文件不存在: %s", opt$tumor_bam))
}

# 所有 BAM 列表（panelcn.MOPS 需要一次性读入全部）
all_bams <- c(opt$tumor_bam, normal_bams)

# 使用 panelcn.MOPS 的 countBamListInBed 统计 read count
# 返回 RangedSummarizedExperiment 对象
log_msg("  正在统计 read counts（可能需要数分钟）...")

rc_obj <- tryCatch({
  countBamListInBed(
    bamfiles    = all_bams,
    bed         = tmp_bed,
    minMapq     = 1,              # 最低 MAPQ（ctDNA 低质量但有效 reads 都保留）
    filteredFlag = 1024           # 过滤 PCR 重复（flag 0x400）
  )
}, error = function(e) {
  log_msg(sprintf("countBamListInBed 出错: %s", conditionMessage(e)), "ERROR")
  stop(e)
})

log_msg(sprintf("  read count 统计完成，共 %d 个 amplicon，%d 个样本",
                nrow(rc_obj), ncol(rc_obj)))

# =============================================================================
# Step 3：区分 tumor / normal，构建分析对象
# =============================================================================

log_msg("Step 3: 构建 panelcn.MOPS 分析对象")

# 样本名称（取 BAM 文件名去后缀）
sample_names <- sapply(all_bams, function(b) {
  nm <- basename(b)
  gsub("\\.dedup\\.bam$|\\.bam$", "", nm)
})
colnames(rc_obj) <- sample_names

tumor_idx  <- 1L                          # 第一个是 tumor
normal_idx <- seq(2L, length(all_bams))   # 其余是 normal

log_msg(sprintf("  肿瘤样本:  %s (idx=%d)", sample_names[tumor_idx], tumor_idx))
log_msg(sprintf("  正常样本:  %s",
                paste(sample_names[normal_idx], collapse=", ")))

# =============================================================================
# Step 4：运行 panelcn.MOPS CNV 检测
# =============================================================================

log_msg("Step 4: 运行 panelcn.MOPS CNV 检测")

result <- tryCatch({
  panelcn.mops(
    input         = rc_obj,
    testiv        = tumor_idx,          # 指定肿瘤样本（列索引）
    I             = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),  # CN 状态先验
    minReadCount  = opt$min_rc,         # 最低 read count 阈值
    parallel      = opt$threads
  )
}, error = function(e) {
  log_msg(sprintf("panelcn.mops() 出错: %s", conditionMessage(e)), "ERROR")
  stop(e)
})

log_msg("  CNV 检测完成")

# =============================================================================
# Step 5：提取逐 amplicon CNV 结果
# =============================================================================

log_msg("Step 5: 提取 CNV 结果")

# 提取每个 amplicon 的详细信息
amplicon_res <- getResults(result, sampleName=sample_names[tumor_idx])

if (is.null(amplicon_res) || nrow(amplicon_res) == 0) {
  log_msg("  未检测到任何 CNV，输出空结果文件", "WARNING")
  # 写空文件（保持流程可继续）
  empty_out <- file.path(opt$output_dir, paste0(opt$sample_name, ".panelcnmops.tsv"))
  write.table(
    data.frame(
      chr=character(), start=integer(), end=integer(),
      gene=character(), CN=integer(), log2ratio=numeric(),
      RC_tumor=numeric(), RC_normal_mean=numeric(),
      RC_normalized=numeric(), pvalue=numeric(),
      cnv_type=character()
    ),
    empty_out, sep="\t", quote=FALSE, row.names=FALSE
  )
  log_msg(sprintf("  已写出空结果: %s", empty_out))
  quit(status=0)
}

log_msg(sprintf("  共 %d 个 amplicon 有 CNV 结果", nrow(amplicon_res)))

# 提取必要字段并标准化
out_df <- data.frame(
  chr              = as.character(seqnames(rowRanges(result))),
  start            = start(rowRanges(result)),
  end              = end(rowRanges(result)),
  gene             = as.character(rowRanges(result)$names),
  CN               = as.integer(amplicon_res$CN),
  log2ratio        = log2(pmax(amplicon_res$CN, 0.01) / 2),   # log2(CN/2) 与 CNVkit 格式对齐
  RC_tumor         = amplicon_res$RC,
  RC_normal_mean   = amplicon_res$medRC,
  RC_normalized    = amplicon_res$RC / pmax(amplicon_res$medRC, 1),
  pvalue           = amplicon_res$pvalue,
  stringsAsFactors = FALSE
)

# 标注 CNV 类型
out_df$cnv_type <- ifelse(out_df$CN > 2, "AMP",
                   ifelse(out_df$CN < 2, "DEL", "NEUTRAL"))

# chr 前缀处理：保留原始格式（与 BED 一致）
# 注：如需去除 chr 前缀可取消下行注释
# out_df$chr <- gsub("^chr", "", out_df$chr)

# =============================================================================
# Step 6：输出结果文件
# =============================================================================

log_msg("Step 6: 输出结果文件")

# 6a. 完整 per-amplicon 结果（所有 amplicon，含 NEUTRAL）
all_out <- file.path(opt$output_dir, paste0(opt$sample_name, ".panelcnmops.all.tsv"))
write.table(out_df, all_out, sep="\t", quote=FALSE, row.names=FALSE)
log_msg(sprintf("  全部 amplicon 结果: %s", all_out))

# 6b. 仅 CNV amplicon 结果（NEUTRAL 已过滤）
cnv_only <- out_df[out_df$cnv_type != "NEUTRAL", ]
cnv_out <- file.path(opt$output_dir, paste0(opt$sample_name, ".panelcnmops.tsv"))
write.table(cnv_only, cnv_out, sep="\t", quote=FALSE, row.names=FALSE)
log_msg(sprintf("  CNV amplicon 结果: %s (%d 个)", cnv_out, nrow(cnv_only)))

# 6c. 输出 CNV 统计摘要
n_amp   <- sum(out_df$cnv_type == "AMP")
n_del   <- sum(out_df$cnv_type == "DEL")
n_total <- nrow(out_df)
log_msg(sprintf("  统计: 总 amplicon=%d, 扩增=%d, 缺失=%d, 中性=%d",
                n_total, n_amp, n_del, n_total - n_amp - n_del))

# 6d. 保存 R 结果对象（用于后续可视化调试）
rds_out <- file.path(opt$output_dir, paste0(opt$sample_name, ".panelcnmops.rds"))
saveRDS(result, rds_out)
log_msg(sprintf("  R 结果对象: %s", rds_out))

# =============================================================================
# Step 7：可视化（每个基因一张图，存入子目录）
# =============================================================================

log_msg("Step 7: 生成可视化图（gene-level）")

plot_dir <- file.path(opt$output_dir, "plots")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

# 获取所有检测到 CNV 的基因（取 gene 字段的唯一值，去掉 NA/空值）
cnv_genes <- unique(na.omit(cnv_only$gene))
cnv_genes <- cnv_genes[cnv_genes != "" & !is.na(cnv_genes)]

if (length(cnv_genes) > 0) {
  log_msg(sprintf("  为 %d 个 CNV 基因生成图...", length(cnv_genes)))
  for (g in cnv_genes) {
    png_path <- file.path(plot_dir, sprintf("%s_%s_CNV.png", opt$sample_name, g))
    tryCatch({
      png(png_path, width=1200, height=600, res=120)
      plotSampleCNVs(result, sampleName=sample_names[tumor_idx], gene=g)
      dev.off()
    }, error = function(e) {
      log_msg(sprintf("  绘图失败 gene=%s: %s", g, conditionMessage(e)), "WARNING")
      if (dev.cur() > 1) dev.off()
    })
  }
  log_msg(sprintf("  图像已保存到: %s", plot_dir))
} else {
  log_msg("  无 CNV 基因，跳过可视化", "WARNING")
}

# =============================================================================
# 完成
# =============================================================================

log_msg("panelcn.MOPS 分析完成")
log_msg(sprintf("结果目录: %s", opt$output_dir))
