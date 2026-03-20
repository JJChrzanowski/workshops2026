# Pancreatic cancer drug repurposing
#
#   This is a TEACHING script. It is deliberately simple.
#   It generates a hypothesis. It is NOT a clinical decision tool.

# -----------------------------
# 1) Setup folders and URLs
# -----------------------------
outdir <- "pancreatic_workshop_output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

gse_file <- file.path(outdir, "GSE15471_series_matrix.txt.gz")
gpl_file <- file.path(outdir, "GPL570.annot.gz")

gse_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15471/matrix/GSE15471_series_matrix.txt.gz"
gpl_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

safe_download <- function(url, destfile) {
  if (!file.exists(destfile)) {
    message("Downloading: ", basename(destfile))
    download.file(url = url, destfile = destfile, mode = "wb", quiet = FALSE)
  } else {
    message("Using cached file: ", basename(destfile))
  }
}

safe_download(gse_url, gse_file)
safe_download(gpl_url, gpl_file)

# -----------------------------
# 2) Read GEO matrix and sample titles
# -----------------------------
extract_sample_titles <- function(matrix_path) {
  lines <- readLines(gzfile(matrix_path), warn = FALSE)
  title_line <- lines[grepl("^!Sample_title", lines)][1]
  if (is.na(title_line)) stop("Could not find !Sample_title line in the GEO matrix file.")
  titles <- strsplit(sub("^!Sample_title\\s*=\\s*", "", title_line), "\t")[[1]]
  titles <- gsub('^"|"$', "", titles)
  titles
}

# Reading GEO series matrix
expr_df <- read.delim(
  gzfile(gse_file),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "!",
  check.names = FALSE
)

sample_titles <- extract_sample_titles(gse_file)[-1]

probe_ids <- expr_df[[1]]
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gsub('"', '', probe_ids)
colnames(expr_mat) <- sample_titles

# Average technical replicate arrays by removing the _rep suffix
sample_keys <- sub("_rep$", "", sample_titles)
unique_keys <- unique(sample_keys)

expr_avg <- sapply(unique_keys, function(k) {
  idx <- which(sample_keys == k)
  if (length(idx) == 1) {
    expr_mat[, idx]
  } else {
    rowMeans(expr_mat[, idx, drop = FALSE])
  }
})

if (is.vector(expr_avg)) {
  expr_avg <- matrix(expr_avg, ncol = 1)
  rownames(expr_avg) <- rownames(expr_mat)
  colnames(expr_avg) <- unique_keys
}

# -----------------------------
# 3) Read GPL annotation and map probes -> gene symbols
# -----------------------------
read_gpl_annotation <- function(gpl_path) {
  lines <- readLines(gzfile(gpl_path), warn = FALSE)
  start_line <- grep("^ID\\t", lines)[1]
  if (is.na(start_line)) stop("Could not find the annotation table header in GPL570.annot.gz")
  annot_text <- lines[start_line:length(lines)]
  annot <- read.delim(
    textConnection(annot_text),
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE
  )
  annot
}

# Reading GPL annotation
annot <- read_gpl_annotation(gpl_file)
required_cols <- c("ID", "Gene symbol")

symbol_map <- annot$`Gene symbol`[match(rownames(expr_avg), annot$ID)]
symbol_map <- sub("///.*", "", symbol_map)

keep <- !is.na(symbol_map) & symbol_map != "" & symbol_map != "---"
expr_avg <- expr_avg[keep, , drop = FALSE]
symbol_map <- symbol_map[keep]

# Collapsing probes to genes by mean expression
expr_gene_sum <- rowsum(expr_avg, group = symbol_map, reorder = TRUE)
gene_counts <- table(symbol_map)
expr_gene <- expr_gene_sum / as.numeric(gene_counts[rownames(expr_gene_sum)])

# -----------------------------
# 4) Build paired tumor vs normal comparison
# -----------------------------
sample_names <- colnames(expr_gene)
group <- ifelse(grepl("^T", sample_names), "Tumor", "Normal")
patient_id <- sub("^[NT]", "", sample_names)

normal_names <- sample_names[group == "Normal"]
tumor_names  <- sample_names[group == "Tumor"]
common_patients <- intersect(sub("^N", "", normal_names), sub("^T", "", tumor_names))
common_patients <- sort(common_patients)

normal_cols <- match(paste0("N", common_patients), sample_names)
tumor_cols  <- match(paste0("T", common_patients), sample_names)

expr_normal <- expr_gene[, normal_cols, drop = FALSE]
expr_tumor  <- expr_gene[, tumor_cols, drop = FALSE]

diff_mat <- expr_tumor - expr_normal
n_pairs <- ncol(diff_mat)

# Vectorized paired t-test
mean_diff <- rowMeans(diff_mat)
sd_diff <- apply(diff_mat, 1, sd)
sd_diff[is.na(sd_diff) | sd_diff == 0] <- 1e-8

t_stat <- mean_diff / (sd_diff / sqrt(n_pairs))
p_value <- 2 * pt(q = -abs(t_stat), df = n_pairs - 1)
adj_p <- p.adjust(p_value, method = "BH")

results <- data.frame(
  gene = rownames(expr_gene),
  logFC = mean_diff,
  t_stat = t_stat,
  p_value = p_value,
  adj_p = adj_p,
  stringsAsFactors = FALSE
)
results <- results[order(results$adj_p, -abs(results$logFC)), ]
results$direction <- ifelse(results$logFC > 0, "Up in tumor", "Down in tumor")

write.csv(results, file = file.path(outdir, "all_gene_results.csv"), row.names = FALSE)

# -----------------------------
# 5) Drug Repurposing
# -----------------------------
# Example hand-built biology modules linked to common drugs.
# Higher score = the genes in that module are upregulated in tumor.
# It does NOT prove the drug will work! It is for hypothesis generation

drug_modules <- list(
  Losartan = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "POSTN", "SPARC", "THBS2", "ACTA2", "TGFB1", "CTGF", "CXCL12"),
  Aspirin = c("PTGS2", "IL6", "CXCL8", "TNF", "CCL2", "ICAM1", "NFKBIA"),
  Simvastatin = c("HMGCR", "HMGCS1", "MVK", "MVD", "FDPS", "FDFT1", "SQLE"),
  Metformin = c("SLC2A1", "HK2", "ENO1", "PGK1", "LDHA", "MTOR", "RPS6KB1"),
  Hydroxychloroquine = c("SQSTM1", "MAP1LC3B", "BECN1", "ATG5", "ATG7", "CTSB", "CTSD")
)

drug_notes <- data.frame(
  drug = c("Losartan", "Aspirin", "Simvastatin", "Metformin", "Hydroxychloroquine"),
  usual_use = c(
    "blood pressure / angiotensin receptor blocker",
    "pain / inflammation",
    "cholesterol lowering statin",
    "type 2 diabetes",
    "malaria / autoimmune disease"
  ),
  teaching_angle = c(
    "fibrosis / ECM / TGF-beta / desmoplasia",
    "inflammation / COX signaling",
    "mevalonate / cholesterol pathway",
    "glycolysis / growth signaling",
    "autophagy / lysosome"
  ),
  stringsAsFactors = FALSE
)

score_module <- function(genes, de_table) {
  genes_present <- intersect(genes, de_table$gene)
  if (length(genes_present) == 0) {
    return(data.frame(
      genes_found = 0,
      coverage = 0,
      raw_mean_logFC = NA_real_,
      module_score = NA_real_,
      best_gene = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  idx <- match(genes_present, de_table$gene)
  vals <- de_table$logFC[idx]
  coverage <- length(genes_present) / length(genes)
  raw_mean <- mean(vals, na.rm = TRUE)
  data.frame(
    genes_found = length(genes_present),
    coverage = coverage,
    raw_mean_logFC = raw_mean,
    module_score = raw_mean * coverage,
    best_gene = genes_present[which.max(vals)],
    stringsAsFactors = FALSE
  )
}

drug_scores <- do.call(
  rbind,
  lapply(names(drug_modules), function(d) {
    out <- score_module(drug_modules[[d]], results)
    data.frame(drug = d, out, stringsAsFactors = FALSE)
  })
)

drug_scores <- merge(drug_scores, drug_notes, by = "drug", all.x = TRUE, sort = FALSE)
drug_scores <- drug_scores[order(drug_scores$module_score, decreasing = TRUE), ]
write.csv(drug_scores, file = file.path(outdir, "drug_scoreboard.csv"), row.names = FALSE)

# -----------------------------
# 6) Visualization
# -----------------------------
known_label_genes <- c("COL1A1", "COL1A2", "POSTN", "THBS2", "SPARC", "S100P", "MSLN", "TGFBI", "CEACAM6", "KRT19")
label_idx <- which(results$gene %in% known_label_genes)
label_idx <- label_idx[order(results$adj_p[label_idx], -abs(results$logFC[label_idx]))]
label_idx <- head(label_idx, 10)

sig_idx <- which(results$adj_p < 0.01 & abs(results$logFC) > 1)

top_up <- head(results$gene[results$logFC > 0], 15)
top_down <- head(results$gene[results$logFC < 0], 15)
top_genes <- unique(c(top_up, top_down))


fibrosis_genes <- intersect(drug_modules$Losartan, rownames(expr_gene))
fibrosis_mat <- expr_gene[fibrosis_genes, , drop = FALSE]
fibrosis_mat_z <- t(scale(t(fibrosis_mat)))
fibrosis_mat_z[is.na(fibrosis_mat_z)] <- 0
fibrosis_score <- colMeans(fibrosis_mat_z, na.rm = TRUE)
normal_score <- fibrosis_score[normal_cols]
tumor_score <- fibrosis_score[tumor_cols]

make_volcano_plot <- function() {
  plot(
    results$logFC,
    -log10(pmax(results$adj_p, 1e-300)),
    pch = 20,
    col = "grey70",
    xlab = "Tumor - Normal average log2 difference",
    ylab = "-log10(FDR)",
    main = "Volcano plot: pancreatic tumor vs matched normal"
  )
  if (length(sig_idx) > 0) {
    points(
      results$logFC[sig_idx],
      -log10(pmax(results$adj_p[sig_idx], 1e-300)),
      pch = 20,
      col = "firebrick"
    )
  }
  abline(v = c(-1, 1), lty = 2, col = "grey40")
  abline(h = -log10(0.01), lty = 2, col = "grey40")
  if (length(label_idx) > 0) {
    text(
      results$logFC[label_idx],
      -log10(pmax(results$adj_p[label_idx], 1e-300)),
      labels = results$gene[label_idx],
      pos = 3,
      cex = 0.75
    )
  }
}

make_fibrosis_plot <- function() {
  ylim <- range(c(normal_score, tumor_score))
  plot(
    c(1, 2),
    ylim,
    type = "n",
    xaxt = "n",
    xlab = "",
    ylab = "Fibrosis / ECM score",
    main = "Losartan teaching signal: fibrosis score"
  )
  axis(1, at = c(1, 2), labels = c("Normal", "Tumor"))
  abline(h = 0, lty = 3, col = "grey50")
  for (i in seq_along(common_patients)) {
    lines(c(1, 2), c(normal_score[i], tumor_score[i]), col = adjustcolor("grey35", 0.4))
  }
  points(1 + runif(length(normal_score), -0.03, 0.03), normal_score, pch = 19, col = "steelblue")
  points(2 + runif(length(tumor_score), -0.03, 0.03), tumor_score, pch = 19, col = "tomato")
}

make_drug_score_plot <- function() {
  ds <- drug_scores
  ds <- ds[order(ds$module_score, decreasing = FALSE), ]
  labels <- paste0(ds$drug, " (n=", ds$genes_found, ")")
  par(mar = c(5, 10, 4, 2))
  mids <- barplot(
    ds$module_score,
    horiz = TRUE,
    las = 1,
    names.arg = labels,
    col = "grey70",
    xlab = "Toy repurposing module score\n(mean tumor-normal log2 difference across module genes)",
    main = "Drug ideas to discuss"
  )
  abline(v = 0, lty = 2, col = "grey40")
  text(
    x = ds$module_score,
    y = mids,
    labels = sprintf(" %.2f", ds$module_score),
    pos = ifelse(ds$module_score >= 0, 4, 2),
    cex = 0.8
  )
}

png(file.path(outdir, "01_volcano.png"), width = 1600, height = 1200, res = 180)
make_volcano_plot()
dev.off()

png(file.path(outdir, "02_fibrosis_score.png"), width = 1600, height = 1200, res = 180)
make_fibrosis_plot()
dev.off()

png(file.path(outdir, "03_drug_scoreboard.png"), width = 1600, height = 1200, res = 180)
make_drug_score_plot()
dev.off()