#!/usr/bin/env Rscript
# =============================================================================
# circos_methylation.R
#
# Circular genome visualisation of methylation, features and (optionally) motifs
# for prokaryotic genomes, built directly from the raw inputs:
#   - reference FASTA
#   - GFF3 annotation
#   - methylation calls (modkit-style bedMethyl)
#   - optional newline-delimited motif list (IUPAC allowed)
#
# Ring layout (outermost -> innermost):
#   [motif rings]   one density ring per motif, OUTSIDE the replicons
#                   (only when --motifs is supplied)
#   replicon        ideogram with oriC (green) and ter (red) marked + Mb scale
#   features        + strand (red, up) and - strand (blue, down) on ONE ring,
#                   the minus strand mirrored ("head down")
#   methylation     density (area) with LOESS smoothing + highlight lines at the
#                   most highly methylated genes
#   GC skew         diverging (G>C green / C>G purple); confirms oriC/ter
#
# oriC/ter are detected from cumulative GC-skew (origin = global minimum;
# terminus = maximum of the skew re-started at the origin, which avoids the
# arbitrary-contig-start artefact), or supplied via --oric/--ter.
# =============================================================================

suppressWarnings(suppressMessages({
  # ---- package bootstrap -----------------------------------------------------
  .need_cran <- c("optparse", "data.table", "circlize")
  .need_bioc <- c("Biostrings")
  .missing_cran <- .need_cran[!vapply(.need_cran, requireNamespace, logical(1),
                                      quietly = TRUE)]
  if (length(.missing_cran)) {
    install.packages(.missing_cran, repos = "https://cloud.r-project.org")
  }
  .missing_bioc <- .need_bioc[!vapply(.need_bioc, requireNamespace, logical(1),
                                      quietly = TRUE)]
  if (length(.missing_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(.missing_bioc, update = FALSE, ask = FALSE)
  }
  library(optparse)
  library(data.table)
  library(circlize)
  library(Biostrings)
}))

# ---------------------------------------------------------------------------- #
# CLI
# ---------------------------------------------------------------------------- #
opt_list <- list(
  make_option("--gff", type = "character", help = "GFF3 annotation [required]"),
  make_option("--fasta", type = "character", help = "Reference FASTA [required]"),
  make_option("--methylation", type = "character",
              help = "bedMethyl methylation file [required]"),
  make_option("--motifs", type = "character", default = NULL,
              help = "Optional motif list (one IUPAC motif per line)."),
  make_option("--mod-codes", type = "character", default = NULL,
              dest = "mod_codes",
              help = "Comma-separated modification codes to keep, e.g. 'a'."),
  make_option("--feature-types", type = "character", default = "CDS",
              dest = "feature_types",
              help = "Comma-separated GFF types to keep [default %default]."),
  make_option("--min-coverage", type = "double", default = 0, dest = "min_cov",
              help = "Minimum valid coverage for a call [default %default]."),
  make_option("--min-frac", type = "double", default = 0.5, dest = "min_frac",
              help = "Fraction-modified threshold for 'methylated' [default %default]."),
  make_option("--bin", type = "double", default = 10000,
              help = "Bin size (bp) for the density rings [default %default]."),
  make_option("--skew-window", type = "double", default = 10000,
              dest = "skew_window",
              help = "Window (bp) for oriC/ter GC-skew detection [default %default]."),
  make_option("--loess-span", type = "double", default = 0.30, dest = "loess_span",
              help = "LOESS span for methylation smoothing [default %default]."),
  make_option("--highlight-top", type = "integer", default = 12,
              dest = "highlight_top",
              help = "Number of most-methylated genes to highlight [default %default]."),
  make_option("--highlight-min-calls", type = "integer", default = 4,
              dest = "highlight_min_calls",
              help = "Min calls in a gene to be eligible for highlight [default %default]."),
  make_option("--oric", type = "character", default = NULL,
              help = "Override oriC, 'seqid:pos[,seqid:pos]' (0-based)."),
  make_option("--ter", type = "character", default = NULL,
              help = "Override ter, 'seqid:pos[,seqid:pos]' (0-based)."),
  make_option("--strip-version", action = "store_true", default = FALSE,
              dest = "strip_version",
              help = "Strip trailing .N from seqids (auto if it reconciles GFF/FASTA)."),
  make_option("--out", type = "character", default = "circos_genome.pdf",
              help = "Output file (.pdf/.png/.svg) [default %default]."),
  make_option("--width", type = "double", default = 11, help = "Figure width (in)."),
  make_option("--height", type = "double", default = 11, help = "Figure height (in)."),
  make_option("--title", type = "character", default = NULL, help = "Plot title.")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (req in c("gff", "fasta", "methylation"))
  if (is.null(opt[[req]])) stop("Missing required --", req)

split_csv <- function(x) if (is.null(x)) NULL else strsplit(x, ",")[[1]]
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

parse_loci <- function(s) {
  out <- list()
  if (is.null(s)) return(out)
  for (tok in strsplit(s, ",")[[1]]) {
    kv <- strsplit(trimws(tok), ":")[[1]]
    if (length(kv) == 2) out[[kv[1]]] <- as.numeric(kv[2])
  }
  out
}

# ---------------------------------------------------------------------------- #
# Read inputs
# ---------------------------------------------------------------------------- #
message("Reading FASTA / GFF / methylation ...")
genome <- readDNAStringSet(opt$fasta)
names(genome) <- sub("\\s.*", "", names(genome))

read_gff <- function(path, ftypes) {
  dt <- fread(path, sep = "\t", header = FALSE, fill = TRUE, quote = "",
              showProgress = FALSE)
  dt <- dt[!startsWith(as.character(V1), "#")]
  dt <- dt[!is.na(V3) & !is.na(V4) & !is.na(V5)]
  if (!is.null(ftypes)) dt <- dt[V3 %in% ftypes]
  dt[, `:=`(seqid = as.character(V1),
            start = as.integer(V4) - 1L,
            end   = as.integer(V5) - 1L,
            strand = as.character(V7),
            type = as.character(V3))]
  dt[, locus := ifelse(grepl("locus_tag=", V9),
                       sub(".*locus_tag=([^;]+).*", "\\1", V9), NA_character_)]
  dt[, .(seqid, start, end, strand, type, locus)]
}
feat <- read_gff(opt$gff, split_csv(opt$feature_types))

read_meth <- function(path, keep_codes, min_cov) {
  dt <- fread(path, sep = "\t", header = FALSE, showProgress = FALSE)
  dt <- dt[, .(seqid = as.character(V1), pos = as.integer(V2),
               code = as.character(V4), strand = as.character(V6),
               cov = as.numeric(V10), frac = as.numeric(V11))]
  dt[frac > 1, frac := frac / 100]               # bedMethyl percent -> fraction
  if (!is.null(keep_codes)) dt <- dt[code %in% keep_codes]
  dt[is.na(cov) | cov >= min_cov]
}
meth <- read_meth(opt$methylation, split_csv(opt$mod_codes), opt$min_cov)

# ---- seqid reconciliation (auto strip a pure .N version mismatch) ----
gff_ids <- unique(feat$seqid); fa_ids <- names(genome)
do_strip <- opt$strip_version
if (!do_strip && length(intersect(gff_ids, fa_ids)) == 0 &&
    length(intersect(strip_ver(gff_ids), strip_ver(fa_ids))) > 0) {
  message("[auto] stripping '.N' version suffix to reconcile seqids.")
  do_strip <- TRUE
}
if (do_strip) {
  names(genome) <- strip_ver(names(genome))
  feat[, seqid := strip_ver(seqid)]
  meth[, seqid := strip_ver(seqid)]
}
if (length(intersect(names(genome), unique(feat$seqid))) == 0)
  stop("GFF and FASTA seqids do not match; try --strip-version.")
glen <- setNames(as.integer(width(genome)), names(genome))   # named lengths

# replicon labels from FASTA descriptions
hdr <- names(readDNAStringSet(opt$fasta))
labels <- setNames(sub("\\s.*", "", hdr), sub("\\s.*", "", hdr))
for (i in seq_along(hdr)) {
  sid <- sub("\\s.*", "", hdr[i]); if (do_strip) sid <- strip_ver(sid)
  d <- tolower(hdr[i])
  lab <- sid
  if (grepl("chromosome", d)) lab <- "chromosome"
  else if (grepl("plasmid", d)) {
    tok <- sub(".*[Pp]lasmid +([^, ]+).*", "\\1", hdr[i])
    lab <- if (nzchar(tok) && tok != hdr[i]) tok else "plasmid"
  }
  labels[sid] <- lab
}

# ---------------------------------------------------------------------------- #
# oriC / ter from cumulative GC skew
# ---------------------------------------------------------------------------- #
gc_oriter <- function(seq, window) {
  L <- length(seq)
  starts <- seq(1, L, by = window)
  widths <- pmin(window, L - starts + 1)
  v <- Views(seq, start = starts, width = widths)
  lf <- letterFrequency(v, c("G", "C"))
  delta <- lf[, "G"] - lf[, "C"]
  cum <- cumsum(delta); n <- length(delta)
  oi <- which.min(cum)
  rot <- cumsum(delta[c(oi:n, seq_len(oi - 1))])
  ti <- ((oi - 1 + which.max(rot) - 1) %% n) + 1
  centers <- (starts - 1) + widths / 2
  list(oric = centers[oi], ter = centers[ti])
}
oric_ov <- parse_loci(opt$oric); ter_ov <- parse_loci(opt$ter)
orit <- list()
for (sid in names(genome)) {
  d <- gc_oriter(genome[[sid]], opt$skew_window)
  oc <- if (!is.null(oric_ov[[sid]])) oric_ov[[sid]] else d$oric
  tr <- if (!is.null(ter_ov[[sid]])) ter_ov[[sid]] else d$ter
  orit[[sid]] <- list(oric = oc, ter = tr)
}

# ---------------------------------------------------------------------------- #
# Binning helpers
# ---------------------------------------------------------------------------- #
bin <- opt$bin
bin_counts <- function(pos, L) {
  nb <- max(1, ceiling(L / bin))
  if (length(pos) == 0) return(rep(0L, nb))
  idx <- pmin(floor(pos / bin) + 1L, nb)
  tabulate(idx, nbins = nb)
}
bin_centers <- function(L) {
  nb <- max(1, ceiling(L / bin))
  (seq_len(nb) - 0.5) * bin
}
skew_bins <- function(seq) {
  L <- length(seq); nb <- max(1, ceiling(L / bin))
  starts <- (seq_len(nb) - 1) * bin + 1
  widths <- pmin(bin, L - starts + 1)
  lf <- letterFrequency(Views(seq, start = starts, width = widths), c("G", "C"))
  g <- lf[, "G"]; c <- lf[, "C"]
  (g - c) / pmax(g + c, 1)
}

# motif occurrences (both strands), IUPAC-aware
find_motif_pos <- function(seq, motif) {
  pat <- DNAString(motif)
  fwd <- start(matchPattern(pat, seq, fixed = "subject"))
  rev <- start(matchPattern(reverseComplement(pat), seq, fixed = "subject"))
  as.integer(c(fwd, rev)) - 1L
}
motifs <- if (!is.null(opt$motifs))
  toupper(trimws(readLines(opt$motifs))) else character(0)
motifs <- motifs[nzchar(motifs) & !startsWith(motifs, "#")]

# ---------------------------------------------------------------------------- #
# Per-sector data + global normalisers
# ---------------------------------------------------------------------------- #
ord <- names(glen)[order(-glen)]
sectors <- ord
feat_dat <- meth_dat <- skew_dat <- list()
motif_dat <- setNames(vector("list", length(motifs)), motifs)
meth_meth <- meth[is.na(frac) | frac >= opt$min_frac]

for (sid in ord) {
  L <- glen[[sid]]
  cen <- bin_centers(L)
  fp <- feat[seqid == sid & strand == "+"]
  fm <- feat[seqid == sid & strand == "-"]
  feat_dat[[sid]] <- data.frame(x = cen,
                                pos = bin_counts(fp$start + (fp$end - fp$start) / 2, L),
                                neg = bin_counts(fm$start + (fm$end - fm$start) / 2, L))
  mm <- meth_meth[seqid == sid]
  meth_dat[[sid]] <- data.frame(x = cen, dens = bin_counts(mm$pos, L))
  skew_dat[[sid]] <- data.frame(x = cen, skew = skew_bins(genome[[sid]]))
  for (mo in motifs) {
    p <- find_motif_pos(genome[[sid]], mo)
    motif_dat[[mo]][[sid]] <- data.frame(x = cen, dens = bin_counts(p, L))
  }
}
feat_max <- max(1, unlist(lapply(feat_dat, function(d) max(d$pos, d$neg))))
meth_max <- max(1, unlist(lapply(meth_dat, function(d) max(d$dens))))
skew_max <- max(1e-9, unlist(lapply(skew_dat, function(d) max(abs(d$skew)))))
motif_max <- sapply(motifs, function(mo)
  max(1, unlist(lapply(motif_dat[[mo]], function(d) max(d$dens)))))

# normalise
for (sid in ord) {
  feat_dat[[sid]]$pos <- feat_dat[[sid]]$pos / feat_max
  feat_dat[[sid]]$neg <- feat_dat[[sid]]$neg / feat_max
  meth_dat[[sid]]$dens <- meth_dat[[sid]]$dens / meth_max
  skew_dat[[sid]]$skew <- skew_dat[[sid]]$skew / skew_max
  for (mo in motifs)
    motif_dat[[mo]][[sid]]$dens <- motif_dat[[mo]][[sid]]$dens / motif_max[[mo]]
}

# ---- highly methylated genes (for the highlight lines) ----
genes <- feat[strand %in% c("+", "-")]
genes[, gid := ifelse(is.na(locus), paste0(seqid, ":", start), locus)]
setkey(genes, seqid, start, end)
mc <- meth[, .(seqid, s = pos, e = pos, frac)]
ov <- foverlaps(mc, genes, by.x = c("seqid", "s", "e"),
                by.y = c("seqid", "start", "end"), nomatch = NULL)
agg <- ov[, .(mean_frac = mean(frac, na.rm = TRUE), n = .N),
          by = .(seqid, gid, start, end)]
agg <- agg[n >= opt$highlight_min_calls]
top <- if (nrow(agg)) {
  agg[order(-mean_frac)][seq_len(min(opt$highlight_top, .N))]
} else agg
highlight <- list()
for (sid in ord) {
  t <- top[seqid == sid]
  highlight[[sid]] <- if (nrow(t))
    data.frame(mid = (t$start + t$end) / 2, gid = t$gid, mf = t$mean_frac) else NULL
}

# ---------------------------------------------------------------------------- #
# Plot
# ---------------------------------------------------------------------------- #
message("Rendering ", opt$out, " ...")
ext <- tolower(tools::file_ext(opt$out))
if (ext == "pdf") {
  pdf(opt$out, width = opt$width, height = opt$height)
} else if (ext == "svg") {
  svg(opt$out, width = opt$width, height = opt$height)
} else {
  png(opt$out, width = opt$width, height = opt$height, units = "in", res = 200)
}

pal <- c("#2c6fbb", "#e08a1e", "#3a9b78", "#b5446e", "#7d6bb0")
repcol <- setNames(pal[(seq_along(ord) - 1) %% length(pal) + 1], ord)

xlim <- cbind(rep(0, length(ord)), as.numeric(glen[ord]))
par(mar = c(1, 1, 2, 1))
gaps <- c(rep(3, length(ord) - 1), 14)             # big gap to host the legend
circos.clear()
circos.par(start.degree = 90 - 7, gap.after = gaps,
           cell.padding = c(0, 0, 0, 0), track.margin = c(0.004, 0.004),
           points.overflow.warning = FALSE)
circos.initialize(factors = factor(ord, levels = ord), xlim = xlim)

# ---- (outer) motif rings -----------------------------------------------------
motif_pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
for (i in seq_along(motifs)) {
  mo <- motifs[i]; col <- motif_pal[(i - 1) %% length(motif_pal) + 1]
  local({
    m <- mo; cc <- col
    circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
      panel.fun = function(x, y) {
        d <- motif_dat[[m]][[CELL_META$sector.index]]
        circos.lines(d$x, d$dens, area = TRUE,
                     col = adjustcolor(cc, 0.55), border = cc, lwd = 0.5)
        if (CELL_META$sector.numeric.index == length(ord))
          circos.text(CELL_META$xlim[2], 0.5, paste0(" ", m), col = cc,
                      cex = 0.5, adj = c(0, 0.5), facing = "downward")
      })
  })
}

# ---- replicon ideogram + oriC/ter -------------------------------------------
circos.track(ylim = c(0, 1), track.height = 0.09, bg.border = NA,
  panel.fun = function(x, y) {
    si <- CELL_META$sector.index; xl <- CELL_META$xlim
    circos.rect(xl[1], 0, xl[2], 0.55, col = repcol[si], border = "white", lwd = 0.5)
    circos.text(mean(xl), 0.93, labels[si], facing = "bending.inside",
                niceFacing = TRUE, cex = 0.8, font = 2, col = repcol[si])
    # Mb axis
    mj <- seq(0, xl[2], by = 5e5)
    circos.axis(h = 0.55, major.at = mj, labels = sprintf("%.1f", mj / 1e6),
                labels.cex = 0.35, major.tick.length = mm_y(1),
                labels.facing = "clockwise", col = "grey40", labels.col = "grey30")
    o <- orit[[si]]$oric; t <- orit[[si]]$ter
    circos.segments(o, 0, o, 0.75, col = "#1e8449", lwd = 2.5)
    circos.text(o, 0.9, "ori", col = "#1e8449", cex = 0.5, font = 2,
                facing = "clockwise", niceFacing = TRUE)
    circos.segments(t, 0, t, 0.75, col = "#c0392b", lwd = 2.5)
    circos.text(t, 0.9, "ter", col = "#c0392b", cex = 0.5, font = 2,
                facing = "clockwise", niceFacing = TRUE)
  })

# ---- features: + (red, up) / - (blue, down) on one mirrored ring -------------
circos.track(ylim = c(-1, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    d <- feat_dat[[CELL_META$sector.index]]
    hw <- bin / 2
    circos.rect(d$x - hw, 0, d$x + hw, d$pos, col = "#d62728", border = NA)
    circos.rect(d$x - hw, -d$neg, d$x + hw, 0, col = "#1f77b4", border = NA)
    circos.segments(CELL_META$xlim[1], 0, CELL_META$xlim[2], 0,
                    col = "grey55", lwd = 0.4)
  })

# ---- methylation density: area + LOESS + highlight lines ---------------------
circos.track(ylim = c(0, 1.18), track.height = 0.17, bg.border = NA,
  panel.fun = function(x, y) {
    si <- CELL_META$sector.index; d <- meth_dat[[si]]
    circos.lines(d$x, pmin(d$dens, 1), area = TRUE,
                 col = "#f5b7b1", border = "#cd6155", lwd = 0.4)
    if (nrow(d) >= 6) {
      lo <- tryCatch(loess(dens ~ x, data = d, span = opt$loess_span),
                     error = function(e) NULL)
      if (!is.null(lo)) {
        xs <- seq(min(d$x), max(d$x), length.out = 250)
        pr <- pmax(0, predict(lo, newdata = data.frame(x = xs)))
        circos.lines(xs, pmin(pr, 1.05), col = "#7b241c", lwd = 2)
      }
    }
    h <- highlight[[si]]
    if (!is.null(h) && nrow(h) > 0) {
      circos.segments(h$mid, 0, h$mid, 1.05, col = "#e67e22", lwd = 1.4)
      circos.points(h$mid, rep(1.05, nrow(h)), pch = 18, col = "#e67e22", cex = 0.6)
      circos.text(h$mid, 1.12, h$gid, col = "#b9690f", cex = 0.32,
                  facing = "clockwise", niceFacing = TRUE)
    }
  })

# ---- GC skew (innermost) -----------------------------------------------------
circos.track(ylim = c(-1, 1), track.height = 0.11, bg.border = NA,
  panel.fun = function(x, y) {
    d <- skew_dat[[CELL_META$sector.index]]; hw <- bin / 2
    circos.rect(d$x - hw, 0, d$x + hw, pmax(d$skew, 0), col = "#27ae60", border = NA)
    circos.rect(d$x - hw, pmin(d$skew, 0), d$x + hw, 0, col = "#8e44ad", border = NA)
    circos.segments(CELL_META$xlim[1], 0, CELL_META$xlim[2], 0,
                    col = "grey55", lwd = 0.4)
  })

circos.clear()

# ---- title + legend ----------------------------------------------------------
ttl <- if (!is.null(opt$title)) opt$title else
  paste0("Genome-wide methylation & feature distribution",
         if (length(motifs)) "  (+ motif rings)" else "")
title(main = ttl, cex.main = 1.05, line = 0)

leg_items <- c("genes + strand", "genes - strand",
               "methylation density", "LOESS", "highly methylated gene",
               "GC skew G>C", "GC skew C>G", "oriC", "ter")
leg_cols  <- c("#d62728", "#1f77b4", "#f5b7b1", "#7b241c", "#e67e22",
               "#27ae60", "#8e44ad", "#1e8449", "#c0392b")
leg_pch   <- c(15, 15, 15, NA, 18, 15, 15, NA, NA)
leg_lty   <- c(NA, NA, NA, 1, NA, NA, NA, 1, 1)
if (length(motifs)) {
  leg_items <- c(leg_items, paste("motif", motifs))
  mc2 <- motif_pal[(seq_along(motifs) - 1) %% length(motif_pal) + 1]
  leg_cols <- c(leg_cols, mc2)
  leg_pch <- c(leg_pch, rep(15, length(motifs)))
  leg_lty <- c(leg_lty, rep(NA, length(motifs)))
}
par(xpd = NA)
legend("bottomleft", legend = leg_items, col = leg_cols, pch = leg_pch,
       lty = leg_lty, bty = "n", cex = 0.62, pt.cex = 1.1,
       inset = c(0.0, 0.0), ncol = 1, seg.len = 1.4, text.col = "grey15")

invisible(dev.off())
message("Done: ", opt$out)