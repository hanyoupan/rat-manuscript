library(BSgenome.Rnorvegicus.UCSC.rn6)
library(GenomicRanges)
library(tidyverse)

#' Return whether each location falls within any relocated segment
moved <- function(chrom, pos, segments) {
    seg_gr <- with(segments, GRanges(tname, IRanges(tstart, tend)))
    loci <- GRanges(chrom, IRanges(pos, pos))
    countOverlaps(loci, seg_gr) > 0
}

#' Return index of the segment containing each location
i_seg_overlap <- function(chrom, pos, segments) {
    ""
    seg_gr <- with(segments, GRanges(tname, IRanges(tstart, tend)))
    loci <- GRanges(chrom, IRanges(pos, pos))
    findOverlaps(loci, seg_gr, select = "first")
}

#' Return the estimated new position after relocation
#' 
#' Note: the alignments aren't perfect, so these are approximate locations.
new_pos <- function(pos, i_seg, segments) {
    if (is.na(i_seg)) return(NA)
    hit <- slice(segments, i_seg)
    if (hit$strand == "+") {
        hit$qstart + (pos - hit$tstart)
    } else {
        hit$qstart + (hit$tend - pos)
    }
}

segments <- pafr::read_paf("data/rn6_mRatBN72_chr_translocation_lt50bp.paf.gz") |>
    as_tibble() |>
    filter(mapq >= 60)
    ## Currently all are mapq == 60:
    # arrange(desc(mapq)) # So first GRanges match is highest quality

genes <- read_tsv("data/genes.txt", col_types = "ccc---i-----") |>
    select(gene_id,
           rn6_tss_chrom = chrom,
           rn6_tss_pos = tss)

# Use NAcc eQTLs from https://doi.org/10.1093/nar/gkac912
eqtls1 <- read_tsv("data/NQCT.trans_qtl_pairs.txt.gz",
                   col_types = "ccd---") |>
    filter(pval < 1e-8)
eqtls2 <- eqtls1 |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("rn6_esnp_chrom", "rn6_esnp_pos"), sep = ":", convert = TRUE,
             remove = FALSE) |>
    mutate(rn6_esnp_chrom = str_replace(rn6_esnp_chrom, "chr", "")) |>
    relocate(rn6_esnp_chrom, rn6_esnp_pos, .after = pval) |>
    left_join(genes, by = "gene_id") |>
    mutate(rn6_type = case_when(
        rn6_esnp_chrom == rn6_tss_chrom & abs(rn6_esnp_pos - rn6_tss_pos) < 1e6 ~ "cis",
        rn6_esnp_chrom == rn6_tss_chrom & abs(rn6_esnp_pos - rn6_tss_pos) < 5e6 ~ "intermediate",
        TRUE ~ "trans"
    ))
eqtls <- eqtls2 |>
    mutate(esnp_move = moved(rn6_esnp_chrom, rn6_esnp_pos, segments),
           tss_move = moved(rn6_tss_chrom, rn6_tss_pos, segments),
           moved = esnp_move | tss_move)

reloc <- eqtls |>
    filter(moved) |>
    select(-moved) |>
    mutate(i_seg_esnp = i_seg_overlap(rn6_esnp_chrom, rn6_esnp_pos, segments),
           i_seg_tss = i_seg_overlap(rn6_tss_chrom, rn6_tss_pos, segments)) |>
    rowwise() |>
    mutate(rn7_esnp_chrom = if (esnp_move) segments$qname[i_seg_esnp] else rn6_esnp_chrom,
           rn7_esnp_pos = if (esnp_move) new_pos(rn6_esnp_pos, i_seg_esnp, segments) else rn6_esnp_pos,
           rn7_tss_chrom = if (tss_move) segments$qname[i_seg_tss] else rn6_tss_chrom,
           rn7_tss_pos = if (tss_move) new_pos(rn6_tss_pos, i_seg_tss, segments) else rn6_tss_pos) |>
    ungroup() |>
    mutate(rn7_type = case_when(
        rn7_esnp_chrom == rn7_tss_chrom & abs(rn7_esnp_pos - rn7_tss_pos) < 1e6 ~ "cis",
        rn7_esnp_chrom == rn7_tss_chrom & abs(rn7_esnp_pos - rn7_tss_pos) < 5e6 ~ "intermediate",
        TRUE ~ "trans"
    ))

eqtls |>
    select(-(rn6_esnp_chrom:rn6_tss_pos),
           -moved) |>
    write_tsv("data/eqtl_moved.txt.gz")

reloc |>
    select(-i_seg_esnp, -i_seg_tss) |>
    write_tsv("data/relocated.txt")

#############################################
## Extract sequences of relocated segments ##
#############################################
# Requested by collaborator for analysis

cis_to_trans <- reloc |>
    filter(rn6_type == "cis",
           rn7_type == "trans")
cis_to_trans |>
    count(gene_id)

i_seg_cis_to_trans <- cis_to_trans |>
    distinct(i_seg_esnp) |>
    pull()

seg_c_t <- segments |>
    slice(i_seg_cis_to_trans) |>
    arrange(as.integer(tname), tstart)

seqs <- seg_c_t |>
    mutate(tname = str_c("chr", tname)) |>
    with(getSeq(Rnorvegicus, tname, tstart, tend)) %>%
    setNames(with(seg_c_t, str_glue("chr{tname}:{tstart}-{tend}")))

writeXStringSet(seqs, "cis_to_trans_5Mb.fa")

# Subset to TSS distance < 1Mb instead of 5Mb:

i_seg_cis_to_trans_1mb <- cis_to_trans |>
    filter(abs(rn6_esnp_pos - rn6_tss_pos) < 1e6) |>
    distinct(i_seg_esnp) |>
    pull()

seg_c_t_1mb <- segments |>
    slice(i_seg_cis_to_trans_1mb) |>
    arrange(as.integer(tname), tstart)

seqs_1mb <- seg_c_t_1mb |>
    mutate(tname = str_c("chr", tname)) |>
    with(getSeq(Rnorvegicus, tname, tstart, tend)) %>%
    setNames(with(seg_c_t_1mb, str_glue("chr{tname}:{tstart}-{tend}")))

writeXStringSet(seqs_1mb, "cis_to_trans_1Mb.fa")
