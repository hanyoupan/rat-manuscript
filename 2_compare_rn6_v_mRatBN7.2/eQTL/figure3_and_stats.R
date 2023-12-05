library(tidyverse)
library(patchwork)

eqtls <- read_tsv("data/eqtl_moved.txt.gz", col_types = "ccdcll") %>%
    mutate(moved = esnp_move | tss_move)

reloc <- read_tsv("data/relocated.txt", col_types = "ccdcicicllcicic")

t_to_c <- reloc |>
    filter(rn6_type == "trans",
           rn7_type == "cis")
t_to_c |>
    distinct(gene_id)
# Number of trans to cis associations:
nrow(t_to_c)
# Number of original trans associations:
with(eqtls, sum(rn6_type == "trans"))
# Fraction of rn6 trans relocating to cis:
nrow(t_to_c) / with(eqtls, sum(rn6_type == "trans"))

c_to_t <- reloc |>
    filter(rn6_type == "cis",
           rn7_type == "trans")
c_to_t |>
    with(n_distinct(gene_id))
# Number of cis to trans associations:
nrow(c_to_t)
# Number of original cis associations:
with(eqtls, sum(rn6_type == "cis"))
# Fraction of rn6 cis relocating to trans:
nrow(c_to_t) / with(eqtls, sum(rn6_type == "cis"))

##############################################
## Examples to visualize with genome viewer ##
##############################################

t_to_c_ex <- t_to_c |>
    group_by(gene_id, rn6_esnp_chrom, rn7_esnp_chrom) |>
    slice_min(pval, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(gene_id, rn6_esnp_chrom, rn6_esnp_pos)

# Line start x-coord:
# 1008.3 + (1277.3 - 1008.3) * ((81309880 - 81200000) / 200000)
# Line end x-coord:
# 1240 + (1498 - 1240) * ((45871675 - 45600000) / 200000)

############
## Figure ##
############

n_trans_egenes <- eqtls |>
    filter(rn6_type == "trans") |>
    with(n_distinct(gene_id))
p1 <- eqtls |>
    filter(rn6_type == "trans") |>
    left_join(select(reloc, variant_id, gene_id, rn7_type),
              by = c("variant_id", "gene_id")) |>
    mutate(
        gene_id = fct_reorder(gene_id, pval, .fun = min),
        eSNPs = case_when(
            moved & rn6_type == "trans" & rn7_type == "cis" ~ "Relocated, resulting in cis-eQTL",
            moved & rn6_type == "trans" & rn7_type == "intermediate" ~ "Relocated, new TSS distance 1-5 Mb",
            moved & rn6_type == "trans" & rn7_type == "trans" ~ "Relocated, still trans-eQTL",
            !moved ~ "Did not relocate"
        ) |>
            fct_relevel("Relocated, resulting in cis-eQTL",
                        "Relocated, new TSS distance 1-5 Mb",
                        "Relocated, still trans-eQTL",
                        "Did not relocate")
    ) |>
    ggplot(aes(x = gene_id, fill = eSNPs)) +
    geom_bar(color = "black", lwd = 0.1) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_fill_manual(values = c("red", "#dd8888", "black", "gray")) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.59, 0.74),
        legend.justification = c(0, 0),
        legend.key.size = unit(11, "pt"),
    ) +
    ylab("trans-eQTL SNPs (Rnor6.0)") +
    xlab(NULL) +
    labs(fill = "SNP and/or TSS:")
p2 <- eqtls |>
    filter(rn6_type == "trans") |>
    group_by(gene_id) |>
    summarise(min_p = min(pval), .groups = "drop") |>
    mutate(gene_id = fct_reorder(gene_id, min_p, .fun = min)) |>
    ggplot(aes(x = gene_id, y = 0, fill = -log10(min_p))) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        legend.position = c(0.59, 31),
        legend.justification = c(0, 1),
        legend.key.width = unit(11, "pt"),
        legend.key.height = unit(12, "pt"),
    ) +
    xlab(str_glue("{n_trans_egenes} genes with Rnor6.0 trans-eQTLs")) +
    ylab(NULL) +
    labs(fill = expression(paste("Top -", log[10](P))))

p1 / plot_spacer() / p2 + plot_layout(heights = c(40, -2.5, 1))

# ggsave("trans_to_cis_plot/eQTL_figure.png", width = 7, height = 4.5, dpi = 600)
ggsave("trans_to_cis_plot/eQTL_figure.pdf", width = 7, height = 4.5)
