# ============================================================
# Corrected ALDEx2 differential abundance analysis
# Fix: collapse repeated per-animal timepoint samples to
# one animal-level count vector BEFORE running ALDEx2,
# since ALDEx2 has no random-effects/blocking option for
# >2-group or repeated-measures designs. Each animal must
# contribute exactly one column to the input table.
# ============================================================

library(ALDEx2)
library(data.table)
library(tidyverse)

# ---- 1. Load raw count table and metadata ----
# counts: rows = taxa, columns = sample_name (one column per sample,
# e.g. "332.T8.S"), values = raw read counts (NOT relative abundance/CLR)
counts <- fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv") %>%
  column_to_rownames(var = colnames(.)[1]) %>%
  as.matrix()

md <- fread("metagenomics/data_files/stool.metagenomic.metadata.with_amyloid.tsv")

# sanity check: each mouse_id should map to ~6 rows (one per ZT)
# and exactly one genotype_diet group
md %>% group_by(mouse_id) %>%
  summarise(n_samples = n(), n_groups = n_distinct(genotype_diet)) %>%
  print()
stopifnot(all((md %>% group_by(mouse_id) %>% summarise(n=n_distinct(genotype_diet)))$n == 1))

# ---- 2. Collapse counts to one column per animal ----
# Sum raw counts across timepoints per animal (NOT mean — ALDEx2 expects
# count data, and summing preserves library-size scaling properly;
# averaging counts directly is not appropriate pre-CLR).
shared_samples <- intersect(colnames(counts), md$sample_name)
length(shared_samples)
counts <- counts[, shared_samples]
md_match <- md[match(colnames(counts), md$sample_name), ]
stopifnot(all(colnames(counts) == md_match$sample_name))  # should pass now

per_animal_counts <- t(counts) %>%
  as.data.frame() %>%
  mutate(mouse_id = md_match$mouse_id) %>%
  group_by(mouse_id) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("mouse_id") %>%
  t() %>%
  as.matrix()
storage.mode(per_animal_counts) <- "integer"

# one-row-per-animal metadata, aligned to per_animal_counts columns
animal_md <- md %>%
  distinct(mouse_id, genotype_diet) %>%
  column_to_rownames("mouse_id")
animal_md <- animal_md[colnames(per_animal_counts), , drop = FALSE]
stopifnot(all(rownames(animal_md) == colnames(per_animal_counts)))

print(table(animal_md$genotype_diet))  # confirm n per group now ~7-9 animals, not 44-53

# ============================================================
# 3. Run ALDEx2 pairwise, matching your original comparisons
#    (Ntg.ALF vs Tg.ALF, Tg.ALF vs Tg.TRF, Ntg.ALF vs Tg.TRF)
# ============================================================

run_aldex_pairwise <- function(counts_mat, group_vec, group_a, group_b, mc.samples = 128) {
  keep <- group_vec %in% c(group_a, group_b)
  sub_counts <- as.matrix(counts_mat[, keep])
  storage.mode(sub_counts) <- "integer"
  sub_groups <- as.character(group_vec[keep])
  
  clr <- aldex.clr(sub_counts, conds = sub_groups, mc.samples = mc.samples, denom = "all", verbose = TRUE)
  tt  <- aldex.ttest(clr, paired.test = FALSE, verbose = TRUE)
  eff <- aldex.effect(clr, verbose = TRUE, include.sample.summary = FALSE)
  
  result <- merge(tt, eff, by = "row.names") %>%
    rename(taxon = Row.names) %>%
    mutate(comparison = paste(group_a, "vs", group_b))
  result
}

comparisons <- list(
  c("Ntg.ALF", "Tg.ALF"),
  c("Tg.ALF",  "Tg.TRF"),
  c("Ntg.ALF", "Tg.TRF")
)

aldex_results <- map(comparisons, function(cmp) {
  run_aldex_pairwise(per_animal_counts, animal_md$genotype_diet, cmp[1], cmp[2])
}) %>% bind_rows()

# ---- 4. Apply same effect size cutoff used in the manuscript ----
significant_taxa <- aldex_results %>%
  filter(abs(effect) >= 0.5)

write_csv(aldex_results, "metagenomics/aldex/all_comparisons.csv")
write_csv(significant_taxa, "metagenomics/aldex/significant_effect0.5.csv")

# also save per-comparison files matching original naming convention,
# so the downstream plotting notebook can read them in with minimal changes
for (cmp in comparisons) {
  label <- paste0(gsub("\\.", "", cmp[1]), "v", gsub("\\.", "", cmp[2]))
  out <- aldex_results %>%
    filter(comparison == paste(cmp[1], "vs", cmp[2])) %>%
    select(taxon, we.ep, we.eBH, wi.ep, wi.eBH, rab.all,
           paste0("rab.win.", cmp[1]),
           paste0("rab.win.", cmp[2]),
           diff.btw, diff.win, effect, overlap, comparison)
  #write_csv(out, paste0("metagenomics/aldex/", label, ".species.csv"))
  write_csv(out, paste0("metagenomics/aldex/all.","species", ".", label, ".csv"))
}
colnames(NtgALFvTgALF_species)
colnames(NtgALFvTgTRF_species)
colnames(TgALFvTgTRF_species)
