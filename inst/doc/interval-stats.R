## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "img/interval-stats-",
  fig.height = 4,
  fig.align = "center",
  fig.width = 4
)

## ----load-data, message = FALSE, warning = FALSE-------------------------
library(valr)
library(dplyr)
library(ggplot2)
library(tidyr)

# load repeats and genes. Data in the valr package is restricted to chr22; the entire
# files can be downloaded from UCSC.
rpts <- read_bed(valr_example('hg19.rmsk.chr22.bed.gz'), n_fields = 6) 
genes <- read_bed(valr_example('hg19.refGene.chr22.bed.gz'), n_fields = 12)

# load chrom sizes
genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))

# create 1bp intervals representing transcription start sites
tss <-
  mutate(genes,
         .start = ifelse(strand == '+', start, end),
         .end = .start + 1) %>%
  select(chrom, start = .start, end = .end, name)

tss

## ----stats---------------------------------------------------------------
distance_stats <- function(x, y, genome, grp, type = NA) {
  group_by_(x, .dots = grp) %>%
    do(
       reldist = bed_reldist(., y, detail = TRUE) %>%
         select(.value = .reldist),
       absdist = bed_absdist(., y, genome) %>%
         select(.value = .absdist)
       ) %>%
    tidyr::gather_('stat', 'value', setdiff(names(.), list(grp))) %>%
    mutate(type = type)
} 

## ----compute_obs---------------------------------------------------------
obs_stats <- distance_stats(rpts, tss, genome, 'name', 'obs')
obs_stats

## ----compute_shf---------------------------------------------------------
seed <- 1010486
shfs <- bed_shuffle(rpts, genome, within = TRUE, seed = seed)
shf_stats <- distance_stats(shfs, tss, genome, 'name', 'shuf')

## ----bind_res------------------------------------------------------------
res <- bind_rows(obs_stats, shf_stats) %>%
  tidyr::unnest(value) %>% 
  group_by(name, stat, type) %>%
  mutate(.id = row_number()) %>%
  tidyr::spread(type, .value) %>%
  na.omit()

res

## ----pvalues-------------------------------------------------------------
library(broom)

# silence the tests with ties ...
ks.test.quiet <- function(...) {suppressWarnings(ks.test(...))}

pvals <- res %>% do(twosided = broom::tidy(ks.test.quiet(.$obs, .$shuf)),
                    less = broom::tidy(ks.test.quiet(.$obs, .$shuf, alternative = 'less')),
                    greater = broom::tidy(ks.test.quiet(.$obs, .$shuf, alternative = 'greater'))) %>%
  tidyr::gather(alt, type, -name, -stat) %>%
  unnest(type) %>%
  select(name:p.value) %>%
  arrange(p.value)

## ----pvalue_viz----------------------------------------------------------
ggplot(pvals, aes(p.value)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(stat ~ alt) + theme_bw()

## ----qvalues-------------------------------------------------------------
pvals <-
  group_by(pvals, stat, alt) %>%
  mutate(q.value = p.adjust(p.value)) %>%
  ungroup() %>%
  arrange(q.value)

## ----ecfs----------------------------------------------------------------
res_gather <- tidyr::gather(res, type, value, -name, -stat, -.id)

signif <- head(pvals, 5)

res_signif <-
  signif %>%
  left_join(res_gather, by = c('name', 'stat'))

ggplot(res_signif, aes(x = value, color = type)) +
  stat_ecdf() + 
  facet_grid(stat ~ name) + theme_classic() + scale_x_log10()

## ----get_promoters-------------------------------------------------------
# create intervals 5kb upstream of tss representing promoters
promoters <-
  bed_flank(genes, genome, left = 5000, strand = TRUE) %>%
  mutate(name = ifelse(grepl('NR_', name), 'non-coding', 'coding')) %>%
  select(chrom:strand)
  
# select coding and non-coding promoters
promoters_coding <- filter(promoters, name == 'coding')
promoters_ncoding <- filter(promoters, name == 'non-coding')

promoters_coding

promoters_ncoding

## ----get_projections-----------------------------------------------------

# function to apply bed_projection to groups
projection_stats <- function(x, y, genome, grp, type = NA) {
  group_by_(x, .dots = grp) %>%
    do(n_repeats = nrow(.),
       projection = bed_projection(., y, genome)) %>%
    mutate(type = type)
} 

pvals_coding <- projection_stats(rpts, promoters_coding, genome, 'name', 'coding')
pvals_ncoding <- projection_stats(rpts, promoters_ncoding, genome, 'name', 'non_coding')

pvals <-
  bind_rows(pvals_ncoding, pvals_coding) %>% 
  tidyr::unnest() %>%
  select(-chrom)

# filter for repeat classes with at least 10 intervals
pvals <- filter(pvals, n_repeats > 10,
                obs_exp_ratio != 0)

# adjust pvalues 
pvals <- mutate(pvals, q.value = p.adjust(p.value))

pvals

## ----table---------------------------------------------------------------
library(knitr)

# find and show top 5 most significant repeats
signif_tests <-
  pvals %>% 
  arrange(q.value) %>%
  group_by(type) %>%
  top_n(-5, q.value) %>%
  arrange(type)

knitr::kable(signif_tests)

