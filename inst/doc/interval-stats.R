## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "img/interval-stats-",
  fig.height = 3,
  fig.align = "center",
  fig.width = 4
)

## ----load-data, message = FALSE, warning = FALSE-------------------------
library(valr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)

# load repeats and genes. Data in the valr package is restricted to chr22; the entire
# files can be downloaded from UCSC.
rpts <- read_bed(valr_example('hg19.rmsk.chr22.bed.gz'), n_fields = 6) 
gens <- read_bed(valr_example('hg19.refGene.chr22.bed.gz'), n_fields = 12)

# load chrom sizes
chrs <- read_genome(valr_example('hg19.chrom.sizes.gz'))

# create 1bp intervals representing transcription start sites
tss <- mutate(gens,
              .start = ifelse(strand == '+', start, end),
              .end = .start + 1) %>%
  select(chrom, start = .start, end = .end, name)

tss

## ----stats---------------------------------------------------------------
gen_stats <- function(x, y, genome, grp, type = NA) {
  group_by_(x, .dots = grp) %>%
    do(reldist = bed_reldist(., y, detail = TRUE) %>%
         select(.value = .reldist),
       absdist = bed_absdist(., y, genome) %>%
         select(.value = .absdist_scaled)
       ) %>%
    gather_('stat', 'value', setdiff(names(.), list(grp))) %>%
    mutate(type = type)
} 

## ----compute-------------------------------------------------------------
obs_stats <- gen_stats(rpts, tss, chrs, 'name', 'obs')

shfs <- bed_shuffle(rpts, chrs, within = TRUE)
shf_stats <- gen_stats(shfs, tss, chrs, 'name', 'shuf')

res <- bind_rows(obs_stats, shf_stats) %>%
  unnest(value) %>% 
  group_by(name, stat, type) %>%
  mutate(.id = row_number()) %>%
  spread(type, .value) %>%
  na.omit()

res

## ----pvalues, warning = FALSE--------------------------------------------
pvals <- res %>% do(twosided = broom::tidy(ks.test(.$obs, .$shuf)),
                    less = broom::tidy(ks.test(.$obs, .$shuf, alternative = 'less')),
                    greater = broom::tidy(ks.test(.$obs, .$shuf, alternative = 'greater'))) %>%
  gather(alt, type, -name, -stat) %>%
  unnest(type) %>%
  select(name:p.value) %>%
  arrange(p.value)

ggplot(pvals, aes(p.value)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(stat ~ alt) + theme_bw()

## ----qvalues-------------------------------------------------------------
pvals <- group_by(pvals, stat, alt) %>%
  mutate(q.value = p.adjust(p.value)) %>%
  ungroup() %>%
  arrange(q.value)

## ----ecfs----------------------------------------------------------------
res_fold <- res %>%
  gather(type, value, -name, -stat, -.id)

signif <- head(pvals, 25)
res_signif <- signif %>% left_join(res_fold, by = c('name','stat'))

ggplot(res_signif, aes(x = value, color = type)) +
  stat_ecdf() + 
  facet_wrap(name ~ stat) + theme_classic() + scale_x_log10()

## ----get_promoters-------------------------------------------------------
# Using the same data as before
gens
rpts


# create intervals 5kb upstream of tss representing promoters
promoters <- mutate(gens,
              .start = ifelse(strand == '+', start - 5000, end - 1),
              .end   = ifelse(strand == '+', start + 1, end + 5000),
              name   = ifelse(grepl("NR_", name), "non-coding", "coding")) %>%
  select(chrom, start = .start, end = .end, name, score, strand)
  
# select coding and non-coding promoters
nc_promoters <- filter(promoters, name == "non-coding")
coding_promoters <- filter(promoters, name == "coding")
nc_promoters
coding_promoters

## ----get_projections-----------------------------------------------------

# function to apply bed_projection to groups
gen_stats <- function(x, y, genome, grp, type = NA) {
  group_by_(x, .dots = grp) %>%
    do(repeat_counts = nrow(.),
       projection = bed_projection(., y, genome) 
       ) %>%
    mutate(type = type)
} 

pvals_nc <- gen_stats(rpts, nc_promoters, chrs, "name", "non_coding")
pvals_cd <- gen_stats(rpts, coding_promoters, chrs, "name", "coding")

pvals <- bind_rows(pvals_nc, pvals_cd) %>% 
  unnest() %>%
  select(-chrom)

#filter for repeat classes with at least 10 intervals
pvals <- filter(pvals, 
                repeat_counts > 10,
                obs_exp_ratio != 0)
# adjust pvalues 
pvals <- pvals %>%
  mutate(q.value = p.adjust(p.value))

pvals

## ------------------------------------------------------------------------

# show top 5 most significant repeats
sig <- pvals %>% 
  arrange(q.value) %>%
  group_by(type) %>%
  top_n(-5, q.value) %>%
  arrange(type)

knitr::kable(sig,
             caption = "The most significant repeats overlapping coding and non-coding gene promoters")

