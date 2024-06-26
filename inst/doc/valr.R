## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## -----------------------------------------------------------------------------
library(valr)
library(dplyr)
library(ggplot2)
library(tibble)

## -----------------------------------------------------------------------------
library(valr)
library(dplyr)

snps <- read_bed(valr_example("hg19.snps147.chr22.bed.gz"))
genes <- read_bed(valr_example("genes.hg19.chr22.bed.gz"))

# find snps in intergenic regions
intergenic <- bed_subtract(snps, genes)
# distance from intergenic snps to nearest gene
nearby <- bed_closest(intergenic, genes)

nearby |>
  select(starts_with("name"), .overlap, .dist) |>
  filter(abs(.dist) < 1000)

## -----------------------------------------------------------------------------
bed_file <- valr_example("3fields.bed.gz")
read_bed(bed_file) # accepts filepaths or URLs

## -----------------------------------------------------------------------------
bed <- tribble(
  ~chrom, ~start,  ~end,
  "chr1", 1657492, 2657492,
  "chr2", 2501324, 3094650
)

bed

## -----------------------------------------------------------------------------
# a chromosome 100 basepairs in length
chrom <- tribble(
  ~chrom, ~start, ~end,
  "chr1", 0,      100
)

chrom

# single base-pair intervals
bases <- tribble(
  ~chrom, ~start, ~end,
  "chr1", 0,      1, # first base of chromosome
  "chr1", 1,      2, # second base of chromosome
  "chr1", 99,     100 # last base of chromosome
)

bases

## -----------------------------------------------------------------------------
#  # access the `refGene` tbl on the `hg38` assembly.
#  if (require(RMariaDB)) {
#    ucsc <- db_ucsc("hg38")
#    tbl(ucsc, "refGene")
#  }

## -----------------------------------------------------------------------------
x <- tribble(
  ~chrom, ~start, ~end,
  "chr1", 25,     50,
  "chr1", 100,    125
)

y <- tribble(
  ~chrom, ~start, ~end,
  "chr1", 30,     75
)

bed_glyph(bed_intersect(x, y))

## -----------------------------------------------------------------------------
x <- tribble(
  ~chrom, ~start, ~end,
  "chr1", 1, 50,
  "chr1", 10, 75,
  "chr1", 100, 120
)

bed_glyph(bed_merge(x))

## -----------------------------------------------------------------------------
x <- tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 1,      100,  "+",
  "chr1", 50,     150,  "+",
  "chr2", 100,    200,  "-"
)

y <- tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 50,     125,  "+",
  "chr1", 50,     150,  "-",
  "chr2", 50,     150,  "+"
)

# intersect tbls by strand
x <- group_by(x, strand)
y <- group_by(y, strand)

bed_intersect(x, y)

## -----------------------------------------------------------------------------
x <- group_by(x, strand)

y <- flip_strands(y)
y <- group_by(y, strand)

bed_intersect(x, y)

## -----------------------------------------------------------------------------
#  # calculate the mean and variance for a `value` column
#  bed_map(a, b, .mean = mean(value), .var = var(value))
#  
#  # report concatenated and max values for merged intervals
#  bed_merge(a, .concat = concat(value), .max = max(value))

## -----------------------------------------------------------------------------
# `valr_example()` identifies the path of example files
bedfile <- valr_example("genes.hg19.chr22.bed.gz")
genomefile <- valr_example("hg19.chrom.sizes.gz")
bgfile <- valr_example("hela.h3k4.chip.bg.gz")

genes <- read_bed(bedfile)
genome <- read_genome(genomefile)
y <- read_bedgraph(bgfile)

## -----------------------------------------------------------------------------
# generate 1 bp TSS intervals, `+` strand only
tss <- genes |>
  filter(strand == "+") |>
  mutate(end = start + 1)

# 1000 bp up and downstream
region_size <- 1000
# 50 bp windows
win_size <- 50

# add slop to the TSS, break into windows and add a group
x <- tss |>
  bed_slop(genome, both = region_size) |>
  bed_makewindows(win_size)

x

## -----------------------------------------------------------------------------
# map signals to TSS regions and calculate summary statistics.
res <- bed_map(x, y, win_sum = sum(value, na.rm = TRUE)) |>
  group_by(.win_id) |>
  summarize(
    win_mean = mean(win_sum, na.rm = TRUE),
    win_sd = sd(win_sum, na.rm = TRUE)
  )

res

## -----------------------------------------------------------------------------
x_labels <- seq(
  -region_size,
  region_size,
  by = win_size * 5
)

x_breaks <- seq(1, 41, by = 5)

sd_limits <- aes(
  ymax = win_mean + win_sd,
  ymin = win_mean - win_sd
)

ggplot(
  res,
  aes(
    x = .win_id,
    y = win_mean
  )
) +
  geom_point() +
  geom_pointrange(sd_limits) +
  scale_x_continuous(
    labels = x_labels,
    breaks = x_breaks
  ) +
  labs(
    x = "Position (bp from TSS)",
    y = "Signal",
    title = "Human H3K4me3 signal near transcription start sites"
  ) +
  theme_classic()

