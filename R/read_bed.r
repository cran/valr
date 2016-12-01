#' @title Read BED and related files.
#'   
#' @description read functions for BED and related formats. Filenames can be 
#'   local file or URLs. The read functions load data into tbls with consistent
#'   \code{chrom}, \code{start} and \code{end} colnames.
#'   
#' @param filename file or URL
#' @param n_fields number fields in the BED file
#' @param col_types column type spec for \code{readr::read_tsv}
#' @param sort sort the tbl by chrom and start
#' @param ... options to pass to \code{readr::read_tsv}
#'   
#' @return \code{data_frame}
#'   
#' @family read-funcs
#'   
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'   
#' @examples
#' # read_bed assumes 3 field BED format. 
#' read_bed(valr_example('3fields.bed.gz'))
#' 
#' read_bed(valr_example('6fields.bed.gz'), n_fields = 6)
#' 
#' # result is sorted by chrom and start unless `sort = FALSE`
#' read_bed(valr_example('3fields.bed.gz'), sort = FALSE)
#' 
#' @export
read_bed <- function(filename, n_fields = 3, col_types = bed12_coltypes,
                     sort = TRUE, ...) {
  
  coltypes <- col_types[1:n_fields]
  colnames <- names(coltypes)
  
  bed_raw <- readr::read_tsv(filename, col_names = colnames, col_types = coltypes, ...)
  bed_tbl <- tibble::as_tibble(bed_raw) 

  if (sort) {
    bed_tbl <- bed_sort(bed_tbl)
  }
 
  bed_tbl
}

#' @rdname read_bed
#' 
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#' 
#' @examples
#' 
#' read_bed12(valr_example('mm9.bed12.gz'))
#' 
#' @export
read_bed12 <- function(filename, ...) {
  bed12_tbl <- read_bed(filename, n_fields = 12)
  bed12_tbl
}

#' @rdname read_bed
#' 
#' @details \url{https://genome.ucsc.edu/goldenPath/help/bedgraph.html}
#' 
#' @examples
#' 
#' read_bedgraph(valr_example('test.bg.gz'))
#' 
#' @export
read_bedgraph <- function(filename, ...) {
  # load as bed4, rename `value` column and covert to double
  bedgraph_tbl <- read_bed(filename, n_fields = 4, sort = FALSE) %>%
    rename(value = name) %>%
    mutate(value = as.double(value))
  bedgraph_tbl <- bed_sort(bedgraph_tbl)
  bedgraph_tbl
}

#' @rdname read_bed
#' 
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format12}
#' 
#' @examples
#' 
#' read_narrowpeak(valr_example('sample.narrowPeak.gz'))
#' 
#' @export
read_narrowpeak <- function(filename, ...) {
  colnames <- names(peak_coltypes)
  x <- readr::read_tsv(filename, col_types = peak_coltypes, col_names = colnames)
  x
}

#' @rdname read_bed
#' 
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format13}
#' 
#' @examples
#' 
#' read_broadpeak(valr_example('sample.broadPeak.gz'))
#' 
#' @export
read_broadpeak <- function(filename, ...) {
  coltypes <- peak_coltypes[1:length(peak_coltypes)-1]
  colnames <- names(coltypes)
  x <- readr::read_tsv(filename, col_names = colnames, col_types = coltypes)
  x
}

peak_coltypes <- list(
  chrom = readr::col_character(),
  start = readr::col_integer(),
  end = readr::col_integer(),
  name = readr::col_character(),
  score = readr::col_integer(),
  strand = readr::col_character(),
  signal = readr::col_double(),
  pvalue = readr::col_double(),
  qvalue = readr::col_double(),
  peak = readr::col_integer()
)

bed12_coltypes <- list(
  chrom = readr::col_character(),
  start = readr::col_integer(),
  end = readr::col_integer(),
  name = readr::col_character(),
  score = readr::col_character(),
  strand = readr::col_character(),
  cds_start = readr::col_integer(),
  cds_end = readr::col_integer(),
  item_rgb = readr::col_character(),
  exon_count = readr::col_integer(),
  exon_sizes = readr::col_character(),
  exon_starts = readr::col_character()
)
