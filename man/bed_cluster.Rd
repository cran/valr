% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_cluster.r
\name{bed_cluster}
\alias{bed_cluster}
\title{Cluster neighboring intervals.}
\usage{
bed_cluster(x, max_dist = 0)
}
\arguments{
\item{x}{\link{ivl_df}}

\item{max_dist}{maximum distance between clustered intervals.}
}
\value{
\link{ivl_df} with \code{.id} column specifying sets of clustered intervals.
}
\description{
The output \code{.id} column can be used in downstream grouping operations. Default
\code{max_dist = 0} means that both overlapping and book-ended intervals will be
clustered.
}
\details{
input tbls are grouped by \code{chrom} by default, and additional
groups can be added using \code{\link[dplyr:group_by]{dplyr::group_by()}}. For example,
grouping by \code{strand} will constrain analyses to the same strand. To
compare opposing strands across two tbls, strands on the \code{y} tbl can
first be inverted using \code{\link[=flip_strands]{flip_strands()}}.
}
\examples{
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,    200,
  "chr1", 180,    250,
  "chr1", 250,    500,
  "chr1", 501,    1000,
  "chr2", 1,      100,
  "chr2", 150,    200
)

bed_cluster(x)

# glyph illustrating clustering of overlapping and book-ended intervals
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1,      10,
  "chr1", 5,      20,
  "chr1", 30,     40,
  "chr1", 40,     50,
  "chr1", 80,     90
)

bed_glyph(bed_cluster(x), label = ".id")

}
\seealso{
\url{https://bedtools.readthedocs.io/en/latest/content/tools/cluster.html}

Other single set operations: 
\code{\link{bed_complement}()},
\code{\link{bed_flank}()},
\code{\link{bed_merge}()},
\code{\link{bed_partition}()},
\code{\link{bed_shift}()},
\code{\link{bed_slop}()}
}
\concept{single set operations}
