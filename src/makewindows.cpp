// makewindows.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

//[[Rcpp::export]]
DataFrame makewindows_impl(DataFrame df, int win_size = 0, int num_win = 0,
                           int step_size = 0, bool reverse = false) {

  NumericVector starts = df["start"] ;
  NumericVector ends = df["end"] ;

  std::vector<int> starts_out ;
  std::vector<int> ends_out ;
  std::vector<int> df_idxs ;
  std::vector<int> win_ids;

  for (int i = 0; i < starts.size(); ++i) {

    auto start = starts[i] ;
    auto end = ends[i] ;

    if (num_win > 0) {
      int ivl_len = end - start ;
      if (ivl_len < num_win) {
        continue ;
      }
      win_size = round((ivl_len) / num_win) ;
    }

    int by = win_size - step_size ;

    // create candidate starts
    std::vector<int> starts_by ;
    if (num_win > 0) {
      int win_idx = 0 ;
      for (int j = start; j < end && win_idx < num_win; j += win_size, ++win_idx) {
        starts_by.push_back(j) ;
      }
    } else {
      for (int j = start; j < end && j + win_size - by <= end; j += by) {
        starts_by.push_back(j) ;
      }
    }

    int nstarts = starts_by.size() ;
    for (int k = 0; k < nstarts; ++k) {

      auto start_by = starts_by[k] ;
      starts_out.push_back(start_by) ;

      // for last iteration through starts also extend to end
      if (start_by + win_size < end && k < nstarts - 1) {
        ends_out.push_back(start_by + win_size) ;
      } else {
        ends_out.push_back(end) ;
      }

      if (reverse) {
        win_ids.push_back(nstarts - k) ;
      } else {
        win_ids.push_back(k + 1) ;
      }

      df_idxs.push_back(i) ;
    }
  }

  DataFrame out = subset_dataframe(df, df_idxs) ;

  // replace original starts, ends, and .win_id
  out["start"] = starts_out ;
  out["end"] = ends_out ;
  out[".win_id"] = win_ids ;

  return out ;
}

/*** R
library(valr)
library(dplyr)

x <- trbl_interval(
  ~chrom, ~start, ~end,
  "chr1", 100,    200
)

bed_makewindows(x, win_size = 10)
bed_makewindows(x, win_size = 10, reverse = TRUE)
*/
