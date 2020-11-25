// [[Rcpp::depends(RSeqAn)]]
// [[Rcpp::plugins(cpp14)]]
#include <iostream>
#include <Rcpp.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace Rcpp;
using namespace seqan;

// [[Rcpp::export]]

Rcpp::StringVector t_coffee(std::vector<std::string> strings)
{
  int n = strings.size();
  Align<DnaString> align;
  resize(rows(align), n);
  
  for (int i = 0; i < n; ++i){
    assignSource(row(align, i), strings[i]);
  }
  globalMsaAlignment(align, SimpleScore(2, -2, 0, -5));
  int aln_len = length(row(align, 0));
  Rcpp::StringVector out(n);
  for (unsigned i = 0; i < n; ++i){
    std::stringstream curout;
    curout << seqan::row(align, i);
    out[i] = curout.str();
  }
  return out;
}