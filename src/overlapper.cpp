// [[Rcpp::depends(RSeqAn)]]

#include <iostream>
#include <Rcpp.h>
#include <seqan/align.h>

using namespace Rcpp;
using namespace seqan;

// [[Rcpp::export]]

Rcpp::DataFrame overlapper(std::vector< std::string > seq1_set, std::vector< std::string > seq2_set)
{
  typedef seqan::String<char> TSequence;                 // sequence type
  typedef Align<TSequence,ArrayGaps> TAlign;      // align type
  typedef Row<TAlign>::Type TRow;                 // gapped sequence type

  // Perform match extension.
  Score<int, Simple> scoringScheme(1, -3, -3);

  int n = seq1_set.size(); // assuming seq1_set and seq2_set have the same number of sequences

  //initialization
  IntegerVector seq1_left_clipped_len_storage(n), seq1_right_clipped_len_storage(n), seq2_left_clipped_len_storage(n), seq2_right_clipped_len_storage(n);
  IntegerVector seq_1_len_storage(n), seq_2_len_storage(n);
  IntegerVector seq_1_ol_start_pos(n), seq_1_ol_end_pos(n), seq_2_ol_start_pos(n), seq_2_ol_end_pos(n);
  NumericVector pid(n);
  IntegerVector aln_ol_len(n), n_match(n);
  CharacterVector seq1_aln(n), seq2_aln(n);

  for(unsigned i = 0; i < n; ++i){
    std::string seq1= seq1_set[i];
    std::string seq2= seq2_set[i];

    // Perform a banded alignment.
    Align<CharString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    // Perform Alignment
    int score = globalAlignment(align, scoringScheme, AlignConfig<true, true, true, true>(), LinearGaps());
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);

    std::stringstream cur_seq1_aln;
    cur_seq1_aln << row1;
    seq1_aln[i] = cur_seq1_aln.str();

    std::stringstream cur_seq2_aln;
    cur_seq2_aln << row2;
    seq2_aln[i] = cur_seq2_aln.str();

    // Initialisation
    int aln_len = length(row1);
    int seq_1_len = seq1.length();
    int seq_2_len = seq2.length();
    int seq_1_left_clipped_len = 0;
    int seq_1_right_clipped_len = 0;
    int seq_2_left_clipped_len = 0;
    int seq_2_right_clipped_len = 0;

    // For sequence 1 clipped length
    for(unsigned i = 0; i < aln_len; ++i){
      if(toSourcePosition(row2, i)==0){
        seq_1_left_clipped_len++;
      }
    }

    for(unsigned i = 0; i < aln_len; ++i){
      if(toSourcePosition(row2, i)==seq_2_len){
        seq_1_right_clipped_len++;
      }
    }


    // For sequence 2 clipped length
    for(unsigned i = 0; i < aln_len; ++i){
      if(toSourcePosition(row1, i)==0){
        seq_2_left_clipped_len++;
      }
    }

    for(unsigned i = 0; i < aln_len; ++i){
      if(toSourcePosition(row1, i)==seq_1_len){
        seq_2_right_clipped_len++;
      }
    }

    seq1_left_clipped_len_storage[i] = seq_1_left_clipped_len - 1;
    seq1_right_clipped_len_storage[i] = seq_1_right_clipped_len;
    seq2_left_clipped_len_storage[i] = seq_2_left_clipped_len - 1;
    seq2_right_clipped_len_storage[i] = seq_2_right_clipped_len;
    seq_1_len_storage[i] = seq_1_len;
    seq_2_len_storage[i] = seq_2_len;
    seq_1_ol_start_pos[i] = toSourcePosition(row1, seq_1_left_clipped_len - 1); //0-based index
    seq_1_ol_end_pos[i] = seq_1_len - seq_1_right_clipped_len - 1;
    seq_2_ol_start_pos[i] = toSourcePosition(row2, seq_2_left_clipped_len - 1);
    seq_2_ol_end_pos[i] = seq_2_len - seq_2_right_clipped_len - 1;

    int total_left_clipped = (seq_1_left_clipped_len - 1) + (seq_2_left_clipped_len - 1);
    int total_right_clipped = seq_1_right_clipped_len + seq_2_right_clipped_len;
    setClippedBeginPosition(row(align,0), total_left_clipped);
    setClippedBeginPosition(row(align,1), total_left_clipped);
    setClippedEndPosition(row(align,0), aln_len - total_right_clipped);
    setClippedEndPosition(row(align,1), aln_len - total_right_clipped);

    aln_ol_len[i] = aln_len - total_right_clipped - total_left_clipped;

    // Getting alignment information
    AlignmentStats stats;
    computeAlignmentStats(stats, align, scoringScheme);
    // std::cout << "Clipping alignment\n" << align;

    pid[i] = stats.alignmentIdentity;
    n_match[i] = stats.numMatches;
  }

  return List::create(Named("seq1_left_clipped_len") = seq1_left_clipped_len_storage,
                           Named("seq1_right_clipped_len") = seq1_right_clipped_len_storage,
                           Named("seq2_left_clipped_len") = seq2_left_clipped_len_storage,
                           Named("seq2_right_clipped_len") = seq2_right_clipped_len_storage,
                           Named("seq1_len") = seq_1_len_storage,
                           Named("seq2_len") = seq_2_len_storage,
                           Named("seq1_ol_start") = seq_1_ol_start_pos,
                           Named("seq1_ol_end") = seq_1_ol_end_pos,
                           Named("seq2_ol_start") = seq_2_ol_start_pos,
                           Named("seq2_ol_end") = seq_2_ol_end_pos,
                           Named("seq1_aln") = seq1_aln,
                           Named("seq2_aln") = seq2_aln,
                           Named("pid") = pid,
                           Named("n_match") = n_match,
                           Named("aln_ol_len") = aln_ol_len);
}
