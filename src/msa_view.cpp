// [[Rcpp::depends(RSeqAn)]]

#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> cxx_query_to_reference_position(std::string raln, std::string qaln){
  int r_pos = 0;
  int aln_len = raln.size();
  std::vector<int> q_in_rpos(aln_len);
  for(unsigned i = 0; i < aln_len; ++i){
    q_in_rpos[i] = r_pos;
    if(raln[i] != '-'){
      r_pos++;
    }
  }
  return(q_in_rpos);
}

std::vector<int> cxx_ref_extension(Rcpp::List vecs){
  int list_size = vecs.size();
  Rcpp::IntegerVector collected_last(list_size);
  int current_max;
  current_max = 0;
  for(unsigned i = 0; i < list_size; ++i){
    SEXP cur_list = vecs[i];
    Rcpp::IntegerVector tmp_vec(cur_list);
    std::vector<int> cur_vec = as<std::vector<int>> (tmp_vec);
    int cur_last = cur_vec.back();
    if(cur_last > current_max){
      current_max = cur_last;
    };
  }

  Rcpp::IntegerVector collected(current_max + 1);

  for(unsigned i = 0; i < list_size; ++i){
    SEXP cur_list = vecs[i];
    Rcpp::IntegerVector tmp_vec(cur_list);
    std::vector<int> cur_vec = as<std::vector<int>> (tmp_vec);
    int counter = 0;
    int previous_pos = 0;
    for(unsigned q = 0; q < cur_vec.size(); ++q){
      int cur_pos = cur_vec[q];
      if(cur_pos == previous_pos){
        counter++;
      }else{
        if(counter > collected[previous_pos]){
          collected[previous_pos] = counter;
        }
        counter = 1;
      }
      previous_pos = cur_pos;
    }
    collected[previous_pos] = counter;
  }

  //int total_collected = std::accumulate(collected.begin(), collected.end(), 0);
  //std::cout << total_collected;
  //for(unsigned i = 0; i < collected.size(); ++i){
  //  std::cout << collected[i];
  //}
  std::vector<int> combined;
  int counter = 0L;
  for(unsigned j = 0; j < collected.size(); ++j){
    for(unsigned p = 0; p < collected[j]; ++p){
      //combined[counter] <- j;
      //++counter;
      //std::cout << j;
      combined.push_back(j);
    }
  }
  return(combined);
}

// [[Rcpp::export]]
Rcpp::StringVector msa_view(std::vector<std::string> raln, std::vector<std::string> qaln){
  int n_aln = raln.size();
  Rcpp::List r_translated(n_aln);
  for(unsigned i = 0; i < n_aln; ++i){
    std::string r_cur = raln[i];
    r_translated[i] = cxx_query_to_reference_position(r_cur, r_cur);
  }
  std::vector<int> extended_ref_pos = cxx_ref_extension(r_translated);
  Rcpp:StringVector collected_extended_q(n_aln);

  int extended_len = extended_ref_pos.size();
  std::string space = "-";

  for(unsigned a = 0; a < n_aln; ++a){
    std::string q_aln = qaln[a];
    std::vector<int> q_pos = cxx_query_to_reference_position(raln[a], q_aln);
    int q_counter = 0;
    std::string extended_q (extended_len, 'X');
    int cur_r_pos = 0;
    int cur_q_pos = 0;
    int n_q = q_pos.size();
    int next_r_pos = 0;
    for(unsigned i = 0; i < (extended_len - 1); ++i){
      //std::cout << "When i = " << i << "\n";
      if(q_counter != (n_q - 1)){
        int next_q_pos = q_pos[q_counter + 1];
        int next_r_pos = extended_ref_pos[i + 1];
        //std::cout << "next_q_pos = " << next_q_pos << "\n";
        //std::cout << "next_r_pos = " << next_r_pos << "\n";
        //std::cout << "cur_q_pos = " << cur_q_pos << "\n";
        //std::cout << "cur_r_pos = " << cur_r_pos << "\n";

        if(cur_r_pos == cur_q_pos & next_r_pos == next_q_pos){
          extended_q[i] = q_aln[q_counter];
          extended_q[i + 1] = q_aln[q_counter + 1];
          cur_q_pos = q_pos[q_counter + 1]; // this actually can be replaced by next_q_pos, but it does not work
          q_counter++;
        }else{
          extended_q[i] = space[0];
        }
      }else{
        extended_q[i + 1] = space[0];
      }
      cur_r_pos = extended_ref_pos[i + 1]; // this actually can be replaced by next_q_pos, but it does not work
      //std::cout << extended_q << "\n";
    }
    collected_extended_q[a] = extended_q;
  }
  return(collected_extended_q);
}

