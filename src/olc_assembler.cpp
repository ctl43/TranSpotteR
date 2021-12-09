// [[Rcpp::depends(RSeqAn)]]
// [[Rcpp::plugins(cpp14)]]

#include <iostream>
#include <Rcpp.h>
#include <seqan/align.h>

using namespace Rcpp;
using namespace seqan;
using namespace std;

// [[Rcpp::export]]
CharacterVector rcpp_split(string &s, char delim)
{
  CharacterVector result;
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim))
  {
    result.push_back(item);
  }
  return result;
}

// [[Rcpp::export]]
List string_split(CharacterVector x, string sep){
  int size = x.size();
  List result(size);
  std::string cur;
  char delim = sep[0];
  for (int i = 0; i < size; i++){
    cur = x[i];
    CharacterVector splittedString = rcpp_split(cur, delim);
    result[i] = splittedString;
  }
  return result;
}


// [[Rcpp::export]]
CharacterVector combine_two(CharacterVector vec1, CharacterVector vec2){
  size_t vec1_len =  vec1.size();
  size_t vec2_len =  vec2.size();
  size_t new_len = vec1_len + vec2_len;
  CharacterVector new_vec(new_len);
  for(size_t j = 0; j < vec1_len; ++j){
    new_vec[j] = vec1[j];
  }
  for(size_t i = 0; i < vec2_len; ++i){
    new_vec[i + vec1_len] = vec2[i];
  }
  return(new_vec);
}

// https://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder
template <int RTYPE>
IntegerVector order_impl(const Vector<RTYPE>& x, bool desc) {
  auto n = x.size();
  IntegerVector idx = no_init(n), out_idx = no_init(n);
  std::iota(idx.begin(), idx.end(), static_cast<size_t>(1));
  if (desc) {
    auto comparator = [&x](size_t a, size_t b){ return x[a - 1] > x[b - 1]; };
    std::stable_sort(idx.begin(), idx.end(), comparator);
  } else {
    auto comparator = [&x](size_t a, size_t b){ return x[a - 1] < x[b - 1]; };
    std::stable_sort(idx.begin(), idx.end(), comparator);
    // simulate na.last
    size_t nas = 0;
    for (size_t i = 0; i < n; ++i, ++nas)
      if (!Vector<RTYPE>::is_na(x[idx[i] - 1])) break;
      std::rotate(idx.begin(), idx.begin() + nas, idx.end());
  }
  for(int j = 0; j<n; ++j){
    out_idx[j] = idx[j] - 1;
  }
  return out_idx;
}

// [[Rcpp::export]]
IntegerVector decreasing_order(SEXP x, bool desc = true) {
  switch(TYPEOF(x)) {
  case INTSXP: return order_impl<INTSXP>(x, desc);
  case REALSXP: return order_impl<REALSXP>(x, desc);
  case STRSXP: return order_impl<STRSXP>(x, desc);
  default: stop("Unsupported type.");
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_which(LogicalVector x) {
  if(x.size()==0){
    return IntegerVector();
  }
  Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
  return v[x];
}


// [[Rcpp::export]]
SEXP subset_df(SEXP x,
               IntegerVector row_indices,
               IntegerVector column_indices) {
  int column_indices_n = column_indices.size();
  int row_indices_n = row_indices.size();
  row_indices = row_indices;
  column_indices = column_indices;
  List output = no_init(column_indices_n);
  CharacterVector x_names = as<CharacterVector>(Rf_getAttrib(x, R_NamesSymbol));
  output.attr("names") = x_names[column_indices];
  
  // Loop and fill!
  for (int j = 0; j < column_indices_n; ++j)
  {
    SEXP element = VECTOR_ELT(x, column_indices[j]);
    SEXP vec = PROTECT(
      Rf_allocVector(TYPEOF(element), row_indices_n)
    );
    for (int i = 0; i < row_indices_n; ++i)
    {
      switch (TYPEOF(vec))
      {
      case REALSXP:
        REAL(vec)[i] =
          REAL(element)[row_indices[i]];
        break;
      case INTSXP:
      case LGLSXP:
        INTEGER(vec)[i] =
          INTEGER(element)[row_indices[i]];
        break;
      case STRSXP:
        SET_STRING_ELT(vec, i, STRING_ELT(element, row_indices[i]));
        break;
      }
    }
    UNPROTECT(1);
    SET_VECTOR_ELT(output, j, vec);
  }
  return output;
}

// [[Rcpp::export]]
List rcpp_rbind(List x, List y) {
  int n = x.size();
  List output = List(n);
  CharacterVector touch_x = x[0];
  CharacterVector touch_y = y[0];
  int nx = touch_x.size();
  int ny  = touch_y.size();
  int new_size = nx + ny;
  CharacterVector list_name = x.names();
  
  for(int i = 0; i<n; ++i){
    SEXP cur_x = x[i];
    SEXP cur_y = y[i];
    
    SEXP vec = Rf_allocVector(TYPEOF(cur_x), new_size);
    int counter = 0;
    for(int j = 0; j<nx; ++j){
      switch (TYPEOF(vec))
      {
      case REALSXP:
        REAL(vec)[counter] = REAL(cur_x)[j];
        break;
      case INTSXP:
      case LGLSXP:
        INTEGER(vec)[counter] = INTEGER(cur_x)[j];
        break;
      case STRSXP:
        SET_STRING_ELT(vec, counter, STRING_ELT(cur_x, j));
        break;
      }
      ++counter;
    }
    
    
    for(int k = 0; k<ny; ++k){
      switch (TYPEOF(vec))
      {
      case REALSXP:
        REAL(vec)[counter] = REAL(cur_y)[k];
        break;
      case INTSXP:
      case LGLSXP:
        INTEGER(vec)[counter] = INTEGER(cur_y)[k];
        break;
      case STRSXP:
        SET_STRING_ELT(vec, counter, STRING_ELT(cur_y, k));
        break;
      }
      ++counter;
    }
    output[i] = vec;
  }
  output.names() = list_name;
  return(output);
}

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

// [[Rcpp::export]]
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
    
    if(counter > collected[previous_pos]){ // for the last element
      collected[previous_pos] = counter;
    }
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

// [[Rcpp::export]]
Rcpp::StringVector replcae_space(CharacterVector x) {
  size_t sz=x.size();
  Rcpp::StringVector out(sz);
  for(unsigned p=0; p<sz;++p){
    std::string current = Rcpp::as<std::string>(x[p]);
    for (unsigned i=0; i<current.length(); ++i){
      if(current[i]=='-'){
        current[i]=' ';
      }else{
        break;
      }
    }
    
    for (unsigned j=current.length()-1; j>=0; j--){
      if(current[j]=='-'){
        current[j]=' ';
      }else{
        break;
      }
    }
    out[p]=current;
  }
  return out;
}

// [[Rcpp::export]]
SEXP rcpp_overlapper(CharacterVector seq_set1,
                CharacterVector seq_set2 = "",
                bool full_set = true, bool flip_side = true)
{
  typedef seqan::String<char> TSequence;                 // sequence type
  typedef Align<TSequence,ArrayGaps> TAlign;      // align type
  typedef Row<TAlign>::Type TRow;                 // gapped sequence type
  
  // Perform match extension.
  Score<int, Simple> scoringScheme(1, -3, -3);
  int z, n, c;
  IntegerVector set1, set2;
  
  // Generate all combination
  if(full_set){
    seq_set2 = seq_set1;
    n = seq_set1.size();
    z = ((n*n)-n)/2; //nCr
    set1 = IntegerVector(z);
    set2 = IntegerVector(z);
    c = 0;
    for(int a = 0; a < (n - 1); ++a){
      for(int b = (a + 1); b < n; ++b){
        set1[c] = a;
        set2[c] = b;
        ++c;
      }
    }
  }else{
    z = seq_set1.size() * seq_set2.size();
    set1 = IntegerVector(z);
    set2 = IntegerVector(z);
    int m = 0;
    for(int k = 0; k < seq_set1.size(); ++k){
      for(int w = 0; w < seq_set2.size(); ++w){
        set1[m] = k;
        set2[m] = w;
        ++m;
      }
    }
  }

  //initialization
  IntegerVector seq1_left_clipped_len_storage(z), seq1_right_clipped_len_storage(z), seq2_left_clipped_len_storage(z), seq2_right_clipped_len_storage(z);
  IntegerVector seq_1_len_storage(z), seq_2_len_storage(z);
  IntegerVector seq_1_ol_start_pos(z), seq_1_ol_end_pos(z), seq_2_ol_start_pos(z), seq_2_ol_end_pos(z);
  NumericVector pid(z), aln_ol_len(z);
  CharacterVector seq1_aln(z), seq2_aln(z), seq1_name(z), seq2_name(z), cur_seq1_name(1), cur_seq2_name(1);
  std::string seq1, seq2;
  
  CharacterVector seqset_name1 = seq_set1.names();
  CharacterVector seqset_name2 = seq_set2.names();
  
  int counter = 0, e, j;
  for(unsigned p = 0; p < set1.size(); ++p){
    e = set1[p];
    cur_seq1_name[0] = seqset_name1[e];
    seq1 = Rcpp::as<std::string>(seq_set1[cur_seq1_name]);
    
    j = set2[p];
    cur_seq2_name[0] = seqset_name2[j];
    seq2 = Rcpp::as<std::string>(seq_set2[cur_seq2_name]);
    
    //Rcout << counter << "\t" << cur_seq1_name << "\t" << cur_seq2_name << "\n";
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
    
    std::stringstream cur_seq2_aln;
    cur_seq2_aln << row2;
    
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
    //Rcout << seq_1_right_clipped_len << "\n";
    if(seq_1_right_clipped_len != 0 & flip_side){
      seq2_name[counter] = cur_seq1_name[0];
      seq1_name[counter] = cur_seq2_name[0];
      seq_2_len_storage[counter] = seq_1_len;
      seq_1_len_storage[counter] = seq_2_len;
      seq2_left_clipped_len_storage[counter] = seq_1_left_clipped_len - 1;
      seq2_right_clipped_len_storage[counter] = seq_1_right_clipped_len;
      seq1_left_clipped_len_storage[counter] = seq_2_left_clipped_len - 1;
      seq1_right_clipped_len_storage[counter] = seq_2_right_clipped_len;
      seq_2_ol_start_pos[counter] = toSourcePosition(row1, seq_1_left_clipped_len - 1); //0-based index
      seq_1_ol_start_pos[counter] = toSourcePosition(row2, seq_2_left_clipped_len - 1);
      seq2_aln[counter] = cur_seq1_aln.str();
      seq1_aln[counter] = cur_seq2_aln.str();
    }else{
      seq1_name[counter] = cur_seq1_name[0];
      seq2_name[counter] = cur_seq2_name[0];
      seq_1_len_storage[counter] = seq_1_len;
      seq_2_len_storage[counter] = seq_2_len;
      seq1_left_clipped_len_storage[counter] = seq_1_left_clipped_len - 1;
      seq1_right_clipped_len_storage[counter] = seq_1_right_clipped_len;
      seq2_left_clipped_len_storage[counter] = seq_2_left_clipped_len - 1;
      seq2_right_clipped_len_storage[counter] = seq_2_right_clipped_len;
      seq_1_ol_start_pos[counter] = toSourcePosition(row1, seq_1_left_clipped_len - 1); //0-based index
      seq_2_ol_start_pos[counter] = toSourcePosition(row2, seq_2_left_clipped_len - 1);
      seq1_aln[counter] = cur_seq1_aln.str();
      seq2_aln[counter] = cur_seq2_aln.str();
    }
    
    int total_left_clipped = (seq_1_left_clipped_len - 1) + (seq_2_left_clipped_len - 1);
    int total_right_clipped = seq_1_right_clipped_len + seq_2_right_clipped_len;
    setClippedBeginPosition(row(align,0), total_left_clipped);
    setClippedBeginPosition(row(align,1), total_left_clipped);
    setClippedEndPosition(row(align,0), aln_len - total_right_clipped);
    setClippedEndPosition(row(align,1), aln_len - total_right_clipped);
    
    // Getting alignment information
    AlignmentStats stats;
    computeAlignmentStats(stats, align, scoringScheme);
    int cur_ol_len = aln_len - total_right_clipped - total_left_clipped;
    aln_ol_len[counter] = cur_ol_len;
    if(cur_ol_len > 0){
      pid[counter] = stats.alignmentIdentity;      
    }else{
      pid[counter] = 0;
    }
    ++counter;
  }
  
  List df = List::create(
    Named("seq1") = seq1_name,
    Named("seq2") = seq2_name,
    Named("seq1_left_clipped_len") = seq1_left_clipped_len_storage,
    Named("seq1_right_clipped_len") = seq1_right_clipped_len_storage,
    Named("seq2_left_clipped_len") = seq2_left_clipped_len_storage,
    Named("seq2_right_clipped_len") = seq2_right_clipped_len_storage,
    Named("seq1_ol_start") = seq_1_ol_start_pos,
    Named("seq2_ol_start") = seq_2_ol_start_pos,
    Named("seq1_aln") = seq1_aln,
    Named("seq2_aln") = seq2_aln,
    Named("seq1_len") = seq_1_len_storage,
    Named("seq2_len") = seq_2_len_storage,
    Named("pid") = pid,
    Named("aln_ol_len") = aln_ol_len);
  return df;
}

// [[Rcpp::export]]
SEXP filter_and_order(List df, double min_pid, int min_ol_len){
  NumericVector pid = df["pid"];
  NumericVector aln_ol_len = df["aln_ol_len"];
  IntegerVector seq1_ol_start = df["seq1_ol_start"];
  IntegerVector seq2_ol_start = df["seq2_ol_start"];
  LogicalVector pass = (pid >= min_pid) & (aln_ol_len >= min_ol_len) & (seq1_ol_start == 0 | seq2_ol_start == 0);
  IntegerVector row_idx = rcpp_which(pass);
  IntegerVector col_idx = Rcpp::seq(0, df.size()-1);
  df = subset_df(df, row_idx, col_idx);
  aln_ol_len = df["aln_ol_len"];
  IntegerVector new_order = decreasing_order(aln_ol_len);
  df = subset_df(df, new_order, col_idx);
  return(df);
}

// [[Rcpp::export]]
SEXP which_is_fully_covered(List info){
  IntegerVector seq1_left_clip = info["seq1_left_clipped_len"];
  IntegerVector seq1_right_clip = info["seq1_left_clipped_len"];
  IntegerVector seq2_left_clip = info["seq2_left_clipped_len"];
  IntegerVector seq2_right_clip = info["seq2_right_clipped_len"];
  LogicalVector fully_covered_1 =  seq1_left_clip == 0 & seq1_right_clip == 0;
  LogicalVector fully_covered_2 = seq2_left_clip == 0 & seq2_right_clip == 0;
  fully_covered_1[fully_covered_1 & fully_covered_2] = false; // head to tail fully overlap
  LogicalVector has_sth = any(fully_covered_1);
  CharacterVector fully_covered;
  CharacterVector seq1 = info["seq1"];
  CharacterVector seq2 = info["seq2"];
  fully_covered = combine_two(seq1[fully_covered_1], seq2[fully_covered_2]);
  return(fully_covered);
}


std::string single_get_consensus(StringVector dna){
  CharacterVector nucleotide{"A", "T", "C", "G", "-"};
  Rcpp::String touch_sv = dna[0];
  std::string touch = touch_sv.get_cstring();
  Rcpp::String cur_seq;
  int m = touch.length();
  int n = dna.size();
  IntegerVector na(m), nt(m), nc(m), ng(m), nh(m), ns(m);
  IntegerVector which_one;
  for(unsigned j = 0; j < n; ++j){
    cur_seq = dna[j];
    std::string cur = cur_seq.get_cstring();
    for(unsigned i = 0; i < m; ++i){
      if(cur[i] == 'A'){ // must use single quote only
        ++na[i];
        continue;
      }
      
      if(cur[i] == 'T'){
        ++nt[i];
        continue;
      }
      
      if(cur[i] == 'C'){
        ++nc[i];
        continue;
      }
      
      if(cur[i] == 'G'){
        ++ng[i];
        continue;
      }
      
      if(cur[i] == '-'){
        ++nh[i];
        continue;
      }
      
      if(cur[i] == ' '){
        ++ns[i];
        continue;
      }
    }
  }
  
  IntegerVector cons(m);
  for(int k = 0; k < m; ++k){
    if(nt[k] > cons[k]){
      cons[k] = 1;
      continue;
    }
    
    if(nc[k] > cons[k]){
      cons[k] = 2;
      continue;
    }
    
    if(ng[k] > cons[k]){
      cons[k] = 3;
      continue;
    }
    
    if(nh[k] > cons[k]){
      cons[k] = 4;
      continue;
    }
  }
  cons = cons[cons != 4];
  CharacterVector seq_con = nucleotide[cons];
  int seq_len = seq_con.size();
  std::string str_con;
  std::string cur_nuc;
  for(int h = 0; h < seq_con.size(); ++h){
    cur_nuc = seq_con[h];
    str_con = str_con + cur_nuc; 
  }
  return(str_con);
}

// [[Rcpp::export]]
CharacterVector get_consensus(List dna){
  int n = dna.size();
  SEXP cur;
  CharacterVector result(dna.size());
  for(int i = 0; i < n; ++i){
    cur = dna[i];
    result[i] = single_get_consensus(cur);
  }
  return(result);
}

// [[Rcpp::export]]
SEXP assembler(CharacterVector seq_set, 
               IntegerVector count = 0,
               double min_pid = 90, 
               int min_ol_len = 8){
  //Rcout << "Performing all to all comparison" << "\n";
  seq_set.names() = seq(0, seq_set.size() - 1);
  if(count[0] == 0){
    count = IntegerVector(seq_set.size(), 1);
  }
  count.names() = seq(0, seq_set.size() - 1);
  CharacterVector original_set = seq_set;
  List ol = rcpp_overlapper(seq_set, "", true);
  ol = filter_and_order(ol, 95, 8);
  CharacterVector fully_covered = which_is_fully_covered(ol);
  
  //Rcout << "Removing sequences that fully represented by other sequences" << "\n";
  // To ol infromation
  CharacterVector ol_seq1 = ol["seq1"];
  CharacterVector ol_seq2 = ol["seq2"];
  LogicalVector keep = !(in(ol_seq1, fully_covered) | in(ol_seq2, fully_covered));
  IntegerVector keep_idx = rcpp_which(keep);
  IntegerVector col_idx = Rcpp::seq(0, ol.size()-1);
  ol = subset_df(ol, keep_idx, col_idx);
  
  // To sequence set
  CharacterVector seq_name = seq_set.names();
  seq_set = seq_set[!in(seq_name, fully_covered)];
  int remaining = ol_seq1.length();
  
  //Rcout << "initialising everything we needed in the loop" << "\n";
  CharacterVector selected_seq1, selected_seq2, used, check;
  std::string left_idx, right_idx, left_seq, right_seq, new_seq, new_seq_name;
  int seq1_ol_start, seq2_ol_start, seq2_len;
  LogicalVector remain;
  List selected, new_ol;
  CharacterVector new_ol_seq1, new_ol_seq2, new_covered, new_fully_covered, discard, new_seq_vec(1), new_seq_name_vec(1);
  int n = seq_set.size();
  CharacterVector cons;
  NumericVector support;
  
  bool has_ol = true;
  check = ol[0]; 
  
  if(check.size() == 0){
    return(List::create(cons, support));
  }
  
  int check_count = 0;
  while(has_ol){
    //Rcout << "selecting the best overlapping pair to combine as a new sequence" << "\n";
    selected = subset_df(ol, 1, col_idx); // 1-based row idx
    selected_seq1 = selected["seq1"];
    selected_seq2 = selected["seq2"];
    left_idx = std::string(selected_seq1[0]);
    right_idx = std::string(selected_seq2[0]);
    left_seq = std::string(seq_set[left_idx]);
    right_seq = std::string(seq_set[right_idx]);
    seq1_ol_start = selected["seq1_ol_start"];
    seq2_ol_start = selected["seq2_ol_start"];
    seq2_len = selected["seq2_len"];
    
    //Rcout << "Combining the sequence" << "\n";
    new_seq = left_seq.substr(0, seq1_ol_start) + right_seq.substr(seq2_ol_start,  seq2_len);
    new_seq_name = left_idx + "_" + right_idx;
    
    // converting std::string to CharacterVector
    new_seq_vec[0] = new_seq;
    new_seq_name_vec[0] = new_seq_name;
    new_seq_vec.names() = new_seq_name_vec;
    
    //Rcout << "Updating the remaining sequence set" << "\n";
    used = combine_two(selected["seq1"], selected["seq2"]);
    seq_name = seq_set.names();
    remain = !in(seq_name, used);
    seq_set = seq_set[remain];
    
    //Rcout << "Overlapping newseq and old seq" << "\n";
    new_ol = rcpp_overlapper(new_seq_vec, seq_set, false); //only find overlap between the new set and old set.
    new_ol = filter_and_order(new_ol, min_pid, min_ol_len);
    
    // Intergrating overlapping data in to the old ol data
    ol = rcpp_rbind(new_ol, ol);
    
    // Removing redundant data from ol
    new_fully_covered = which_is_fully_covered(new_ol);
    ol_seq1 = ol["seq1"];
    ol_seq2 = ol["seq2"];
    discard = combine_two(new_fully_covered, used);
    keep = !((in(ol_seq1, discard) | in(ol_seq2, discard)));
    keep_idx = rcpp_which(keep);
    ol = subset_df(ol, keep_idx, col_idx);
    
    // Removing sequences that have been covered by others to avoid repeated calculation
    seq_name = seq_set.names();
    seq_set = seq_set[!in(seq_name, new_fully_covered)];
    
    // Updating sequence information
    seq_set.push_back(new_seq, new_seq_name);
    
    IntegerVector ol_len = ol["aln_ol_len"];
    IntegerVector new_order = decreasing_order(ol_len);
    ol = subset_df(ol, new_order, col_idx);
    check = ol["seq1"]; 
    if(check.size()==0){
      has_ol = false;
    }
    ++check_count;
  }
  seq_name = seq_set.names();
  List members = string_split(seq_name, "_");
  int n_member = members.size();
  CharacterVector cur_member, cur_seqset, cur_ass_name, cur_ass, cur_msaview;
  List cur_ol;
  List msa_out(n_member);
  IntegerVector cur_count;
  support = NumericVector(n_member);
  
  for(int h = 0; h < n_member; ++h){
    cur_ass_name = seq_name[h];
    cur_ass = seq_set[cur_ass_name];
    cur_member = members[h];
    cur_seqset = original_set[cur_member];
    cur_count = count[cur_member];
    support[h] = Rcpp::sum(cur_count);
    cur_ol = rcpp_overlapper(cur_ass, cur_seqset, false, false);
    cur_msaview = msa_view(cur_ol["seq1_aln"], cur_ol["seq2_aln"]);
    msa_out[h] = replcae_space(cur_msaview);
  }
  
  cons = get_consensus(msa_out);
  cons.names() = seq_name;
  support.names() = seq_name;
  List out = List::create(cons, support);
  return(out);
}