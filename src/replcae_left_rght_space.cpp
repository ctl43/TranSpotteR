// [[Rcpp::depends(RSeqAn)]]
// [[Rcpp::plugins(cpp14)]]

#include <iostream>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::StringVector replcae_space(std::vector<std::string>  x) {
  size_t sz=x.size();
  Rcpp::StringVector out(sz);
  for(unsigned p=0; p<sz;++p){
    std::string current = x[p];
    
    //std::string str ("Test string");
    for (unsigned i=0; i<current.length(); ++i)
    {
      if(current[i]=='-'){
        //x.replace(i, 1, " ");
        current[i]=' ';
      }else{
        break;
      }
    }
    
    for (unsigned j=current.length()-1; j>=0; j--)
    {
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

/*** R

*/
