#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <cmath>
#include <stdexcept>

// [[Rcpp::export]]
Rcpp::List calculate_thermodynamics_cpp(std::string seq, double temperature, double Na, double oligo_conc) {
  
  // 1. Add a "guard clause" to handle empty or very short strings gracefully.
  if (seq.length() < 2) {
    return Rcpp::List::create(
      Rcpp::Named("Tm") = NA_REAL,
      Rcpp::Named("dG") = NA_REAL,
      Rcpp::Named("dH") = NA_REAL,
      Rcpp::Named("dS") = NA_REAL
    );
  }
  
  std::unordered_map<std::string, double> delH_map = {
    {"AA", -7.8}, {"AC", -5.9}, {"AG", -9.1}, {"AT", -8.3},
    {"CA", -9.0}, {"CC", -9.3}, {"CG", -16.3}, {"CT", -7.0},
    {"GA", -5.5}, {"GC", -8.0}, {"GG", -12.8}, {"GT", -7.8},
    {"TA", -7.8}, {"TC", -8.6}, {"TG", -10.4}, {"TT", -11.5}
  };
  
  std::unordered_map<std::string, double> delS_map = {
    {"AA", -21.9}, {"AC", -12.3}, {"AG", -23.5}, {"AT", -23.9},
    {"CA", -26.1}, {"CC", -23.2}, {"CG", -47.1}, {"CT", -19.7},
    {"GA", -13.5}, {"GC", -17.1}, {"GG", -31.9}, {"GT", -21.6},
    {"TA", -23.2}, {"TC", -22.9}, {"TG", -28.4}, {"TT", -36.4}
  };
  
  double total_delH = 0.0;
  double total_delS = 0.0;
  int len = seq.length();
  
  for (int i = 0; i < len - 1; ++i) {
    std::string dinucleotide = seq.substr(i, 2);
    
    // 2. Use the .at() method in a try-catch block for safe lookups.
    try {
      total_delH += delH_map.at(dinucleotide);
      total_delS += delS_map.at(dinucleotide);
    } catch (const std::out_of_range& oor) {
      Rcpp::stop("Invalid dinucleotide found: '" + dinucleotide + "'. Please ensure the input sequence only contains A, T, G, C.");
    }
  }
  
  // Add initiation parameters
  total_delH += 1.9;
  total_delS += -3.9;
  
  double dG = (total_delH * 1000 - (temperature + 273.15) * total_delS) / 1000;
  double Tm = (total_delH * 1000 / (total_delS + (1.9872 * log(oligo_conc / 4)))) - 273.15 + (16.6 * log10(Na));
  
  return Rcpp::List::create(
    Rcpp::Named("Tm") = Tm,
    Rcpp::Named("dG") = dG,
    Rcpp::Named("dH") = total_delH,
    Rcpp::Named("dS") = total_delS
  );
}