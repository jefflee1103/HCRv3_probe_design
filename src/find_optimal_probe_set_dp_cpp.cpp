#include <Rcpp.h>
#include <algorithm> // For std::stable_sort

using namespace Rcpp;

//' Find Optimal Probe Set using Dynamic Programming in C++
 //'
 //' This is a C++ implementation of the find_optimal_probe_set_dp function,
 //' designed for speed. It now includes the previously discussed stability
 //' and initialization fixes, and adds the 'score' column to the output.
 //'
 //' @param probes_df A DataFrame of candidate probes for a single contiguous region.
 //' @param probe_spacing The minimum number of nucleotides between probes.
 //' @return A DataFrame containing the optimal set of non-overlapping probes.
 // [[Rcpp::export]]
 DataFrame find_optimal_probe_set_dp_cpp(DataFrame probes_df, int probe_spacing) {
   // --- Data Extraction from DataFrame (Extract ALL columns) ---
   CharacterVector unique_id_vec = probes_df["unique_id"];
   NumericVector centre_pos_vec = probes_df["centre_pos"];
   IntegerVector start_vec = probes_df["start"];
   IntegerVector end_vec = probes_df["end"];
   CharacterVector target_sequence_vec = probes_df["target_sequence"];
   IntegerVector length_vec = probes_df["length"];
   CharacterVector rev_comp_vec = probes_df["rev_comp"];
   NumericVector GC_content_vec = probes_df["GC_content"];
   NumericVector A_content_vec = probes_df["A_content"];
   NumericVector C_content_vec = probes_df["C_content"];
   NumericVector dG_vec = probes_df["dG"];
   NumericVector Tm_vec = probes_df["Tm"];
   NumericVector dG_1st_half_vec = probes_df["dG_1st_half"];
   NumericVector Tm_1st_half_vec = probes_df["Tm_1st_half"];
   NumericVector dG_2nd_half_vec = probes_df["dG_2nd_half"];
   NumericVector Tm_2nd_half_vec = probes_df["Tm_2nd_half"];
   LogicalVector passed_a_comp_vec = probes_df["passed_a_comp"];
   LogicalVector passed_c_comp_vec = probes_df["passed_c_comp"];
   LogicalVector passed_a_stack_vec = probes_df["passed_a_stack"];
   LogicalVector passed_c_stack_vec = probes_df["passed_c_stack"];
   LogicalVector passed_c_spec_stack_vec = probes_df["passed_c_spec_stack"];
   NumericVector dG_dev_vec = probes_df["dG_deviation"];
   NumericVector dG_deviation_halves_vec = probes_df["dG_deviation_halves"];
   CharacterVector chrom_vec = probes_df["chrom"];
   CharacterVector region_id_vec = probes_df["region_id"];
   
   int n = probes_df.nrows();
   
   // --- Guard Clause ---
   if (n < 2) {
     if (n == 1) {
       // If there's one probe, return it with the score column added
       probes_df["score"] = 5000 - dG_dev_vec[0];
     }
     return probes_df;
   }
   
   // --- Indexing for Sorting ---
   IntegerVector p_idx(n);
   for(int i = 0; i < n; ++i) {
     p_idx[i] = i;
   }
   
   // Use std::stable_sort to match R's `arrange` behavior
   std::stable_sort(p_idx.begin(), p_idx.end(), [&](int i, int j) {
     return end_vec[i] < end_vec[j];
   });
   
   // Create sorted vectors based on p_idx
   IntegerVector s_start(n), s_end(n);
   NumericVector s_score(n);
   for(int i = 0; i < n; ++i) {
     s_start[i] = start_vec[p_idx[i]];
     s_end[i] = end_vec[p_idx[i]];
     s_score[i] = 5000 - dG_dev_vec[p_idx[i]];
   }
   
   // --- Predecessor Calculation (p array) ---
   // Initialize with -1 to correctly handle no-predecessor case
   IntegerVector p(n, -1);
   for (int i = 1; i < n; ++i) {
     int required_end = s_start[i] - probe_spacing;
     if (required_end < s_end[0]) {
       p[i] = -1; // Explicitly set to -1
       continue;
     }
     
     int low = 0, high = i - 1;
     int j = -1;
     while (low <= high) {
       int mid = low + (high - low) / 2;
       if (s_end[mid] <= required_end) {
         j = mid;
         low = mid + 1;
       } else {
         high = mid - 1;
       }
     }
     p[i] = j;
   }
   
   // --- Dynamic Programming ---
   NumericVector dp_scores(n + 1, 0.0);
   IntegerVector dp_counts(n + 1, 0);
   
   for (int i = 0; i < n; ++i) {
     int pred_idx = p[i] + 1;
     double score_with_i = s_score[i] + dp_scores[pred_idx];
     int count_with_i = 1 + dp_counts[pred_idx];
     
     double score_without_i = dp_scores[i];
     int count_without_i = dp_counts[i];
     
     if (count_with_i > count_without_i) {
       dp_counts[i + 1] = count_with_i;
       dp_scores[i + 1] = score_with_i;
     } else if (count_with_i == count_without_i && score_with_i > score_without_i) {
       dp_counts[i + 1] = count_with_i;
       dp_scores[i + 1] = score_with_i;
     } else {
       dp_counts[i + 1] = count_without_i;
       dp_scores[i + 1] = score_without_i;
     }
   }
   
   // --- Backtracking ---
   IntegerVector selected_indices;
   int i = n - 1;
   while (i >= 0) {
     int pred_idx = p[i] + 1;
     double score_with_i = s_score[i] + dp_scores[pred_idx];
     int count_with_i = 1 + dp_counts[pred_idx];
     
     if (dp_counts[i + 1] == count_with_i && std::abs(dp_scores[i + 1] - score_with_i) < 1e-9) {
       selected_indices.push_back(p_idx[i]);
       i = p[i];
     } else {
       i = i - 1;
     }
   }
   
   std::reverse(selected_indices.begin(), selected_indices.end());
   
   // --- Build Final DataFrame ---
   // Create new vectors to hold the data for the selected rows
   CharacterVector out_unique_id, out_target_sequence, out_rev_comp, out_chrom, out_region_id;
   IntegerVector out_start, out_end, out_length;
   NumericVector out_centre_pos, out_GC_content, out_A_content, out_C_content, out_dG, out_Tm;
   NumericVector out_dG_1st_half, out_Tm_1st_half, out_dG_2nd_half, out_Tm_2nd_half;
   NumericVector out_dG_deviation, out_dG_deviation_halves, out_score; // Added out_score
   LogicalVector out_passed_a_comp, out_passed_c_comp, out_passed_a_stack, out_passed_c_stack, out_passed_c_spec_stack;
   
   // Loop through the indices of the selected probes and populate the output vectors
   for (int index : selected_indices) {
     out_unique_id.push_back(unique_id_vec[index]);
     out_centre_pos.push_back(centre_pos_vec[index]);
     out_start.push_back(start_vec[index]);
     out_end.push_back(end_vec[index]);
     out_target_sequence.push_back(target_sequence_vec[index]);
     out_length.push_back(length_vec[index]);
     out_rev_comp.push_back(rev_comp_vec[index]);
     out_GC_content.push_back(GC_content_vec[index]);
     out_A_content.push_back(A_content_vec[index]);
     out_C_content.push_back(C_content_vec[index]);
     out_dG.push_back(dG_vec[index]);
     out_Tm.push_back(Tm_vec[index]);
     out_dG_1st_half.push_back(dG_1st_half_vec[index]);
     out_Tm_1st_half.push_back(Tm_1st_half_vec[index]);
     out_dG_2nd_half.push_back(dG_2nd_half_vec[index]);
     out_Tm_2nd_half.push_back(Tm_2nd_half_vec[index]);
     out_passed_a_comp.push_back(passed_a_comp_vec[index]);
     out_passed_c_comp.push_back(passed_c_comp_vec[index]);
     out_passed_a_stack.push_back(passed_a_stack_vec[index]);
     out_passed_c_stack.push_back(passed_c_stack_vec[index]);
     out_passed_c_spec_stack.push_back(passed_c_spec_stack_vec[index]);
     out_dG_deviation.push_back(dG_dev_vec[index]);
     out_dG_deviation_halves.push_back(dG_deviation_halves_vec[index]);
     out_chrom.push_back(chrom_vec[index]);
     out_region_id.push_back(region_id_vec[index]);
     out_score.push_back(5000 - dG_dev_vec[index]); // Add the score
   }
   
   // Create and return the new DataFrame with all original columns plus the score
   return DataFrame::create(
     Named("unique_id") = out_unique_id,
     Named("centre_pos") = out_centre_pos,
     Named("start") = out_start,
     Named("end") = out_end,
     Named("target_sequence") = out_target_sequence,
     Named("length") = out_length,
     Named("rev_comp") = out_rev_comp,
     Named("GC_content") = out_GC_content,
     Named("A_content") = out_A_content,
     Named("C_content") = out_C_content,
     Named("dG") = out_dG,
     Named("Tm") = out_Tm,
     Named("dG_1st_half") = out_dG_1st_half,
     Named("Tm_1st_half") = out_Tm_1st_half,
     Named("dG_2nd_half") = out_dG_2nd_half,
     Named("Tm_2nd_half") = out_Tm_2nd_half,
     Named("passed_a_comp") = out_passed_a_comp,
     Named("passed_c_comp") = out_passed_c_comp,
     Named("passed_a_stack") = out_passed_a_stack,
     Named("passed_c_stack") = out_passed_c_stack,
     Named("passed_c_spec_stack") = out_passed_c_spec_stack,
     Named("dG_deviation") = out_dG_deviation,
     Named("dG_deviation_halves") = out_dG_deviation_halves,
     Named("chrom") = out_chrom,
     Named("region_id") = out_region_id,
     Named("score") = out_score // Appended the score column
   );
 }
 