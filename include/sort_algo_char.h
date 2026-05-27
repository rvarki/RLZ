/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef SORT_ALGO_CHAR_H
#define SORT_ALGO_CHAR_H

#include <string>
#include "utility.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <fstream>
#include <system_error>
#include <filesystem>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

template<typename int_t>
class RLZ_CHAR_SORT
{
    public:
        struct RLZ_Factor { int_t p; int_t l; };

        // Represents a suffix starting at a specific character in the parsed text
        // factor_idx: The 0-based index of the factor within the sequence of factors where the suffix begins
        // offset: How far within the factor where the suffix starts
        struct SuffixID { size_t factor_idx; int_t offset; };

        // Represents a resynchronized state for the sort
        // id: The underlying SuffixID coordinate
        // first_factor: The factor after resynchronization or the original factor if option not used
        // is_ind: Whether first_factor is indicative or not
        struct SortableSuffix { SuffixID id; RLZ_Factor first_factor; bool is_ind; };

        std::string ref_content;
        sort_csa_index_t csa_ref; 
        sort_lcp_index_t lcp_ref; 
        sort_rmq_index_t lcp_rmq; 
        std::vector<RLZ_Factor> rlz_factors;
        std::vector<SuffixID> sa_T; // Stores the suffix array as 2D representation 

        size_t original_size = 0; 

        RLZ_CHAR_SORT(const std::string ref_file, const std::string parse_file);
        ~RLZ_CHAR_SORT();

        void sort_naive(bool apply_resync = false);
        void sort_lcp_interval(bool apply_resync = false);
        void sort_induced(bool apply_resync = false);
        void stream_sa_to_file(const std::string out_file);
    
    private:
        int_t get_lce(int_t i, int_t j); 
        bool compare_suffixes(const SortableSuffix& a, const SortableSuffix& b);
        bool is_indicative(const RLZ_Factor& f);
        std::pair<int_t, int_t> get_sa_range(const RLZ_Factor& f);
        RLZ_Factor apply_resynchronization(const RLZ_Factor& f_i, const RLZ_Factor& f_next);
        std::vector<SortableSuffix> generate_all_suffixes(bool apply_resync);
};

/**
 * @brief Initializes the RLZ sequence sorting framework.
 * Loads the reference text and the binary sequence of RLZ factors. 
 * Automatically constructs the necessary succinct data structures 
 * (Compressed Suffix Array, LCP Array, and RMQ support) over the 
 * reference string to enable O(1) LCE queries.
 * @param [in] ref_file   [std::string] Path to the raw reference text file.
 * @param [in] parse_file [std::string] Path to the binary RLZ parse file.
 */

template <typename int_t>
RLZ_CHAR_SORT<int_t>::RLZ_CHAR_SORT(const std::string ref_file, const std::string parse_file)
{
    spdlog::debug("Reading in reference file");
    spdlog::stopwatch sw_convert;

    std::error_code ec;
    uintmax_t ref_size = std::filesystem::file_size(ref_file, ec);
    if (ec) {
        spdlog::error("Error getting file size for {}: {}", ref_file, ec.message());
        std::exit(EXIT_FAILURE);
    }

    std::ifstream ref(ref_file, std::ios::binary);
    if (!ref) {
        spdlog::error("Error opening {}", ref_file);
        std::exit(EXIT_FAILURE);
    }

    ref_content.resize(ref_size);

    if (!ref.read(&ref_content[0], ref_size)) {
        spdlog::error("Error reading data from {}", ref_file);
        std::exit(EXIT_FAILURE);
    }
    ref.close();

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished reading reference file in {:.3} seconds", sw_convert_elapsed.count());
    
    spdlog::debug("Constructing ISA, LCP, and RMQ from reference");
    spdlog::stopwatch sw_build;

    sdsl::construct_im(csa_ref, ref_content, 1);
    sdsl::construct_im(lcp_ref, ref_content, 1);
    lcp_rmq = sdsl::rmq_succinct_sada<>(&lcp_ref);

    auto sw_build_elapsed = sw_build.elapsed();
    spdlog::debug("Finished building data structures in {:.3} seconds", sw_build_elapsed.count());
    
    spdlog::debug("Reading in RLZ parse");
    spdlog::stopwatch sw_parse;

    uint64_t num_pairs; 

    std::ifstream parse(parse_file, std::ios::binary);
    if (!parse) {
        spdlog::error("Error opening {}", parse_file);
        std::exit(EXIT_FAILURE);
    }

    parse.read(reinterpret_cast<char*>(&num_pairs), sizeof(num_pairs));

    spdlog::debug("The RLZ parse contains {} factors", num_pairs);

    uintmax_t parse_size = std::filesystem::file_size(parse_file);
    assert(parse_size - sizeof(uint64_t) == num_pairs * sizeof(RLZ_Factor));

    rlz_factors.resize(num_pairs);
    parse.read(reinterpret_cast<char*>(rlz_factors.data()), num_pairs * sizeof(RLZ_Factor));

    auto sw_parse_elapsed = sw_parse.elapsed();
    spdlog::debug("Finished reading RLZ parse in {:.3} seconds", sw_parse_elapsed.count());
}

/**
 * @brief Default destructor for the RLZ_CHAR_SORT class.
 */

template<typename int_t>
RLZ_CHAR_SORT<int_t>::~RLZ_CHAR_SORT(){}

/**
 * @brief Computes the Longest Common Extension (LCE) between two reference positions.
 * Implements Equation 2.2 from the paper. Utilizes the Inverse Suffix Array (ISA) 
 * and a Range Minimum Query (RMQ) over the LCP array to determine the length of 
 * the longest common prefix between two suffixes of R in O(1) time (assuming no compression).
 * @param [in] i [int_t] First position in the reference text.
 * @param [in] j [int_t] Second position in the reference text.
 * @return [int_t] The LCE length.
 */

template<typename int_t>
int_t RLZ_CHAR_SORT<int_t>::get_lce(int_t i, int_t j) {
    if (i == j) return ref_content.size() - i;
    
    int_t ri = csa_ref.isa[i];
    int_t rj = csa_ref.isa[j];
    
    if (ri > rj) std::swap(ri, rj);
    
    size_t min_idx = lcp_rmq(ri + 1, rj);
    return lcp_ref[min_idx];
}

/**
 * @brief Lexicographically compares two factor-aligned suffixes of the text.
 * Implements Algorithm 3.1. Evaluates the lexicographical order of two suffixes 
 * by performing LCE queries on their corresponding reference segments. This allows 
 * the comparison to "jump" over matching segments of the reference text in O(1) time 
 * per factor transition, avoiding decompression and character-by-character scans.
 * @param [in] a [SortableSuffix] The first suffix object.
 * @param [in] b [SortableSuffix] The second suffix object.
 * @return [bool] True if Suffix A is strictly less than Suffix B.
 */

template<typename int_t>
bool RLZ_CHAR_SORT<int_t>::compare_suffixes(const SortableSuffix& a, const SortableSuffix& b) {
    int_t pa = a.first_factor.p;
    int_t pb = b.first_factor.p;
    int_t la = a.first_factor.l;
    int_t lb = b.first_factor.l;
    
    size_t idx_a = a.id.factor_idx + 1; // The next factor following a 
    size_t idx_b = b.id.factor_idx + 1; // The next factor following b
    
    while (true) {

        spdlog::trace("Comparing Factor A: ({},{}) with Factor B: ({},{})", pa, la, pb, lb);

        int_t k = get_lce(pa, pb);

        spdlog::trace("LCE value: {} ---- min len: {}", k, std::min(la,lb));

        if (k < std::min(la, lb)) {
            spdlog::trace("Mismatch: A: ['{}'] ---- B: ['{}']", ref_content[pa +k], ref_content[pb + k]);
            return ref_content[pa + k] < ref_content[pb + k];
        } 
        
        if (la < lb) {
            spdlog::trace("Factor A fully consumed --- Factor B partially consumed");
            if (idx_a >= rlz_factors.size()){
                spdlog::trace("Suffix A has no more factors"); 
                return true; 
            }
            // Adjust Factor B by la
            pb += la; 
            lb -= la;
            // Factor A is now next factor
            pa = rlz_factors[idx_a].p;
            la = rlz_factors[idx_a].l;
            idx_a++;
        } else if (lb < la) {
            spdlog::trace("Factor A partially consumed --- Factor B fully consumed");
            if (idx_b >= rlz_factors.size()){
                spdlog::trace("Suffix B has no more factors");  
                return false;
            }
            // Adjust Factor A by lb
            pa += lb; 
            la -= lb;
            // Factor B is now next factor 
            pb = rlz_factors[idx_b].p;
            lb = rlz_factors[idx_b].l;
            idx_b++;
        } else {
            spdlog::trace("Factor A fully consumed --- Factor B fully consumed");
            if (idx_a >= rlz_factors.size() && idx_b >= rlz_factors.size()) {
                spdlog::error("Suffixes A and B have no more factors! Should not be possible.");
                return false;
            }
            if (idx_a >= rlz_factors.size()){
                spdlog::trace("Suffix A has no more factors"); 
                return true;
            }
            if (idx_b >= rlz_factors.size()){
                spdlog::trace("Suffix B has no more factors"); 
                return false;
            }
            
            // Both factors are the next factor
            pa = rlz_factors[idx_a].p;
            la = rlz_factors[idx_a].l;
            idx_a++;
            pb = rlz_factors[idx_b].p;
            lb = rlz_factors[idx_b].l;
            idx_b++;
        }
    }
}

/**
 * @brief Determines if an RLZ factor uniquely identifies a match (Algorithm 3.2).
 * An indicative factor has a length strictly greater than the longest common prefix 
 * it shares with its immediate lexicographical neighbors in the reference. 
 * Therefore, it corresponds to a unique match (a single SA range) in R.
 * @param [in] f [RLZ_Factor] The RLZ factor to evaluate.
 * @return [bool] True if the factor is indicative, false otherwise.
 */

template<typename int_t>
bool RLZ_CHAR_SORT<int_t>::is_indicative(const RLZ_Factor& f) {
    size_t i = csa_ref.isa[f.p];
    
    int_t v_prev = (i > 0) ? lcp_ref[i] : 0;
    int_t v_next = (i < csa_ref.size() - 1) ? lcp_ref[i + 1] : 0;
    
    return f.l > std::max(v_prev, v_next);
}

/**
 * @brief Finds the Suffix Array (SA) range for a non-indicative factor (Algorithm 3.3).
 * Non-indicative factors correspond to matches that occur multiple times in the 
 * reference. This method scans the LCP array to find the consecutive bounds in 
 * the SA that contain all alternative valid occurrences of the factor.
 * @param [in] f [RLZ_Factor] The non-indicative RLZ factor.
 * @return [std::pair<int_t, int_t>] The [start, end] indices in the Suffix Array.
 */

template<typename int_t>
std::pair<int_t, int_t> RLZ_CHAR_SORT<int_t>::get_sa_range(const RLZ_Factor& f) {
    size_t i = csa_ref.isa[f.p];
    size_t sp = i;
    size_t ep = i;

    while (sp > 0 && lcp_ref[sp] >= f.l) sp--;
    while (ep < csa_ref.size() - 1 && lcp_ref[ep + 1] >= f.l) ep++;
    
    return {sp, ep};
}

/**
 * @brief Restores right-maximality to incomplete factors (Algorithm 3.4).
 * Incomplete factors lack right-maximality due to artificial boundaries created 
 * by the greedy parsing strategy. This method searches the SA range for an 
 * overlapping occurrence that allows the factor to seamlessly consume characters 
 * from the succeeding factor, maximizing its length and indicativity.
 * @param [in] f_i    [RLZ_Factor] The incomplete factor.
 * @param [in] f_next [RLZ_Factor] The factor immediately following f_i.
 * @return [RLZ_Factor] The updated right-maximal factor (or original if no overlap).
 */

template<typename int_t>
typename RLZ_CHAR_SORT<int_t>::RLZ_Factor RLZ_CHAR_SORT<int_t>::apply_resynchronization(const RLZ_Factor& f_i, const RLZ_Factor& f_next) {
    std::pair<int_t, int_t> range_i = get_sa_range(f_i);
    
    RLZ_Factor f_next_char = {f_next.p, 1};
    std::pair<int_t, int_t> range_next = get_sa_range(f_next_char);
    
    bool overlap_exists = false;
    int_t best_p = f_i.p; 
    int_t max_k = 0;
    
    for (size_t j = range_i.first; j <= range_i.second; ++j) {
        size_t target_pos = csa_ref[j] + f_i.l;
        
        if (target_pos < csa_ref.size()) {
            size_t rank = csa_ref.isa[target_pos];
            
            if (rank >= range_next.first && rank <= range_next.second) {
                int_t candidate_p = csa_ref[j];
                int_t current_k = get_lce(candidate_p + f_i.l, f_next.p);
                
                if (current_k > max_k) {
                    max_k = current_k;
                    best_p = candidate_p;
                    overlap_exists = true;
                    if (max_k >= f_next.l) break;
                }
            }
        }
    }
    
    if (overlap_exists) {
        int_t k_prime = std::min(max_k, f_next.l);
        return {best_p, static_cast<int_t>(f_i.l + k_prime)};
    }
    
    return f_i;
}

/**
 * @brief Generates the full set of character-level suffixes for the original text.
 * Unrolls the 1D array of RLZ factors into individual Suffix objects, assigning 
 * an internal offset to every character position. If resynchronization is enabled, 
 * it dynamically applies Algorithm 3.4 to incomplete suffixes (offset > 0) to 
 * bypass artificial boundaries before the sort begins.
 * @param [in] apply_resync [bool] Toggle to enable or bypass boundary repair.
 * @return [std::vector<SortableSuffix>] The pre-processed array of all suffixes.
 */

template<typename int_t>
std::vector<typename RLZ_CHAR_SORT<int_t>::SortableSuffix> RLZ_CHAR_SORT<int_t>::generate_all_suffixes(bool apply_resync) {
    spdlog::debug("Generating all character-level suffixes (Resync: {})", apply_resync);
    
    std::vector<SortableSuffix> all_suffixes;
    
    for (size_t i = 0; i < rlz_factors.size(); ++i) {
        for (size_t offset = 0; offset < rlz_factors[i].l; ++offset) {
            
            SortableSuffix suf;
            suf.id = {static_cast<size_t>(i), static_cast<int_t>(offset)};
            
            RLZ_Factor effective_first = {
                static_cast<int_t>(rlz_factors[i].p + offset), 
                static_cast<int_t>(rlz_factors[i].l - offset)
            };
            
            if (offset > 0 && apply_resync && i + 1 < rlz_factors.size()) {
                effective_first = apply_resynchronization(effective_first, rlz_factors[i+1]);
            }
            
            suf.first_factor = effective_first;
            suf.is_ind = is_indicative(effective_first);
            
            spdlog::trace("Considering factor: ({},{}) --- Indicative: {}", suf.first_factor.p, suf.first_factor.l, suf.is_ind);

            all_suffixes.push_back(suf);
        }
    }
    
    return all_suffixes;
}

/**
 * @brief Baseline sorting strategy (O(nlogn) factor comparisons).
 * Evaluates all suffixes indiscriminately using introsort combined with Algorithm 3.1. 
 * Does not utilize LCP interval pruning or induced sorting mechanics. Acts as the 
 * performance control for the benchmark.
 * @param [in] apply_resync [bool] Toggle to enable resynchronization preprocessing.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_naive(bool apply_resync) {
    spdlog::info("Executing Naive Sort (Resync: {})", apply_resync);
    spdlog::stopwatch sw_sort;
    
    std::vector<SortableSuffix> all_suffixes = generate_all_suffixes(apply_resync);

    std::sort(all_suffixes.begin(), all_suffixes.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        return compare_suffixes(a, b);
    });

    sa_T.reserve(all_suffixes.size());
    for(const auto& suf : all_suffixes) sa_T.push_back(suf.id);

    spdlog::info("Naive Sort completed in {:.3} seconds", sw_sort.elapsed().count());
}

/**
 * @brief Sorting strategy utilizing LCP Interval Pruning.
 * Sorts suffixes by first checking their precomputed indicativity flags. If two 
 * indicative suffixes are compared, their lexicographical order is resolved 
 * immediately in O(1) time via their ISA ranks, bypassing LCE jumps entirely.
 * @param [in] apply_resync [bool] Toggle to enable resynchronization preprocessing.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_lcp_interval(bool apply_resync) {
    spdlog::info("Executing LCP Interval Sort (Resync: {})", apply_resync);
    spdlog::stopwatch sw_sort;
    
    std::vector<SortableSuffix> all_suffixes = generate_all_suffixes(apply_resync);

    std::sort(all_suffixes.begin(), all_suffixes.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        if (a.is_ind && b.is_ind) {
            int_t rank_a = csa_ref.isa[a.first_factor.p];
            int_t rank_b = csa_ref.isa[b.first_factor.p];
            if (rank_a != rank_b) return rank_a < rank_b; 
        }
        return compare_suffixes(a, b);
    });

    sa_T.reserve(all_suffixes.size());
    for(const auto& suf : all_suffixes) sa_T.push_back(suf.id);

    spdlog::info("LCP Interval Sort completed in {:.3} seconds", sw_sort.elapsed().count());
}

/**
 * @brief Highly optimized sorting strategy leveraging Lemma 3.2.
 * Bins suffixes into two categories: Complete (offset == 0) and Incomplete (offset > 0).
 * Sorts the complete, highly indicative backbone first. It then leverages this stable 
 * structural backbone to dramatically reduce the search space and LCE jump overhead 
 * when resolving the more difficult incomplete suffixes.
 * @param [in] apply_resync [bool] Toggle to enable resynchronization preprocessing.
 * @return void
 */

template<typename int_t>
void RLZ_CHAR_SORT<int_t>::sort_induced(bool apply_resync) {
    spdlog::info("Executing Induced Sort (Resync: {})", apply_resync);
    spdlog::stopwatch sw_sort;
    
    std::vector<SortableSuffix> all_suffixes = generate_all_suffixes(apply_resync);
    
    std::vector<SortableSuffix> complete_bin;
    std::vector<SortableSuffix> incomplete_bin;
    
    for (const auto& suf : all_suffixes) {
        if (suf.id.offset == 0) complete_bin.push_back(suf);
        else incomplete_bin.push_back(suf);
    }

    std::sort(complete_bin.begin(), complete_bin.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        int_t rank_a = csa_ref.isa[a.first_factor.p];
        int_t rank_b = csa_ref.isa[b.first_factor.p];
        return (rank_a != rank_b) ? (rank_a < rank_b) : compare_suffixes(a, b);
    });

    std::sort(incomplete_bin.begin(), incomplete_bin.end(), [&](const SortableSuffix& a, const SortableSuffix& b) {
        return compare_suffixes(a, b);
    });

    sa_T.reserve(all_suffixes.size());
    
    size_t i = 0, j = 0;
    while (i < complete_bin.size() && j < incomplete_bin.size()) {
        if (compare_suffixes(complete_bin[i], incomplete_bin[j])) {
            sa_T.push_back(complete_bin[i++].id);
        } else {
            sa_T.push_back(incomplete_bin[j++].id);
        }
    }
    while (i < complete_bin.size()) sa_T.push_back(complete_bin[i++].id);
    while (j < incomplete_bin.size()) sa_T.push_back(incomplete_bin[j++].id);

    spdlog::info("Induced Sort completed in {:.3} seconds", sw_sort.elapsed().count());
}

/**
 * @brief Streams the 2D Suffix Array to disk as a 1D Suffix Array in text format.
 * Bypasses full 1D vector allocation to save memory (preventing OOM on large datasets).
 * Calculates the absolute 1D character index on the fly using a prefix sum of factor lengths.
 * Writes each absolute index as a standard text string followed by a newline.
 * * @param [in] parse_file [std::string] The input parse file path used to generate output path.
 */
template<typename int_t>
void RLZ_CHAR_SORT<int_t>::stream_sa_to_file(const std::string parse_file) {

    std::string out_file = parse_file + ".sa";

    spdlog::info("Streaming 1D suffix array directly to text file {}", out_file);
    spdlog::stopwatch sw_stream;

    std::ofstream out(out_file);
    if (!out) {
        spdlog::error("Error opening output file {}", out_file);
        std::exit(EXIT_FAILURE);
    }

    // Precompute absolute starting positions of each factor in T
    std::vector<size_t> factor_starts;
    factor_starts.reserve(rlz_factors.size());
    
    size_t current_text_pos = 0;
    for (const auto& factor : rlz_factors) {
        factor_starts.push_back(current_text_pos);
        current_text_pos += factor.l;
    }

    // Stream elements directly to the file, one per line.
    for (const auto& suf : sa_T) {
        size_t absolute_index = factor_starts[suf.factor_idx] + suf.offset;
        out << absolute_index << '\n';
    }

    out.close();
    spdlog::info("Finished streaming suffix array to file in {:.3f} seconds", sw_stream.elapsed().count());
}




#endif  // SORT_ALGO_CHAR_H