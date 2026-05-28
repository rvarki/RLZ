/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <fstream>
#include <filesystem>
#include <string>
#include "spdlog/spdlog.h"

/**
 * @brief Appends sorting benchmark stats to a JSON Lines (.jsonl) file.
 * Writes standard hardware and execution metrics for all sorting algorithms. 
 * If the run utilizes the RLZ architecture, it optionally appends a sparse 
 * sub-object containing optimization benchmarks.
 *
 * @param input_path         Base path of the input text (used to name the output file).
 * @param config_name        Identifier for the algorithm/configuration (e.g., "text_naive", "rlz_greedy").
 * @param text_size          Total size of the uncompressed reference text in bytes.
 * @param sort_time          Total wall-clock execution time of the sorting phase (seconds).
 * @param suffix_comps       Total number of suffix-level comparisons evaluated.
 * @param unit_comps         Total number of individual character/boundary comparisons (LCE cost).
 * @param avg_unit_per_comp  Average char/boundary depth per suffix comparison .
 * @param is_rlz             Flag indicating if RLZ-specific stats should be attached.
 * @param total_factors      [RLZ] The total number of parsed macro-blocks/factors.
 * @param interval_hits      [RLZ] Times the O(1) disjoint interval shortcut bypassed LCE jumps.
 * @param backbone_hits      [RLZ] Times incomplete suffixes were resolved via the O(1) backbone trick.
 */


void write_sort_benchmark_jsonl(const std::string& input_path, 
                                const std::string& config_name,
                                size_t text_size,
                                double sort_time,
                                size_t suffix_comps, 
                                size_t unit_comps,
                                double avg_unit_per_comp,
                                bool is_rlz,
                                size_t total_factors,
                                size_t interval_hits,
                                size_t backbone_hits)
{
    // Switch to .jsonl to indicate line-delimited JSON
    std::string jsonl_path = input_path + ".sort.jsonl";

    std::ofstream json_file(jsonl_path);

    if (!json_file) {
        spdlog::error("Failed to open JSONL file for writing: {}", jsonl_path);
        return; 
    }

    // Construct the standard baseline metrics
    json_file << "{"
              << "\"Config\": \"" << config_name << "\", " 
              << "\"Text_Size\": " << text_size << ", "
              << "\"Sort_Time\": " << sort_time << ", "
              << "\"Suffix_Cmps\": " << suffix_comps << ", "
              << "\"Unit_Cmps\": " << unit_comps << ", " 
              << "\"Avg_Units_Per_Suffix_Cmp\": " << avg_unit_per_comp;

    // Append sparse RLZ metrics only if this was an RLZ run
    if (is_rlz) {
        json_file << ", \"RLZ_Metrics\": {"
                  << "\"Total_Factors\": " << total_factors << ", "
                  << "\"Interval_Shortcut_Hits\": " << interval_hits << ", "
                  << "\"Backbone_Shortcut_Hits\": " << backbone_hits
                  << "}";
    }
    
    // Close the JSON object and add the required newline for the JSONL specification
    json_file << "}\n"; 
    
    json_file.close();

    spdlog::info("Successfully wrote benchmark metrics to: {}", jsonl_path);
}