/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include "sort_algo_text.h"
#include <iostream>
#include <vector>
#include <string>
#include <string_view> //GCC 7
#include <algorithm>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <fstream>
#include <system_error>
#include <filesystem>


/**
 * @brief Constructor of the TEXT_SORT
 * @param [in] seq_file [string] Path to the sequence file from which to construct the Suffix Array
 */

TEXT_SORT::TEXT_SORT(const std::string seq_file)
{
    spdlog::debug("Reading in sequence file");

    spdlog::stopwatch sw_convert;

    // Getting size of sequence file
    std::error_code ec;
    uintmax_t seq_size = std::filesystem::file_size(seq_file, ec);
    if (ec) {
        spdlog::error("Error getting file size for {}: {}", seq_file, ec.message());
        std::exit(EXIT_FAILURE);
    }

    // Opening sequence file
    std::ifstream seq(seq_file, std::ios::binary);
    if (!seq) {
        spdlog::error("Error opening {}", seq_file);
        std::exit(EXIT_FAILURE);
    }

    // Preloading size of sequence buffer
    seq_content.resize(seq_size);

    // Directly loading sequence into buffer
    if (!seq.read(&seq_content[0], seq_size)) {
        spdlog::error("Error reading data from {}", seq_file);
        std::exit(EXIT_FAILURE);
    }
    seq.close();

    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished reading sequence file in {:.3} seconds", sw_convert_elapsed.count());
}

/**
* @brief Destructor of TEXT_SORT class.
*
* Currently does nothing.
*
*/
TEXT_SORT::~TEXT_SORT(){}


/**
 * @brief Comparison function to use for sorting
 * 
 * @param [in] a [string_view] suffix a
 * @param [in] b [string_view] suffix b
 * 
 * @return a less than b
 */

bool TEXT_SORT::comparator(std::string_view a, std::string_view b){ return a < b; } 


/**
 * @brief Creates the Suffix Array of a text with naive method using custom comparison operator
 * 
 * Sorts the suffixes using a custom comparison operator. Worst case is O(N^2log(N)) since 
 * there are O(Nlog(N)) comparisons and each comparison takes O(N) time worst case. There 
 * are more efficient manners to create the suffix array, but the purpose is to compare 
 * sorting time with RLZ version, therefore an inefficient but compatible metho with RLZ was chosen.
 * 
 * @param [in] csv [bool] Whether to produce csv containing sort stats
 */

void TEXT_SORT::buildSuffixArray(bool csv) 
{
    spdlog::stopwatch sw_sort;

    size_t n = seq_content.size();
    suffix_array.resize(n);
    
    // Fill the SA with the indices to start
    for (size_t i = 0; i < n; i++) {
        suffix_array[i] = i;
    }

    // Read only view into original text
    std::string_view view(seq_content);

    // Sort the indices using the custom comparator
    std::sort(suffix_array.begin(), suffix_array.end(), [&](size_t a, size_t b) {
        return comparator(view.substr(a), view.substr(b)); // Extract the suffix in O(1) time 
    });

    auto sw_sort_elapsed = sw_sort.elapsed();
    spdlog::info("Finished sorting suffixes in {:.3} seconds", sw_sort_elapsed.count());
}


/**
 * @brief Writes the Suffix Array of the text to file
 * 
 * @param [in] seq_file [string] The sequence file name is used to create output filename
 */

void TEXT_SORT::writeSuffixArray(const std::string seq_file)
{
    spdlog::stopwatch sw_write;

    std::string out_file = seq_file + ".sa";

    std::ofstream out(out_file);
    if (!out) {
        spdlog::error("Error opening {}", out_file);
        std::exit(EXIT_FAILURE);
    }

    for (size_t offset : suffix_array){
        out << offset << "\n";
    }

    out.close();

    auto sw_write_elapsed = sw_write.elapsed();
    spdlog::info("Finished writing suffix array in {:.3} seconds", sw_write_elapsed.count());
}