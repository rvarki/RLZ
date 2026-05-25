/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef SORT_ALGO_TEXT_H
#define SORT_ALGO_TEXT_H

#include <vector>
#include <string>
#include <string_view>

class TEXT_SORT 
{
    public:
        std::string seq_content;
        std::vector<size_t> suffix_array;

        TEXT_SORT(const std::string seq_file);
        ~TEXT_SORT();

        bool comparator(std::string_view a, std::string_view b);
        void buildSuffixArray(bool csv);
        void writeSuffixArray(const std::string seq_file);
};

#endif  // SORT_ALGO_TEXT_H