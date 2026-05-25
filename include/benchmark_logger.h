/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef BENCHMARK_LOGGER_H
#define BENCHMARK_LOGGER_H

#include <string>

void write_sort_benchmark_csv(const std::string& input_path, 
                         const std::string& config_name,
                         size_t text_size,
                         double sort_time,
                         size_t suffix_comps, 
                         size_t unit_comps,
                         double avg_unit_per_comp);

#endif