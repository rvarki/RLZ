/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <fstream>
#include <filesystem>
#include "spdlog/spdlog.h"

void write_sort_benchmark_csv(const std::string& input_path, 
                         const std::string& config_name,
                         size_t text_size,
                         double sort_time,
                         size_t suffix_comps, 
                         size_t unit_comps,
                         double avg_unit_per_comp)
{
    
    std::string csv_path = input_path + ".sort.csv";

    std::ofstream csv_file(csv_path);

    if (!csv_file) {
        spdlog::error("Failed to open CSV file for writing: {}", csv_path);
        return; 
    }

    csv_file << "Config,Text_Size,Sort_Time,Suffix_Cmps,Unit_Cmps,Avg_Units_Per_Suffix_Cmp\n";

    csv_file << config_name << ","
             << text_size << ","
             << sort_time << ","
             << suffix_comps << ","
             << unit_comps << ","
             << avg_unit_per_comp << "\n";
    
    csv_file.close();

    spdlog::info("Successfully wrote CSV at: {}", csv_path);
}