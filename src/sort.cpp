/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <CLI11.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <cstdint>
#include <filesystem> // Note that this requires at least gcc 9
#include <limits>

int main(int argc, char **argv) 
{
    CLI::App app("sort - Sorting RLZ factors.\n\nImplemented by Rahul Varki");

    std::string ref_file;
    std::string seq_file;
    std::string parse_file;
    std::string output;
    bool bit = false;
    bool rlz_repair = false;
    bool match_limit; // If max_len was set, we do not know whether the factors are maximal
    int verbosity = 0;
    bool csv = false; // Write out CSV of stats
    std::string version = "Version: 1.1.0";

    // RLZ sorting
    auto* rlz_cmd = app.add_subcommand("rlz", "Sorting suffixes directly from RLZ factors");
    rlz_cmd->add_option("-r,--ref", ref_file, "Reference file")->required();
    rlz_cmd->add_option("-p,--parse", parse_file, "RLZ parse file to sort")->required();
    rlz_cmd->add_option("-o,--output", parse_file, "Output prefix")->required();
    rlz_cmd->add_flag("--bit", bit, "Set if used during compression");
    rlz_cmd->add_flag("--repair", rlz_repair, "Set if used during compression");
    rlz_cmd->add_flag("--limit", match_limit, "Set if a match limit was specified during compression");
    rlz_cmd->add_flag("--csv", csv, "Output CSV file containing sorting statistics");
    rlz_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);
    
    // Text sorting
    auto* text_cmd = app.add_subcommand("text", "Sorting suffixes directly from text");
    text_cmd->add_option("-s,--seq", seq_file, "Sequence file to sort")->required();
    text_cmd->add_option("-o,--output", parse_file, "Output prefix")->required();
    text_cmd->add_flag("--csv", csv, "Output CSV file containing sorting statistics");
    text_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);

    // Choose between rlz or text sorting
    app.require_subcommand(1, 1); 

    // Set version flag
    app.set_version_flag("--version", version);

    // Footer Updates
    app.footer("Example usage:\n"
               "  ./sort rlz -r reference.fasta -p sequence.fasta.rlz -o path/to/output/prefix [--bit] [--repair] [--limit] [--csv]\n"
               "  ./sort -s sequence.fasta -o path/to/output/prefix [--csv]\n");

    CLI11_PARSE(app, argc, argv);

    if (verbosity == 2) {
        spdlog::set_level(spdlog::level::trace);
    }
    else if (verbosity == 1){
        spdlog::set_level(spdlog::level::debug);
    }
    else if (verbosity == 0){ 
        spdlog::set_level(spdlog::level::info);
    }

    if (rlz_cmd->parsed())
    {
        if (bit)
        {
            if (rlz_repair)
            {
                return 0;
            }
        }
        else
        {
            if (rlz_repair)
            {
                return 0;
            }
        }
    }
    else if (text_cmd->parsed())
    {

    }
    else{ spdlog::error("Neither rlz or text sorting. This condition should be impossible!"); }


    return 0;
}