/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <CLI11.hpp>
#include "sort_algo_bit.h"
#include "sort_algo_char.h"
#include "sort_algo_text.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <cstdint>
#include <filesystem> // Note that this requires at least gcc 9
#include <limits>

template <typename int_t>
void run_rlz_bit_sort(const std::string& ref_file, const std::string& parse_file, bool match_limit, bool csv)
{
    
}

template <typename int_t>
void run_rlz_char_sort(const std::string& ref_file, const std::string& parse_file, bool match_limit, bool csv)
{
    
}


void run_text_sort(const std::string& seq_file, bool csv)
{
    TEXT_SORT main_parser(seq_file);
    main_parser.buildSuffixArray(csv);
    main_parser.writeSuffixArray(seq_file);
}



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
    // rlz_cmd->add_option("-o,--output", parse_file, "Output prefix")->required();
    rlz_cmd->add_flag("--bit", bit, "Set if used during compression");
    rlz_cmd->add_flag("--repair", rlz_repair, "Set if used during compression");
    rlz_cmd->add_flag("--limit", match_limit, "Set if a match limit was specified during compression");
    rlz_cmd->add_flag("--csv", csv, "Output CSV file containing sorting statistics");
    rlz_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);
    
    // Text sorting
    auto* text_cmd = app.add_subcommand("text", "Sorting suffixes directly from text");
    text_cmd->add_option("-s,--seq", seq_file, "Sequence file to sort")->required();
    // text_cmd->add_option("-o,--output", parse_file, "Output prefix")->required();
    text_cmd->add_flag("--csv", csv, "Output CSV file containing sorting statistics");
    text_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);

    // Choose between rlz or text sorting
    app.require_subcommand(1, 1); 

    // Set version flag
    app.set_version_flag("--version", version);

    // Footer Updates
    app.footer("Example usage:\n"
               "  ./sort rlz -r reference.fasta -p sequence.fasta.rlz [--bit] [--repair] [--limit] [--csv]\n"
               "  ./sort -s sequence.fasta [--csv]\n");

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
            spdlog::info("Bit alphabet sorting enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference 
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
            uintmax_t ref_size_bits = ref_size * 8;

            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair)
            {
                if (ref_size_bits < std::numeric_limits<int>::max()){
                    spdlog::info("Assuming entries encoded with int");
                    run_rlz_bit_sort<int>(ref_file, parse_file, match_limit, csv);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    exit(1);
                }
                return 0;
            }
            // Entry size is determined by the size of the reference
            if (ref_size_bits <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_rlz_bit_sort<uint8_t>(ref_file, parse_file, match_limit, csv); }
            else if (ref_size_bits <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_rlz_bit_sort<uint16_t>(ref_file, parse_file, match_limit, csv); }
            else if (ref_size_bits <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_rlz_bit_sort<uint32_t>(ref_file, parse_file, match_limit, csv); }
            else if (ref_size_bits <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_rlz_bit_sort<uint64_t>(ref_file, parse_file, match_limit, csv); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                exit(1);
            }
        }
        else
        {
            spdlog::info("Original alphabet sorting enabled");

            // Cannot use max len to determine size because position of match can be anywhere on the reference
            spdlog::info("Using the reference size to determine entry size");
            uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes

            // Solely for RLZ-RePair which should actually takes entries as int
            if (rlz_repair)
            {
                if (ref_size < std::numeric_limits<int>::max()){
                    spdlog::info("Assuming entries encoded with int");
                    run_rlz_char_sort<int>(ref_file, parse_file, match_limit, csv);
                }
                else{
                    spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                    exit(1);
                }
                return 0;
            }
            // Entries is determined by the size of the reference
            if (ref_size <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_rlz_char_sort<uint8_t>(ref_file, parse_file, match_limit, csv); }
            else if (ref_size <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_rlz_char_sort<uint16_t>(ref_file, parse_file, match_limit, csv); }
            else if (ref_size <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_rlz_char_sort<uint32_t>(ref_file, parse_file, match_limit, csv); }
            else if (ref_size <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_rlz_char_sort<uint64_t>(ref_file, parse_file, match_limit, csv); }
            else{
                spdlog::error("Determined reference size is too large! Choose a smaller reference file.");
                exit(1);
            }
        }
    }
    else if (text_cmd->parsed())
    {
        spdlog::info("Text sorting enabled");
        run_text_sort(seq_file, csv);
    }
    else{ spdlog::error("Neither rlz or text sorting. This condition should be impossible!"); }


    return 0;
}