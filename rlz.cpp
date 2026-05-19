/*
* RLZ - Compute the RLZ parse of a sequence file using a reference file
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#include <CLI11.hpp>
#include "rlz_algo_bit.h"
#include "rlz_algo_char.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <cstdint>
#include <filesystem> // Note that this requires at least gcc 9

template <typename int_t>
void run_bit_decompression(const std::string& ref_file, const std::string& parse_file)
{
    spdlog::debug("Starting to decompress the compressed sequence file");
    spdlog::stopwatch sw;
    spdlog::stopwatch sw_parser;
    RLZ_BIT<int_t> main_parser(ref_file);
    auto sw_parser_elapsed = sw_parser.elapsed();
    spdlog::debug("Built main parser in {:.3} seconds", sw_parser_elapsed.count());
    spdlog::debug("Starting to store the reference file as a bit vector");
    spdlog::stopwatch sw_ref;
    main_parser.load_file_to_bit_vector(ref_file, main_parser.ref_bit_array);
    auto sw_ref_elapsed = sw_ref.elapsed();
    spdlog::debug("Loaded file in {:.3} seconds", sw_ref_elapsed.count());
    spdlog::stopwatch sw_decompress;
    main_parser.decompress(parse_file);
    auto sw_decompress_elapsed = sw_decompress.elapsed();
    auto elapsed = sw.elapsed();
    spdlog::debug("Decompression function finished in {:.3} seconds", sw_decompress_elapsed.count());
    spdlog::debug("Finished decompressing the compressed sequence file");
    spdlog::info("Decompressed in {:.3} seconds", elapsed.count());
}

template <typename int_t>
void run_bit_compression(const std::string& ref_file, const std::string& seq_file, int threads, int max_len)
{
    spdlog::debug("Starting to compress the sequence file");
    spdlog::stopwatch sw;
    spdlog::debug("The reference file provided: {}", ref_file);
    spdlog::debug("The sequence file provided: {}", seq_file);
    spdlog::stopwatch sw_parser;
    // Stream the sequence file
    RLZ_BIT<int_t> main_parser(ref_file);
    auto sw_parser_elapsed = sw_parser.elapsed();
    spdlog::debug("Built main parser in {:.3} seconds", sw_parser_elapsed.count());
    spdlog::stopwatch sw_ref;
    spdlog::debug("Starting to store the reference file as a bit vector");
    main_parser.load_reverse_file_to_bit_vector(ref_file, main_parser.ref_bit_array);
    auto sw_ref_elapsed = sw_ref.elapsed();
    spdlog::debug("Loaded file in {:.3} seconds", sw_ref_elapsed.count());
    spdlog::stopwatch sw_compress;
    main_parser.compress(threads, seq_file);
    auto sw_compress_elapsed = sw_compress.elapsed();
    spdlog::debug("Compression function finished in {:.3} seconds", sw_compress_elapsed.count());
    auto elapsed = sw.elapsed();
    spdlog::debug("Finished compressing the sequence file");
    spdlog::info("Compressed in {:.3} seconds", elapsed.count());
    spdlog::info("#############################################################");
    spdlog::info("File Size Statistics:");
    uintmax_t ref_size = std::filesystem::file_size(ref_file); //bytes
    uintmax_t seq_size = std::filesystem::file_size(seq_file); //bytes
    uintmax_t parse_size = std::filesystem::file_size(seq_file + ".rlz"); //bytes
    double comp_ratio = static_cast<double>(ref_size + parse_size) / 
                static_cast<double>(ref_size + seq_size) * 100;
    spdlog::info("The reference (ref) file provided: {} is {} bytes", ref_file, ref_size);
    spdlog::info("The sequence (seq) file provided: {} is {} bytes", seq_file, seq_size);
    spdlog::info("The parse (parse) file created: {} is {} bytes", seq_file + ".rlz", parse_size);
    spdlog::info("Compression ratio [((ref + parse)/(ref + seq)) * 100]: {:.3}%", comp_ratio);
}

template <typename int_t>
void run_char_decompression(const std::string& ref_file, const std::string& parse_file)
{
    spdlog::debug("Starting to decompress the compressed sequence file");
    spdlog::stopwatch sw;
    spdlog::stopwatch sw_parser;
    RLZ_CHAR<int_t> main_parser(ref_file);
    auto sw_parser_elapsed = sw_parser.elapsed();
    spdlog::debug("Built main parser in {:.3} seconds", sw_parser_elapsed.count());
    spdlog::debug("Starting to read the reference file");
    spdlog::stopwatch sw_ref;
    main_parser.load_file_to_string(ref_file, main_parser.ref_content);
    auto sw_ref_elapsed = sw_ref.elapsed();
    spdlog::debug("Loaded file in {:.3} seconds", sw_ref_elapsed.count());
    spdlog::stopwatch sw_decompress;
    main_parser.decompress(parse_file);
    auto sw_decompress_elapsed = sw_decompress.elapsed();
    auto elapsed = sw.elapsed();
    spdlog::debug("Decompression function finished in {:.3} seconds", sw_decompress_elapsed.count());
    spdlog::debug("Finished decompressing the compressed sequence file");
    spdlog::info("Decompressed in {:.3} seconds", elapsed.count());
}

template <typename int_t>
void run_char_compression(const std::string& ref_file, const std::string& seq_file, int threads, int max_len)
{
    spdlog::debug("Starting to compress the sequence file");
    spdlog::stopwatch sw;
    spdlog::debug("The reference file provided: {}", ref_file);
    spdlog::debug("The sequence file provided: {}", seq_file);
    spdlog::stopwatch sw_parser;
    // Stream the sequence file
    RLZ_CHAR<int_t> main_parser(ref_file);
    auto sw_parser_elapsed = sw_parser.elapsed();
    spdlog::debug("Built main parser in {:.3} seconds", sw_parser_elapsed.count());
    spdlog::stopwatch sw_ref;
    spdlog::debug("Starting to read the reference file");
    main_parser.load_reverse_file_to_string(ref_file, main_parser.ref_content);
    auto sw_ref_elapsed = sw_ref.elapsed();
    spdlog::debug("Finished reading file in {:.3} seconds", sw_ref_elapsed.count());
    spdlog::stopwatch sw_compress;
    main_parser.compress(threads, seq_file);
    auto sw_compress_elapsed = sw_compress.elapsed();
    spdlog::debug("Compression function finished in {:.3} seconds", sw_compress_elapsed.count());
    auto elapsed = sw.elapsed();
    spdlog::debug("Finished compressing the sequence file");
    spdlog::info("Compressed in {:.3} seconds", elapsed.count());
    spdlog::info("#############################################################");
    spdlog::info("File Size Statistics:");
    uintmax_t ref_size = std::filesystem::file_size(ref_file); //bytes
    uintmax_t seq_size = std::filesystem::file_size(seq_file); //bytes
    uintmax_t parse_size = std::filesystem::file_size(seq_file + ".rlz"); //bytes
    double comp_ratio = static_cast<double>(ref_size + parse_size) / 
                static_cast<double>(ref_size + seq_size) * 100;
    spdlog::info("The reference (ref) file provided: {} is {} bytes", ref_file, ref_size);
    spdlog::info("The sequence (seq) file provided: {} is {} bytes", seq_file, seq_size);
    spdlog::info("The parse (parse) file created: {} is {} bytes", seq_file + ".rlz", parse_size);
    spdlog::info("Compression ratio [((ref + parse)/(ref + seq)) * 100]: {:.3}%", comp_ratio);
}

int main(int argc, char **argv) 
{
    CLI::App app("rlz - An implementation of RLZ that compresses a sequence file using a reference file.\n\nImplemented by Rahul Varki");

    std::string ref_file;
    std::string seq_file;
    std::string parse_file;
    bool decompress = false;
    int verbosity = 0;
    bool clean = false;
    bool bit = false;
    uint64_t max_len = 0; // 0 means not set
    int threads = 1;
    std::string version = "Version: 1.0.0";
    
    // Compress Subcommand
    auto* compress_cmd = app.add_subcommand("compress", "Compress a sequence file using RLZ");
    compress_cmd->add_option("-r,--ref", ref_file, "Reference file")->required();
    compress_cmd->add_option("-s,--seq", seq_file, "Sequence file to compress")->required();
    compress_cmd->add_option("-t,--threads", threads, "Number of threads to use")->default_val(1);
    compress_cmd->add_option("-l, --len", max_len, "Maximum length a match can span")->check(CLI::Range(1UL, UINT64_MAX));
    compress_cmd->add_flag("--bit", bit, "Experimental: Set if ref lacks unique sequence chars");
    compress_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);   

    // Decompress Subcommand
    auto* decompress_cmd = app.add_subcommand("decompress", "Decompress an RLZ parse file");
    decompress_cmd->add_option("-r,--ref", ref_file, "Reference file")->required();
    decompress_cmd->add_option("-p,--parse", parse_file, "RLZ parse file to decompress")->required();
    decompress_cmd->add_option("-l, --len", max_len, "Maximum length a match can span (must be specified if used for compression)")->check(CLI::Range(1UL, UINT64_MAX));
    decompress_cmd->add_flag("--bit", bit, "Experimental: Set if ref lacks unique sequence chars (must be specified if used for compression)");
    decompress_cmd->add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = info, 1 = debug, 2 = trace)")->check(CLI::Range(0, 2))->default_val(0);   

    // Choose between compression or decompression subcommand
    app.require_subcommand(1, 1); 

    // Set version flag
    app.set_version_flag("--version", version);

    // Footer Updates
    app.footer("Example usage:\n"
               "  ./rlz compress -r reference.fasta -s sequence.fasta [--bit] [--len [int] ]\n"
               "  ./rlz decompress -r reference.fasta -p sequence.fasta.rlz [--bit] [--len [int] ]\n");

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
    
    if (compress_cmd->parsed())
    {
        // Encode with "bit" level compression
        if (bit)
        {
            spdlog::info("Bit alphabet compression enabled");
            uintmax_t upper_bound_bits = 0;

            // Determine size of parse entries
            if (max_len > 0) {
                spdlog::info("Using the specified match-length constraint to determine entry size");
                upper_bound_bits = max_len * 8;
            } else {
                spdlog::info("Using the reference size constraint to determine entry size");
                uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
                upper_bound_bits = ref_size * 8;
            }

            if (upper_bound_bits <= UINT8_MAX) { spdlog::info("Encoding entries with uint8_t"); run_bit_compression<uint8_t>(ref_file, seq_file, threads, max_len * 8); }
            else if (upper_bound_bits <= UINT16_MAX) { spdlog::info("Encoding entries with uint16_t"); run_bit_compression<uint16_t>(ref_file, seq_file, threads, max_len * 8); }
            else if (upper_bound_bits <= UINT32_MAX) { spdlog::info("Encoding entries with uint32_t"); run_bit_compression<uint32_t>(ref_file, seq_file, threads, max_len * 8); }
            else if (upper_bound_bits <= UINT64_MAX) { spdlog::info("Encoding entries with uint64_t"); run_bit_compression<uint64_t>(ref_file, seq_file, threads, max_len * 8); }
            else{
                spdlog::error("Determined size is too large! Check your reference file or maximum match length parameter.");
                exit(1);
            }
        }
        // Encode with regular alphabet compression
        else
        {
            spdlog::info("Original alphabet compression enabled");
            uintmax_t upper_bound = 0;

            // Determine size of parse entries
            if (max_len > 0) {
                spdlog::info("Using the specified match-length constraint to determine entry size");
                upper_bound = max_len;
            } else {
                spdlog::info("Using the reference size constraint to determine entry size");
                upper_bound = std::filesystem::file_size(ref_file); // bytes
            }

            if (upper_bound <= UINT8_MAX) { spdlog::info("Encoding entries with uint8_t"); run_bit_compression<uint8_t>(ref_file, seq_file, threads, max_len); }
            else if (upper_bound <= UINT16_MAX) { spdlog::info("Encoding entries with uint16_t"); run_char_compression<uint16_t>(ref_file, seq_file, threads, max_len); }
            else if (upper_bound <= UINT32_MAX) { spdlog::info("Encoding entries with uint32_t"); run_char_compression<uint32_t>(ref_file, seq_file, threads, max_len); }
            else if (upper_bound <= UINT64_MAX) { spdlog::info("Encoding entries with uint64_t"); run_char_compression<uint64_t>(ref_file, seq_file, threads, max_len); }
            else{
                spdlog::error("Determined size is too large! Check your reference file or maximum match length parameter.");
                exit(1);
            }
        }
    }
    else if (decompress_cmd->parsed())
    {
        // Encoded with "bit" level compression 
        if (bit)
        {
            spdlog::info("Bit alphabet decompression enabled");
            uintmax_t upper_bound_bits = 0;

            // Determine size of parse entries
            if (max_len > 0) {
                spdlog::info("Using the specified match-length constraint to determine entry size");
                upper_bound_bits = max_len * 8;
            } else {
                spdlog::info("Using the reference size constraint to determine entry size");
                uintmax_t ref_size = std::filesystem::file_size(ref_file); // bytes
                upper_bound_bits = ref_size * 8;
            }

            // Entries are decoded dynamically by upper bound specified
            if (upper_bound_bits <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_bit_decompression<uint8_t>(ref_file, parse_file); }
            else if (upper_bound_bits <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_bit_decompression<uint16_t>(ref_file, parse_file); }
            else if (upper_bound_bits <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_bit_decompression<uint32_t>(ref_file, parse_file); }
            else if (upper_bound_bits <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_bit_decompression<uint64_t>(ref_file, parse_file); }
            else{
                spdlog::error("Determined size is too large! Check your reference file or maximum match length parameter.");
                exit(1);
            }
        }
        // Encoded with regular alphabet compression
        else
        {
            spdlog::info("Original alphabet decompression enabled");
            uintmax_t upper_bound = 0;

            // Determine size of parse entries
            if (max_len > 0) {
                spdlog::info("Using the specified match-length constraint to determine entry size");
                upper_bound = max_len;
            } else {
                spdlog::info("Using the reference size constraint to determine entry size");
                upper_bound = std::filesystem::file_size(ref_file); // bytes
            }

            // Entries are decoded dynamically by upper bound specified
            if (upper_bound <= UINT8_MAX) { spdlog::info("Assuming entries were encoded with uint8_t"); run_char_decompression<uint8_t>(ref_file, parse_file); }
            else if (upper_bound <= UINT16_MAX) { spdlog::info("Assuming entries were encoded with uint16_t"); run_char_decompression<uint16_t>(ref_file, parse_file); }
            else if (upper_bound <= UINT32_MAX) { spdlog::info("Assuming entries were encoded with uint32_t"); run_char_decompression<uint32_t>(ref_file, parse_file); }
            else if (upper_bound <= UINT64_MAX) { spdlog::info("Assuming entries were encoded with uint64_t"); run_char_decompression<uint64_t>(ref_file, parse_file); }
            else{
                spdlog::error("Determined size is too large! Check your reference file or maximum match length parameter.");
                exit(1);
            }
        }
    }
    else{ spdlog::error("Neither compression or decompression. This condition should be impossible!"); }
    
    return 0;
}