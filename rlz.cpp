#include <CLI11.hpp>
#include "rlz_algo.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <filesystem> // Note that this requires at least gcc 9

int main(int argc, char **argv) 
{
    CLI::App app("rlz - An implementation of RLZ that compresses a sequence file using a reference file.\n\nImplemented by Rahul Varki");

    std::string ref_file;
    std::string seq_file;
    bool decompress = false;
    bool verbose = false;
    bool clean = false;
    int threads = 1;
    std::string version = "Version: 1.0.0";

    app.add_option("-r,--ref", ref_file, "The reference file to be used for compression")->configurable()->required();
    app.add_option("-s,--seq", seq_file, "The sequence file to compress")->configurable()->required();
    app.add_option("-t,--threads", threads, "Number of threads")->configurable();
    app.add_flag("-d,--decompress", decompress, "Decompress the RLZ parse into the original sequence file")->configurable();
    app.add_flag("--verbose", verbose, "Verbose output")->configurable();
    app.add_flag("--clean", clean, "Remove RLZ output files")->configurable();
    app.set_version_flag("-v,--version", version);
    app.footer("Example usage:\n"
           "  Compress: ./rlz --ref reference.fasta --seq sequence.fasta\n"
           "  Decompress: ./rlz --decompress --ref reference.fasta --seq sequence.fasta");
    app.description("RLZ compression tool");
    CLI11_PARSE(app, argc, argv);

    if (verbose) { spdlog::set_level(spdlog::level::debug); }

    if (clean) {RLZ main_parser(ref_file, seq_file); main_parser.clean(); return 0;}

    if (decompress)
    {
        spdlog::info("Starting to decompress the compressed sequence file");
        spdlog::stopwatch sw;
        spdlog::stopwatch sw_parser;
        RLZ main_parser(ref_file, seq_file);
        auto sw_parser_elapsed = sw_parser.elapsed();
        spdlog::debug("Built main parser in {:.3} seconds", sw_parser_elapsed.count());
        spdlog::stopwatch sw_decompress;
        main_parser.decompress();
        auto sw_decompress_elapsed = sw_decompress.elapsed();
        auto elapsed = sw.elapsed();
        spdlog::debug("Decompression function finished in {:.3} seconds", sw_decompress_elapsed.count());
        spdlog::info("Finished decompressing the compressed sequence file");
        spdlog::info("Decompressed in {:.3} seconds", elapsed.count());

    }
    else
    {
        spdlog::info("Starting to compress the sequence file");
        spdlog::stopwatch sw;
        spdlog::info("The reference file provided: {}", ref_file);
        spdlog::info("The sequence file provided: {}", seq_file);
        spdlog::stopwatch sw_parser;
        RLZ main_parser(ref_file, seq_file);
        auto sw_parser_elapsed = sw_parser.elapsed();
        spdlog::debug("Built main parser in {:.3} seconds", sw_parser_elapsed.count());
        spdlog::stopwatch sw_load_bit;
        main_parser.load_bit_vectors();
        auto sw_load_bit_elapsed = sw_load_bit.elapsed();
        spdlog::debug("Loaded bit vectors in {:.3} seconds", sw_load_bit_elapsed.count());
        spdlog::stopwatch sw_compress;
        main_parser.compress(threads);
        auto sw_compress_elapsed = sw_compress.elapsed();
        auto elapsed = sw.elapsed();
        spdlog::debug("Compression function finished in {:.3} seconds", sw_compress_elapsed.count());
        spdlog::info("Finished compressing the sequence file");
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

    return 0;
}