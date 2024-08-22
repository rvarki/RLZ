#include <CLI11.hpp>
#include "rlz_algo.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

int main(int argc, char **argv) 
{
    CLI::App app("rlz - An implementation of RLZ that compresses a sequence file using a reference file.\n\nImplemented by Rahul Varki and Dhruv R. Makwana");

    std::string ref_file;
    std::string seq_file;
    bool decompress = false;
    bool verbose = false;
    int threads = 1;
    std::string version = "Version: 1.0.0";

    app.add_option("-r,--ref", ref_file, "The reference file to be used for compression")->configurable()->required();
    app.add_option("-s,--seq", seq_file, "The sequence file to compress")->configurable()->required();
    app.add_option("-t,--threads", threads, "Number of threads")->configurable();
    app.add_flag("-d,--decompress", decompress, "Decompress the RLZ parse into the original sequence file")->configurable();
    app.add_flag("--verbose", verbose, "Verbose output")->configurable();
    app.set_version_flag("-v,--version", version);
    app.footer("Example usage:\n"
           "  Compress: ./rlz --ref reference.fasta --seq sequence.fasta\n"
           "  Decompress: ./rlz --decompress --ref reference.fasta --seq sequence.fasta");
    app.description("RLZ compression tool");
    CLI11_PARSE(app, argc, argv);

    if (verbose) { spdlog::set_level(spdlog::level::debug); }

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
    }

    return 0;
}