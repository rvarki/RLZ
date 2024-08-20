#include <CLI/CLI.hpp>
#include "rlz_algo.h"

int main(int argc, char **argv) 
{
    CLI::App app("rlz - An implementation of RLZ that compresses a sequence file using a reference file.\n\nImplemented by Rahul Varki and Dhruv R. Makwana");

    std::string ref_file;
    std::string seq_file;
    bool decompress;
    std::string version = "Version: 1.0.0";

    app.add_option("-r,--ref", ref_file, "The reference file to be used for compression")->configurable()->required();
    app.add_option("-s,--seq", seq_file, "The sequence file to compress")->configurable()->required();
    app.add_flag("-d,--decompress", decompress, "Decompress the RLZ parse into the original sequence file")->configurable();
    app.set_version_flag("-v,--version", version);
    app.footer("Example usage:\n"
           "  Compress: ./rlz --ref reference.fasta --seq sequence.fasta\n"
           "  Decompress: ./rlz --ref reference.fasta --seq sequence.fasta --decompress");
    app.description("RLZ compression tool");
    CLI11_PARSE(app, argc, argv);

    if (decompress)
    {
        printf("Decompress the files\n");
        RLZ_Algo::RLZ main_parser(ref_file, seq_file);
        main_parser.decompress();
        //main_parser.clean();
    }
    else
    {
        printf("Compress the files\n");
        printf("Ref file: %s\n", ref_file.c_str());
        printf("Seq file: %s\n", seq_file.c_str());
        RLZ_Algo::RLZ main_parser(ref_file, seq_file);
        //main_parser.save_as_binary_file(seq_file);
        main_parser.load_bit_vectors();
        main_parser.compress();
    }

    return 0;
}