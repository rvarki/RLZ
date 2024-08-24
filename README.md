# Relative Lempel-Ziv (RLZ)

## Description

This software computes the Relative Lempel Ziv (RLZ) parse of the target sequence file using a reference file. The software should work for any type of file whether that be FASTA files, English files, etc... However, currently the decompression expects that the files were originally ASCII (8 bit) encoded.

## Algorithm Workflow

To compress the target sequence file in relation to a reference file, the software performs the following steps:

1. Converts both the reference and sequence file into its binary representation.
2. Build an FM-index using the reference bit vector.
3. In reverse, attempt to match bit by bit of the sequence bit vector against the reference bit representation using the backwards match of the FM-index. 
    - 3a. If match, check if next bit also matches (ex. 001: if 1 matches then check if 01 matches etc...)
    - 3b. If match and at end of sequence file, record (pos,len) pair
    - 3c. If mismatch, record (prev pos, len - 1) pair. Reset search from bit that caused mismatch.
4. All the (pos, len) pairs are then written to a file in order. This is the RLZ parse.

> [!NOTE]
> The pairs encoded in the RLZ parse are in reference to the reference file.

## Prerequisites

- [CMake](https://cmake.org/) 3.15 or higher.
- GCC 9+
- C++17-compatible compiler.
- OpenMP

## Getting Started

### Building the Project

```
git clone https://github.com/rvarki/RLZ.git
cd RLZ
mkdir build
cd build
cmake ..
make -j
```

### Running the project

After building the project, an executable named rlz will be created in the build directory. Run it with:
```
./rlz -r [reference file] -s [file to compress] [options] 
```

### Example

In this section, we will go through a small example. In the data directory, we have provided an example reference and target sequence file that were derived from the English text in the [Pizza&Chili Corpus](https://pizzachili.dcc.uchile.cl/texts/nlang/).

1. To compress the sequence file, run the following command from the build directory

```
./rlz -r ../data/english_ref.txt -s ../data/english_seq.txt
```

This command will produce the following files in the data directory: `english_ref.txt.sdsl`, `english_seq.txt.sdsl`, and `english_seq.txt.rlz`. The .sdsl files are the bit vectors of the reference and sequence files that are written to file. The .rlz file contains the RLZ parse.

> [!NOTE]
> Multithreading is supported in the compression step with the -t [num. of threads] option which can significantly make this step faster. However, the RLZ parse is slightly different since we cannot identify phrases that potentially span where the file was split. Potentially might have an additional thread number of parse entries.

> [!NOTE]
> The compression ratio is quite high in this example. That is due to the example being small and not optimizing how the parse is written to file.

2. To decompress the file, run the following command

```
./rlz -r ../data/english_ref.txt -s ../data/english_seq.txt -d
```
This command should produce a file called `english_seq.txt.out` in the data directory. This is the decompressed sequence file.

3. Check to see if the file decompressed correctly
```
diff ../data/english_seq.txt ../data/english_seq.txt.out
```

There should be no output from this command if compressed and decompressed correctly. 

> [!NOTE]
> To get more information from the tool. Run the command with --verbose flag.

### License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/rvarki/rlz/blob/main/LICENSE) file for details

## Dependencies
- [SDSL](https://github.com/simongog/sdsl-lite)
- [CLI11](https://github.com/CLIUtils/CLI11)
- [spdlog](https://github.com/gabime/spdlog)

## Acknowledgements

- [Dhruv R. Makwana](https://github.com/Dhruv-mak) [Helped develop the code]

- [S. Kuruppu, S. J. Puglisi and J. Zobel, Relative Lempel-Ziv Compression of Genomes for Large-Scale Storage and Retrieval](http://dx.doi.org/10.1007/978-3-642-16321-0_20). Proc. 17th International Symposium on String Processing and Information Retrieval (SPIRE 2010) Lecture Notes in Computer Science, Volume 6393, (2010) pp. 201-206. 