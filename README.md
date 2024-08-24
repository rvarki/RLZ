# Relative Lempel-Ziv (RLZ)

## Description

This software computes the Relative Lempel Ziv (RLZ) parse of the target file using a reference file. The software should work for any type of file whether that be FASTA files, English files, etc... However, currently the decompression expects that the files were originally ASCII (8 bit) encoded.

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

### License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/rvarki/rlz/blob/main/LICENSE) file for details

## Dependencies
- [CLI11](https://github.com/CLIUtils/CLI11)
- [SDSL](https://github.com/simongog/sdsl-lite)

## Acknowledgements

- [Dhruv R. Makwana](https://github.com/Dhruv-mak) [Helped develop the code]

- [S. Kuruppu, S. J. Puglisi and J. Zobel, Relative Lempel-Ziv Compression of Genomes for Large-Scale Storage and Retrieval](http://dx.doi.org/10.1007/978-3-642-16321-0_20). Proc. 17th International Symposium on String Processing and Information Retrieval (SPIRE 2010) Lecture Notes in Computer Science, Volume 6393, (2010) pp. 201-206. 