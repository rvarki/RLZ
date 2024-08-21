# RLZ

## Description

This software computes the Relative Lempel Ziv (RLZ) parse of the target file using a reference file. 

## Prerequisites

- [CMake](https://cmake.org/) 3.15 or higher.
- C++17-compatible compiler.

## Getting Started

### Building the Project

```
git clone https://github.com/Dhruv-mak/rlz
cd rlz
mkdir build
cd build
cmake ..
make
```

### Running the project

After building the project, an executable named rlz will be created in the build directory. Run it with:
```
./rlz [options] [args]
```

### License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Dhruv-mak/rlz/blob/main/LICENSE) file for details

## Dependencies
- [CLI11](https://github.com/CLIUtils/CLI11)
- [SDSL](https://github.com/simongog/sdsl-lite)

### Acknowledgements

- [S. Kuruppu, S. J. Puglisi and J. Zobel, Relative Lempel-Ziv Compression of Genomes for Large-Scale Storage and Retrieval](http://dx.doi.org/10.1007/978-3-642-16321-0_20). Proc. 17th International Symposium on String Processing and Information Retrieval (SPIRE 2010) Lecture Notes in Computer Science, Volume 6393, (2010) pp. 201-206. 