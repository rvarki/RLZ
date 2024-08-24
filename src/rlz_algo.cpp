#include "rlz_algo.h"
#include "fm_wrapper.h"
#include <fstream>
#include <cstdlib>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/int_vector.hpp>
#include <vector>
#include <tuple>
#include <cstdio>
#include <stack>
#include <omp.h>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"


/**
* @brief Constuctor of RLZ class.
*
* Assigns the ref_file and seq_file variables with the file paths
*
* @param[in] ref_file [string] Path to reference file 
* @param[in] seq_file [string] Path to sequence file
*/

RLZ::RLZ(const std::string ref_file, const std::string seq_file): ref_file(ref_file), seq_file(seq_file){}

/**
* @brief Destructor of RLZ class.
*
* Currently does nothing.
*
*/

RLZ::~RLZ(){}

/**
* @brief Loads both the reference and sequence content into sdsl bit vectors.
*
* Wrapper function that calls the load_file_to_bit_vector function 
* on both the reference and sequence file to load them into bit vectors. 
*
* @return void
*/

void RLZ::load_bit_vectors()
{
    spdlog::stopwatch sw_ref;
    spdlog::debug("Starting to store the reference file as a bit vector");
    load_file_to_bit_vector(ref_file, ref_bit_array);
    auto sw_ref_elapsed = sw_ref.elapsed();
    spdlog::debug("Loaded file in {:.3} seconds", sw_ref_elapsed.count());
    spdlog::stopwatch sw_seq;
    spdlog::debug("Starting to store the sequence file as a bit vector");
    load_file_to_bit_vector(seq_file, seq_bit_array);
    auto sw_seq_elapsed = sw_ref.elapsed();
    spdlog::debug("Loaded file in {:.3} seconds", sw_seq_elapsed.count());
}


/**
* @brief Loads the file content into a bit vector.
*
* Loads the file content directly into a sdsl bit vector. Opens the input file in binary mode
* and moves pointer at end of file to get file size quickly. We then resize the bit vector
* to be large enough to hold the file content in bits. Read file byte by byte and store
* in bit vector.   
*
* @param[in] input_file [string] Path to either the reference or sequence file 
* @param[in] bit_array [sdsl::bit_vector] The corresponding bit array to store the file
* @return void
*/

void RLZ::load_file_to_bit_vector(const std::string& input_file, sdsl::bit_vector& bit_array)
{
    spdlog::stopwatch sw_convert;
    spdlog::debug("Reading file and creating bit array");

    // std::ios::ate moves cursor to end of file
    std::ifstream file(input_file, std::ios::binary | std::ios::ate);
    if (!file) {
        spdlog::error("Error opening {}", input_file);
        std::exit(EXIT_FAILURE);
    }

    // Get the file size in bytes
    std::streamsize file_size = file.tellg();
    // std::ios::beg moves cursor to beginning
    file.seekg(0, std::ios::beg);

    // Resize the bit array to hold the number of bits required
    bit_array.resize(file_size * 8);

    // Read the file and populate the bit vector
    char byte;
    std::size_t bit_index = 0;
    while (file.get(byte)) {
        // Starts from most sig bit of byte becoming least significant bit and then mask other bits and store in bit array.
        // This stores the bits of the byte in order.
        for (int i = 7; i >= 0; --i) {
            bit_array[bit_index++] = (byte >> i) & 1;
        }
    }

    file.close();
    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished creating bit array in {:.3} seconds", sw_convert_elapsed.count());

    spdlog::stopwatch sw_save;
    spdlog::debug("Saving bit array sdsl object to file");
    // Save the bit vector to a file
    sdsl::store_to_file(bit_array, input_file + ".sdsl");
    auto sw_save_elapsed = sw_save.elapsed();
    spdlog::debug("Finished saving in {:.3} seconds", sw_save_elapsed.count());
}

/**
* @brief Parses the sequence file in relation to the reference file
*
* This function does the RLZ parsing of the sequence file. It currently works at the
* "psuedo" bit level. To clarify, it currently processes the string representation of the bits of the 
* sequence file. Working at bit level allows us to compress all types of files.
*
* RLZ algorithm tries to greedily find the longest substring match within the reference starting from the
* first character in the sequence character such that the RLZ parse in the end contains (pos, len) pairs
* in relation to the reference such that the sequence can be reconstructed from only the RLZ parse and the reference
* file. It is a O(n) algorithm. The size is the reference file + the RLZ parse.
*
* The algorithm implemented here is as follows.
* 1. Starting from the last bit of the sequence file or sequence file chunk, check if bit matches the reference
* (via backwards match with FM-index) 
* 2a. If match, check if next bit also matches (ex. 001. I know that 1 matches then check if 01 matches etc...)
* 2b. If match and end of sequence file or sequence file chunk, push current (pos,len) pair to parse stack
* 2c. If mismatch, push (prev pos, len - 1) to parse stack. Reset search from bit that caused mismatch.
*
* Push to parse stack since we process the string in reverse. Popping from stack gives correct order.
*
* @param [in] fm_index [sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>] the fm-index of the reference
* @param [in] fm_support [FM_Wrapper] Utility object that allows us to do search and locate queries with fm-index.
* @param [in] seq_parse_stack_vec [std::vector<std::stack<std::tuple<uint64_t, size_t>>>] empty RLZ parse stacks equal to number of threads
* @param [in] num_bits_to_process [size_t] the number of bits that should be processed. Useful for the OpenMP parallelization.
* @param [in] loop_iter [size_t] the loop iteration. Useful for OpenMP and making sure we are thread-safe.
* @param [in] num_threads [size_t] the total number of threads allocated.
*
* @return void
*/

void RLZ::parse(const sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index,
        FM_Wrapper& fm_support,
        const sdsl::bit_vector& seq_bit_array,
        std::vector<std::stack<std::tuple<uint64_t, size_t>>>& seq_parse_stack_vec,
        size_t num_bits_to_process,
        size_t loop_iter,
        size_t num_threads)
{
    std::string pattern = "";
    size_t prev_left = 0;
    size_t prev_right = fm_index.bwt.size();
    size_t next_left = 0;
    size_t next_right = fm_index.bwt.size();
    long long int seq_size = static_cast<long long int>(seq_bit_array.size());

    long long int start_loc = (seq_size - 1) - (loop_iter * num_bits_to_process);
    long long int end_loc;

    // Last chunk so process remaining bits
    if (loop_iter == num_threads - 1)
        end_loc = -1;
    // Process bits up to next chunk
    else
        end_loc = (seq_size - 1) - ((loop_iter + 1) * num_bits_to_process);

    // Process the file in reverse for backwards matching with FM-index.
    for (long long int i = start_loc; i > end_loc; i--) 
    {
        char next_char = seq_bit_array[i] ? '1' : '0';

        pattern = next_char + pattern;

        std::tuple<size_t,size_t> previous_ranges = std::make_tuple(prev_left, prev_right);
        std::tuple<size_t,size_t> next_ranges = fm_support.backward_match(fm_index, previous_ranges, next_char);
        next_left = std::get<0>(next_ranges);
        next_right = std::get<1>(next_ranges);

        // If same then that means no perfect match so we reset.
        if (next_left == next_right){
            seq_parse_stack_vec[loop_iter].push(std::make_tuple(fm_support.get_suffix_array_value(fm_index, prev_left), pattern.size()-1));
            prev_left = 0;
            prev_right = fm_index.bwt.size()-1;
            next_left = 0;
            next_right = fm_index.bwt.size()-1;
            pattern = "";
            ++i;
        }
        // If at the end we are still in a perfect match, we save what we have. 
        else if (i == end_loc + 1)
        {
            seq_parse_stack_vec[loop_iter].push(std::make_tuple(fm_support.get_suffix_array_value(fm_index, next_left), pattern.size()));
        }
        // Currently in a perfect match
        else{
            prev_left = next_left;
            prev_right = next_right;
        }
    }
}

/**
* @brief Compresses the sequence file in relation to the reference file.
*
* Creates a FM-index from the reference bit array which we query using the sequence bit array.
* We first convert the reference bit array into its string representation so that we can create the FM-index.
* We query the index one bit at time from the sequence bit array. When the sequence bit does not have a match, 
* we add the last matching ref position of the sequence and the length of the match to the parse. Then we 
* restart the match at the last mismatch position. The parse is stored on a stack due to processing the bits 
* in reverese (backwards match with FM-index). Popping from stack gives correct order. The parse is ultimately
* stored in a vector in the correct order. The parse at the end is serialized to a file.
*
* @param [in] threads [int] The number of threads provided by the user.
*
* @return void
*
* @warning Providing multiple threads changes the output of the RLZ parse slightly. 
* Might create two phrases at chunk boundaries if phrase spans chunk boundary. For proper RLZ parse should run with 1 thread.
* 
* @note Supposedly cannot create a FM-index directly from bit array.
* Have to first convert into the reference bits into their string representation and then create the FM-index.
* Likely a bottleneck in the code as have to store a bit as a byte. [check if there is a way to build bit level FM-index]
*
* @note Might be more efficient to serialize the stack then create the vector [implementation detail]
*/

void RLZ::compress(int threads)
{
    sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024> fm_index;
    std::string binary_reference_text;

    // Convert the reference bit array into its string representation
    for (size_t i = 0; i < ref_bit_array.size(); ++i) {
        binary_reference_text += (ref_bit_array[i] ? '1' : '0');
    }

    // Creates the FM-index
    construct_im(fm_index, binary_reference_text, 1);

    FM_Wrapper fm_support;
    fm_support.rank0 = fm_index.bwt.rank(fm_index.size(), '0');
    fm_support.rank1 = fm_index.bwt.rank(fm_index.size(), '1');

    std::vector<std::stack<std::tuple<uint64_t, size_t>>> seq_parse_stack_vec(threads);
    size_t num_bits_to_process = seq_bit_array.size() / threads;  // Integer division

    // Comment (Testing only)
    // bits_to_str(seq_bit_array, ".orig.bits");

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < threads; i++)
    {
        parse(fm_index, fm_support, seq_bit_array, seq_parse_stack_vec, num_bits_to_process, i, threads);
    }

    // Pop from stack and store in seq_parse vector
    size_t bits_stored = 0;
    std::vector<std::tuple<uint64_t, size_t>> seq_parse;
    // Have to process the parse stacks in reverse order since the first stack contains the parse of the end of the sequence.
    for (int i = threads - 1; i >= 0; i--)
    {
        while (!seq_parse_stack_vec[i].empty())
        {
            bits_stored += std::get<1>(seq_parse_stack_vec[i].top());
            seq_parse.emplace_back(seq_parse_stack_vec[i].top());
            seq_parse_stack_vec[i].pop();
        }
    }

    spdlog::debug("The sequence was encoded in {} bits", seq_bit_array.size());
    spdlog::debug("The rlz parse encodes for {} bits", bits_stored);

    serialize(seq_parse);

    // Comment (Testing only)
    // print_serialize(seq_parse);
}

/**
* @brief Serializes the parse of the sequence file
*
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We serialize the parse vector into binary file called seq_file_name.rlz
* 
* @param[in] seq_parse [std::vector<std::tuple<uint64_t, size_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
*
* @return void
*/

void RLZ::serialize(const std::vector<std::tuple<uint64_t, size_t>>& seq_parse)
{
    std::ofstream ofs(seq_file + ".rlz", std::ios::binary);
    if (!ofs) {
        spdlog::error("Error opening {}", seq_file + ".rlz");
        std::exit(EXIT_FAILURE);
    }
    size_t size = seq_parse.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    ofs.write(reinterpret_cast<const char*>(seq_parse.data()), size * sizeof(std::tuple<uint64_t, size_t>));
    ofs.close();
}


/**
* @brief Deserializes the parse of the sequence file
*
* Decompress seq_file_name.rlz into tuple vector <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... > . 
* Return the vector.
*
* @return Return the vector.
*/

std::vector<std::tuple<uint64_t, size_t>> RLZ::deserialize()
{
    std::ifstream ifs(seq_file + ".rlz", std::ios::binary);
    if (!ifs) {
        spdlog::error("Error opening {}", seq_file + ".rlz");
        std::exit(EXIT_FAILURE);
    }
    size_t size;
    std::vector<std::tuple<uint64_t, size_t>> seq_parse;

    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    seq_parse.resize(size);
    ifs.read(reinterpret_cast<char*>(seq_parse.data()), size * sizeof(std::tuple<uint64_t, size_t>));
    ifs.close();

    return seq_parse;
}

/**
* @brief Decompresses the sequence parse back to the original sequence file
*
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We read the ref position and length from each tuple in the sequence parse and get the corresponding bits from the reference bit array
* We then convert the bits into bytes and write the ASCII equivalent to file which sbould be the original sequence file.
*
* @warning might have to change int to long long int depending on size
*/

void RLZ::decompress()
{
    std::vector<std::tuple<uint64_t, size_t>> seq_parse = deserialize();
    sdsl::load_from_file(ref_bit_array, ref_file + ".sdsl");
    
    int bit_size = 0;
    for (const auto& [pos, len] : seq_parse){
        bit_size += len;
    }

    spdlog::debug("The compessed sequence file had {} bits", bit_size);

    // Resize the array to be equal to the number of bits to be stored
    seq_bit_array.resize(bit_size);

    // Get the sequence bits from the parse + reference bits
    int prev_pos = 0;
    for (const auto& [pos, len] : seq_parse){
        int curr_pos = prev_pos + len;
        int len_count = 0;
        for (int i = prev_pos; i < curr_pos; ++i){
            seq_bit_array[i] = ref_bit_array[pos + len_count];
            len_count++;
        }
        prev_pos = curr_pos;
    }

    // Comment (Testing only)
    // bits_to_str(seq_bit_array, ".decompress.bits");

    std::ofstream output_file(seq_file + ".out");
    if (!output_file) {
        spdlog::error("Error opening {}", seq_file + ".out");
        std::exit(EXIT_FAILURE);
    }

    // Convert the sequence bits back into string
    std::string uncompressed_seq;

    uint8_t byte = 0;
    size_t bit_count = 0;

    // Write the sequence bits back into characters
    for (size_t i = 0; i < seq_bit_array.size(); ++i) {
        byte <<= 1; // Shift byte left by 1 bit
        byte |= seq_bit_array[i]; // Add the current bit to the byte
        ++bit_count;

        if (bit_count == 8) { // If 8 bits have been processed
            uncompressed_seq += static_cast<char>(byte); // Convert byte to char and append to result
            byte = 0; // Reset byte
            bit_count = 0; // Reset bit count
        }
    }

    output_file << uncompressed_seq;
    output_file.close();
}


/**
* @brief Converts the bits to its corresponding character representation
*
* Prints out debug info about how many 0 and 1s are in the bit array.
* Writes the string representation of the bits to a file.
* Testing purposes only.
* 
* @param [in] bit_array [sdsl::bit_vector] the bit array of interest
* @param [in] ext [std::string] the extension of the file
*
*
*/

void RLZ::bits_to_str(sdsl::bit_vector bit_array, std::string ext)
{
    std::string bitstr;
    for (size_t i = 0; i < bit_array.size(); ++i) {
        bitstr += (bit_array[i] ? '1' : '0');
    }

    sdsl::rank_support_v<1> b_rank1(&bit_array);
    sdsl::rank_support_v<0> b_rank0(&bit_array);
    spdlog::debug("^^^^^^^^^^^^^^^^^^^^^^^^^");
    spdlog::debug("Number of 1s: {}", b_rank1(bit_array.size()));
    spdlog::debug("Number of 0s: {}", b_rank0(bit_array.size()));
    spdlog::debug("^^^^^^^^^^^^^^^^^^^^^^^^^");

    std::ofstream ofs(seq_file + ext);
    if (!ofs) {
        spdlog::error("Error opening {}", seq_file + ext);
        std::exit(EXIT_FAILURE);
    }
    ofs << bitstr;
    ofs.close();
}

/**
* @brief Write the non-binary serialization of the sequence parse to a file.
*
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We write the non-binary serialization to a file. Testing purposes only.
* 
* @param[in] seq_parse [std::vector<std::tuple<uint64_t, size_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
*
* @return void
*/

void RLZ::print_serialize(const std::vector<std::tuple<uint64_t, size_t>>& seq_parse)
{
    std::ofstream ofs(seq_file + ".readable.rlz");
    if (!ofs) {
        spdlog::error("Error opening {}", seq_file + ".readable.rlz");
        std::exit(EXIT_FAILURE);
    }

    for (const auto& [pos, len] : seq_parse){
        ofs << "Position: " << pos << "  Length: " << len << std::endl;
    }

    ofs.close();
}


/**
* @brief Cleans some intermediate file
*
* Removes intermediate files for faster testing
*
* @return void
*/

void RLZ::clean()
{
    std::remove((ref_file + ".sdsl").c_str());
    spdlog::info("Removed {} if it existed", ref_file + ".sdsl");
    std::remove((seq_file + ".sdsl").c_str());
    spdlog::info("Removed {} if it existed", seq_file + ".sdsl");
    std::remove((seq_file + ".rlz").c_str());
    spdlog::info("Removed {} if it existed", seq_file + ".rlz");
    std::remove((seq_file + ".readable.rlz").c_str());
    spdlog::info("Removed {} if it existed", seq_file + ".readable.rlz");
    std::remove((seq_file + ".orig.bits").c_str());
    spdlog::info("Removed {} if it existed", seq_file + ".orig.bits");
    std::remove((seq_file + ".decompress.bits").c_str());
    spdlog::info("Removed {} if it existed", seq_file + ".decompress.bits");
    std::remove((seq_file + ".out").c_str());
    spdlog::info("Removed {} if it existed", seq_file + ".out");
}


