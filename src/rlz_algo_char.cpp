#include "rlz_algo_char.h"
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
#include <map>
#include <omp.h>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"


/**
* @brief Constuctor of RLZ_CHAR class.
*
* Assigns the ref_file and seq_file variables with the file paths
*
* @param[in] ref_file [string] Path to reference file 
* @param[in] seq_file [string] Path to sequence file
*/

RLZ_CHAR::RLZ_CHAR(const std::string ref_file, const std::string seq_file): ref_file(ref_file), seq_file(seq_file){}

/**
* @brief Destructor of RLZ_CHAR class.
*
* Currently does nothing.
*
*/

RLZ_CHAR::~RLZ_CHAR(){}

/**
* @brief Loads the file content into a string.
*
* Loads the file content directly into a string. Opens the input file in binary mode
* and moves pointer at end of file to get file size quickly. We then resize the string
* to be large enough to hold the file content in bytes.   
*
* @param[in] input_file [string] Path to either the reference or sequence file 
* @param[in] content [string] Where the string representation of the file is located.
* @return void
*/

void RLZ_CHAR::load_file_to_string(const std::string& input_file, std::string& content)
{
    spdlog::stopwatch sw_convert;
    spdlog::debug("Storing file as string");

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

    // Resize the content field
    content.resize(file_size);

    // Load the file as string
    if (file.read(&content[0], file_size)) {
        spdlog::info("File read successfully.");
    } else {
        spdlog::error("Error reading file: {}", input_file);
        std::exit(EXIT_FAILURE);
    }

    file.close();
    auto sw_convert_elapsed = sw_convert.elapsed();
    spdlog::debug("Finished storing file in {:.3} seconds", sw_convert_elapsed.count());

    // spdlog::stopwatch sw_save;
    // spdlog::debug("Saving bit array sdsl object to file");
    //// Save the bit vector to a file
    // sdsl::store_to_file(bit_array, input_file + ".sdsl");
    // auto sw_save_elapsed = sw_save.elapsed();
    // spdlog::debug("Finished saving in {:.3} seconds", sw_save_elapsed.count());
}

/**
* @brief Parses the sequence file in relation to the reference file
*
* This function does the RLZ_CHAR parsing of the sequence file. It parses relative to the original file content of the reference.
* Fails if there is a character in the sequence that is not present in the reference.
*
* RLZ algorithm tries to greedily find the longest substring match within the reference starting from the
* first character in the sequence character such that the RLZ parse in the end contains (pos, len) pairs
* in relation to the reference such that the sequence can be reconstructed from only the RLZ parse and the reference
* file. It is a O(n) algorithm. The size is the reference file + the RLZ parse.
*
* The algorithm implemented here is as follows.
* 1. Starting from the last char of the sequence file or sequence file chunk, check if char matches the reference
* (via backwards match with FM-index) 
* 2a. If match, check if next char also matches (ex. aab. I know that b matches then check if ab matches etc...)
* 2b. If match and end of sequence file or sequence file chunk, push current (pos,len) pair to parse stack
* 2c. If mismatch, push (prev pos, len - 1) to parse stack. Reset search from bit that caused mismatch.
*
* Push to parse stack since we process the string in reverse. Popping from stack gives correct order.
*
* @param [in] fm_index [sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>] the fm-index of the reference
* @param [in] fm_support [FM_Wrapper] Utility object that allows us to do search and locate queries with fm-index.
* @param [in] seq_parse_stack_vec [std::vector<std::stack<std::tuple<uint64_t, uint64_t>>>] empty RLZ_CHAR parse stacks equal to number of threads
* @param [in] num_char_to_process [size_t] the number of chars that should be processed. Useful for the OpenMP parallelization.
* @param [in] loop_iter [size_t] the loop iteration. Useful for OpenMP and making sure we are thread-safe.
* @param [in] num_threads [size_t] the total number of threads allocated.
*
* @return void
*/

void RLZ_CHAR::parse(const sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index,
        FM_Wrapper& fm_support,
        const std::map<char, uint64_t>& occs, 
        const std::string& seq_content,
        std::vector<std::stack<std::tuple<uint64_t, uint64_t>>>& seq_parse_stack_vec,
        size_t num_char_to_process,
        size_t loop_iter,
        size_t num_threads)
{
    std::string pattern = "";
    size_t prev_left = 0;
    size_t prev_right = fm_index.bwt.size();
    size_t next_left = 0;
    size_t next_right = fm_index.bwt.size();
    long long int seq_size = static_cast<long long int>(seq_content.size());

    long long int start_loc = (seq_size - 1) - (loop_iter * num_char_to_process);
    long long int end_loc;

    // Last chunk so process remaining bits
    if (loop_iter == num_threads - 1)
        end_loc = -1;
    // Process bits up to next chunk
    else
        end_loc = (seq_size - 1) - ((loop_iter + 1) * num_char_to_process);

    // Process the file in reverse for backwards matching with FM-index.
    for (long long int i = start_loc; i > end_loc; i--) 
    {
        char next_char = seq_content[i];

        pattern = next_char + pattern;

        std::tuple<size_t,size_t> previous_ranges = std::make_tuple(prev_left, prev_right);
        std::tuple<size_t,size_t> next_ranges = fm_support.backward_match(fm_index, occs, previous_ranges, next_char);
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
* @brief Calculates the occurances of each char in the provided text in lexicographical order
*
* @param [in] content [string] The string which we are deriving the occurances from
*
* @return void
*/

void RLZ_CHAR::calculate_occs(std::string content, std::map<char, uint64_t>& occs)
{
    // Sort the string lexicographically
    std::sort(content.begin(), content.end());
    uint64_t count = 0;
    char prev_char = '\0';

    for (size_t i = 0; i < content.size(); i++)
    {
        if (prev_char != content[i])
        {
            occs[content[i]] = count;
        }

        prev_char = content[i];
        count++;
    }
}


/**
* @brief Compresses the sequence file in relation to the reference file.
*
* Creates a FM-index from the reference string which we query using the sequence string.
* We first create the FM-index from the string representation of the reference.
* We query the index one char at time from the sequence string. When the sequence char does not have a match, 
* we add the last matching ref position of the sequence and the length of the match to the parse. Then we 
* restart the match at the last mismatch position. The parse is stored on a stack due to processing the chars 
* in reverse (backwards match with FM-index). Popping from stack gives correct order. The parse is ultimately
* stored in a vector in the correct order. The parse at the end is serialized to a file.
*
* @param [in] threads [int] The number of threads provided by the user.
*
* @return void
*
* @warning Providing multiple threads changes the output of the RLZ_CHAR parse slightly. 
* Might create two phrases at chunk boundaries if phrase spans chunk boundary. For proper RLZ_CHAR parse should run with 1 thread.
* 
* @warning Will fail if the sequence file contains a char not present in the reference file
*
* @note Supposedly cannot create a FM-index directly from bit array.
* Have to first convert into the reference bits into their string representation and then create the FM-index.
* Likely a bottleneck in the code as have to store a bit as a byte. [check if there is a way to build bit level FM-index]
*
* @note Might be more efficient to serialize the stack then create the vector [implementation detail]
*/

void RLZ_CHAR::compress(int threads)
{
    sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024> fm_index;
    
    // Creates the FM-index
    construct_im(fm_index, ref_content, 1);

    // Get the number of occurances of each char in lexicographical order
    std::map<char, uint64_t> occs;
    calculate_occs(ref_content, occs);

    FM_Wrapper fm_support;

    std::vector<std::stack<std::tuple<uint64_t, uint64_t>>> seq_parse_stack_vec(threads);
    size_t num_char_to_process = seq_content.size() / threads;  // Integer division

    // Comment (Testing only)
    // bits_to_str(seq_bit_array, ".orig.bits");

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < threads; i++)
    {
        parse(fm_index, fm_support, occs, seq_content, seq_parse_stack_vec, num_char_to_process, i, threads);
    }

    // Pop from stack and store in seq_parse vector
    size_t chars_stored = 0;
    std::vector<std::tuple<uint64_t, uint64_t>> seq_parse;
    // Have to process the parse stacks in reverse order since the first stack contains the parse of the end of the sequence.
    for (int i = threads - 1; i >= 0; i--)
    {
        while (!seq_parse_stack_vec[i].empty())
        {
            chars_stored += std::get<1>(seq_parse_stack_vec[i].top());
            seq_parse.emplace_back(seq_parse_stack_vec[i].top());
            seq_parse_stack_vec[i].pop();
        }
    }

    spdlog::debug("The sequence was encoded in {} chars", seq_content.size());
    spdlog::debug("The rlz parse encodes for {} chars", chars_stored);

    serialize(seq_parse);

    // Comment (Testing only)
    // print_serialize(seq_parse);
}

/**
* @brief Serializes the parse of the sequence file
*
* The sequence parse contains tuples (ref pos, size) that can reconstruct the sequence file given the reference.
* We serialize the parse vector into binary file called seq_file_name.rlz
*
* File content of the .rlz file
* (uint64_t byte: size num of pair) (uint64_t byte: size pos) (uint64_t byte: size len) (uint64_t byte: size pos) (uint64_t byte: size len) ...
*  
* @param[in] seq_parse [std::vector<std::tuple<uint64_t, uint64_t>>] The parse of the seq <(ref pos,len),(ref pos,len),(ref pos,len)... >
*
* @return void
*/

void RLZ_CHAR::serialize(const std::vector<std::tuple<uint64_t, uint64_t>>& seq_parse)
{
    std::ofstream ofs(seq_file + ".rlz", std::ios::binary);
    if (!ofs) {
        spdlog::error("Error opening {}", seq_file + ".rlz");
        std::exit(EXIT_FAILURE);
    }
    size_t size = seq_parse.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (size_t i = 0; i < size; i++)
    {
        ofs.write(reinterpret_cast<const char*>(&std::get<0>(seq_parse[i])), sizeof(uint64_t));
        ofs.write(reinterpret_cast<const char*>(&std::get<1>(seq_parse[i])), sizeof(uint64_t));
    }
    ofs.close();
}


/**
* @brief Deserializes the parse of the sequence file
*
* Decompress seq_file_name.rlz into tuple vector <(ref pos,len),(ref pos,len),(ref pos,len)... > . 
* Return the vector.
*
* @return Return the vector.
*/

std::vector<std::tuple<uint64_t, uint64_t>> RLZ_CHAR::deserialize()
{
    std::ifstream ifs(seq_file + ".rlz", std::ios::binary);
    if (!ifs) {
        spdlog::error("Error opening {}", seq_file + ".rlz");
        std::exit(EXIT_FAILURE);
    }
    uint64_t size;
    std::vector<std::tuple<uint64_t, uint64_t>> seq_parse;

    ifs.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    seq_parse.reserve(size);
    std::tuple<uint64_t, uint64_t> elem;
    uint64_t val;

    for (size_t i = 0; i < 2 * size; i++)
    {
        ifs.read(reinterpret_cast<char*>(&val), sizeof(uint64_t));
        if (i % 2 == 0){
            std::get<0>(elem) = val;
        }
        else{
            std::get<1>(elem) = val;
            seq_parse.emplace_back(elem);
        }
    }
    ifs.close();

    return seq_parse;
}

/**
* @brief Decompresses the sequence parse back to the original sequence file
*
* The sequence parse contains tuples (ref pos, size) that can reconstruct the sequence file given the reference.
* We read the ref position and length from each tuple in the sequence parse and get the corresponding chars from the reference to reconstruct the sequence
*
* @warning might have to change int to long long int depending on size
*/

void RLZ_CHAR::decompress()
{
    std::vector<std::tuple<uint64_t, uint64_t>> seq_parse = deserialize();
    
    size_t char_size = 0;
    for (const auto& [pos, len] : seq_parse){
        char_size += len;
    }

    spdlog::debug("The compessed sequence file had {} chars", char_size);

    // Resize the array to be equal to the number of bits to be stored
    seq_content.resize(char_size);

    // Get the sequence chars from the parse + reference chars
    size_t prev_pos = 0;
    for (const auto& [pos, len] : seq_parse){
        size_t curr_pos = prev_pos + len;
        size_t len_count = 0;
        for (size_t i = prev_pos; i < curr_pos; ++i){
            seq_content[i] = ref_content[pos + len_count];
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

    output_file << seq_content;
    output_file.close();
}


/**
* @brief Write the non-binary serialization of the sequence parse to a file.
*
* The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
* We write the non-binary serialization to a file. Testing purposes only.
* 
* @param[in] seq_parse [std::vector<std::tuple<uint64_t, uint64_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
*
* @return void
*/

void RLZ_CHAR::print_serialize(const std::vector<std::tuple<uint64_t, uint64_t>>& seq_parse)
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

void RLZ_CHAR::clean()
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


