#include "rlz_algo.h"
#include <fstream>
#include <cstdlib>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/int_vector.hpp>
#include <vector>
#include <tuple>
#include <cstdio>
#include <stack>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

namespace RLZ_Algo {

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
    * @brief Loads both the reference and sequence content into bit vectors.
    *
    * Calls load_file_to_bit_vector
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
    * @brief Loads the file content into bit vectors.
    *
    * Loads the file content directly into sdsl bit vectors. Opens the input file in binary mode
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
            std::cerr << "Error opening file!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // Get the file size
        std::streamsize file_size = file.tellg();
        // std::ios::beg moves cursor to beginning
        file.seekg(0, std::ios::beg);

        // Effectively resizes the bit array 
        bit_vector_resize(bit_array, file_size);

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
        sdsl::store_to_file(bit_array, input_file + ".bit_array.sdsl");
        auto sw_save_elapsed = sw_save.elapsed();
        spdlog::debug("Finished saving in {:.3} seconds", sw_save_elapsed.count());

    }

    /**
    * @brief "Resizes" the bit vector.
    *
    * sdsl bit vectors are immutable. Therefore resizing means creating a new bit vector and
    * overwriting the previous bit member which is a member of our RLZ class. Should be safe.
    *
    * @param[in] bit_array [sdsl::bit_vector] Path to the bit vector to resize
    * @param[in] bit_array [std::streamsize] The size (bytes) of the file which the bit vector should be resized. 
    * @return void
    * @note the bit array size passed should be in bits so multiple the bytes by 8.
    * @warning Assumes that the original bit vector is empty.
    */

    void RLZ::bit_vector_resize(sdsl::bit_vector& bit_array, std::streamsize file_size)
    {
        sdsl::bit_vector new_bit_array(file_size * 8);
        // This is safe.
        bit_array = std::move(new_bit_array);
    }


    /**
    * @brief Compresses the sequence file in relation to the reference file.
    *
    * Creates a FM-index from the reference bit array which we query using the query bit array.
    * We first convert the reference bit array into a string so that we can create the FM-index.
    * We query the index with the sequence. When the sequence does not have a match, we add the
    * last matching ref position of the sequence and the length of the match to the parse. Then we 
    * restart the match at the last mismatch position. The parse gets serialized afterwards. 
    * 
    * @warning Supposedly cannot create a FM-index directly from bit array.
    * Have to first convert into string and then create the FM-index.
    * Likely a bottleneck in the code. 
    *
    * @warning The FM-index of sdsl has the limitation that it does not give any way to find position of 
    * last char match of pattern in the event that the pattern does not perfectly match a substring in the index.
    * This means that a lot of redundant work is going to happen since we have to keep rechecking previous matches.
    * Have to write a wrapper around sdsl FM-index or maybe look at r-index implementation and see if that does it, idk. [maybe addressed]
    *
    * @warning might be more efficient to serialize the stack if this work than create the vector from stack and then serialize
    */

    void RLZ::compress()
    {
        sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024> fm_index;
        std::string binary_reference_text;

        for (size_t i = 0; i < ref_bit_array.size(); ++i) {
            binary_reference_text += (ref_bit_array[i] ? '1' : '0');
        }
        // Creates the FM-index
        construct_im(fm_index, binary_reference_text, 1);

        RLZ_Algo::FM_Wrapper fm_support;
        fm_support.rank0 = fm_index.bwt.rank(fm_index.size(), '0');
        fm_support.rank1 = fm_index.bwt.rank(fm_index.size(), '1');

        std::string pattern = "";
        size_t prev_left = 0;
        size_t prev_right = fm_index.bwt.size();
        size_t next_left = 0;
        size_t next_right = fm_index.bwt.size();
        
        std::vector<std::tuple<uint64_t, size_t>> seq_parse;
        // Since we process the phrases backwards we store in stack initially
        std::stack<std::tuple<uint64_t, size_t>> seq_parse_stack;

        // Static cast size_t into a signed integer (hopefully large enough)
        long long int seq_size = static_cast<long long int>(seq_bit_array.size());

        // How to process the bits in reverse for backwards matching
        for (long long int i = seq_size - 1 ; i >= 0; i--) 
        {
            char next_char = seq_bit_array[i] ? '1' : '0';

            pattern = next_char + pattern;

            std::tuple<size_t,size_t> previous_ranges = std::make_tuple(prev_left, prev_right);
            std::tuple<size_t,size_t> next_ranges = fm_support.backward_match(fm_index, previous_ranges, next_char);
            next_left = std::get<0>(next_ranges);
            next_right = std::get<1>(next_ranges);

            // If same then that means no perfect match so we reset.
            if (next_left == next_right){
                seq_parse_stack.push(std::make_tuple(fm_support.get_suffix_array_value(fm_index, prev_left), pattern.size()-1));
                prev_left = 0;
                prev_right = fm_index.bwt.size()-1;
                next_left = 0;
                next_right = fm_index.bwt.size()-1;
                pattern = "";
                ++i;
            }
            // If at the end we are still in a perfect match, we save what we have. 
            else if (i == 0)
            {
                seq_parse_stack.push(std::make_tuple(fm_support.get_suffix_array_value(fm_index, next_left), pattern.size()));
            }
            // Currently in a perfect match
            else{
                prev_left = next_left;
                prev_right = next_right;
            }
        }

        // Pop from stack and store in seq_parse vector
        int count = 0;
        while (!seq_parse_stack.empty()) {
            count += std::get<1>(seq_parse_stack.top());
            seq_parse.emplace_back(seq_parse_stack.top());
            seq_parse_stack.pop();
        }

        serialize(seq_parse);
        print_serialize(seq_parse);
    }


    /**
    * @brief Serializes the parse of the sequence file
    *
    * The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
    * We serialize the parse vector into binary file called seq_file_name.rlz
    * 
    * @param[in] seq_parse [std::vector<std::tuple<uint64_t, size_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
    *
    */

    void RLZ::serialize(const std::vector<std::tuple<uint64_t, size_t>>& seq_parse)
    {
        std::ofstream ofs(seq_file + ".rlz", std::ios::binary);
        if (!ofs) {
            std::cerr << "Error opening file!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        size_t size = seq_parse.size();
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
        ofs.write(reinterpret_cast<const char*>(seq_parse.data()), size * sizeof(std::tuple<uint64_t, size_t>));
        ofs.close();
    }


    /**
    * @brief Write the non-binary serialization of the sequence parse to a file.
    *
    * The sequence parse contains tuples (binary ref pos, size) that can reconstruct the sequence file given the reference.
    * We write the non-binary serialization to a file. Testing purposes.
    * 
    * @param[in] seq_parse [std::vector<std::tuple<uint64_t, size_t>>] The parse of the seq <(binary ref pos,len),(binary ref pos,len),(binary ref pos,len)... >
    *
    */

    void RLZ::print_serialize(const std::vector<std::tuple<uint64_t, size_t>>& seq_parse)
    {
        std::ofstream ofs(seq_file + ".readable.rlz");
        if (!ofs) {
            std::cerr << "Error opening file!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (const auto& [pos, len] : seq_parse){
            ofs << "Position: " << pos << "  Length: " << len << std::endl;
        }

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
            std::cerr << "Error opening file!" << std::endl;
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
        sdsl::load_from_file(ref_bit_array, ref_file + ".bit_array.sdsl");
        
        int bit_size = 0;
        for (const auto& [pos, len] : seq_parse){
            bit_size += len;
        }

        // Should be wholly divisible by 8
        bit_vector_resize(seq_bit_array, bit_size/8);

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

        std::ofstream output_file(seq_file + ".out");
        if (!output_file) {
            std::cerr << "Error opening file for writing!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::ofstream output_file_bin(seq_file + ".out.bin", std::ios::binary);

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
                output_file_bin.write(reinterpret_cast<const char*>(&byte), sizeof(byte));
                uncompressed_seq += static_cast<char>(byte); // Convert byte to char and append to result
                byte = 0; // Reset byte
                bit_count = 0; // Reset bit count
            }
        }

        output_file << uncompressed_seq;
        output_file.close();
        output_file_bin.close();
    }

    /**
    * @brief Saves regular file as binary
    *
    * Saves input file as binary file. Used for testing purposes
    *
    */

    // void RLZ::save_as_binary_file(const std::string& infile)
    // {
    //     std::ifstream input_file(infile, std::ios::binary);
    //     std::ofstream output_file(infile + ".bin", std::ios::binary);
    //     output_file << input_file.rdbuf();
    //     input_file.close();
    //     output_file.close();
    // }


    /**
    * @brief Cleans some intermediate file
    *
    * Removes intermediate files during compression for faster testing
    *
    * @warning not safe but yolo
    */

    // void RLZ::clean()
    // {
    //     std::remove((ref_file + ".bit_array.sdsl").c_str());
    //     std::remove((seq_file + ".bit_array.sdsl").c_str());
    //     std::remove((seq_file + ".rlz").c_str());
    // }


    /**
    * @brief FM Wrapper Constructor
    *
    * Does nothing 
    *
    */

    FM_Wrapper::FM_Wrapper(){}

     /**
    * @brief FM Wrapper Destructor
    *
    * Does nothing 
    *
    */

    FM_Wrapper::~FM_Wrapper(){}

    /**
    * @brief Extend backward match with FM-index
    *
    * A wrapper around sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024> fm_index.
    * Calculates the new backwards search range after trying to match the next char (backwards)
    * Can continue from previous character match so do not have to keep redoing previously done backwards matches.
    *
    * @param[in] fm_index [sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>] the fm_index queried
    * @param[in] prev_backward_range [sdsl::range] the previous backwards search range.
    * @param[in] next_char [char] the next character to match 
    * 
    * @return the backwards search range of next_char 
    */

    std::tuple<size_t, size_t> FM_Wrapper::backward_match(sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index, 
                                                        std::tuple<size_t, size_t>& prev_backward_range,
                                                        char next_char)
    {
        std::tuple<size_t, size_t> next_char_backward_range;
        std::size_t next_left = fm_index.bwt.rank(std::get<0>(prev_backward_range), next_char);
        std::size_t next_right = fm_index.bwt.rank(std::get<1>(prev_backward_range), next_char);

        // There is a special (smaller) character appended so have to add 1 to the offset for both
        // For bit 1 have to offset the range by the number of 0 bits. 
        if (next_char == '0'){
            next_left += 1;
            next_right += 1;
        }
        if (next_char == '1'){
            next_left += this->rank0 + 1;
            next_right += this->rank0 + 1;
        }

        return std::make_tuple(next_left,next_right);
    }

    /**
    * @brief Get corresponding location in reference via suffix array
    *
    * A wrapper around sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024> fm_index.
    * When the backwards search range becomes empty. Means there is no perfect match with the current pattern that is being processed.
    * With the prior non-empty range we find one location of the previous perfect pattern match. This should always be next_left of backward_match
    *
    * @param[in] fm_index [sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>] the fm_index queried
    * @param[in] location [size_t] I think should always be next_left of backward_match.
    * 
    * @return the suffix array index
    */

    size_t FM_Wrapper::get_suffix_array_value(sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index,
                                            size_t location)
    {
        return fm_index[location];
    }

}
