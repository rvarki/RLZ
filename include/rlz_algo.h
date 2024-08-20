#ifndef RLZ_ALGO_H
#define RLZ_ALGO_H

#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <fstream>
#include <vector>
#include <tuple>

namespace RLZ_Algo {

    class RLZ {
        public:
            std::string ref_file;
            std::string seq_file;
            sdsl::bit_vector ref_bit_array;
            sdsl::bit_vector seq_bit_array;

            RLZ(const std::string ref_file, const std::string seq_file);
            ~RLZ();
            void compress();
            void decompress();
            void load_bit_vectors();
            //void clean();
            void save_as_binary_file(const std::string& infile);
    
        private:
            void load_file_to_bit_vector(const std::string& input_file, sdsl::bit_vector& bit_array);
            void bit_vector_resize(sdsl::bit_vector& bit_array, std::streamsize file_size);
            void serialize(const std::vector<std::tuple<uint64_t, size_t>>& seq_parse);
            void print_serialize(const std::vector<std::tuple<uint64_t, size_t>>& seq_parse);
            std::vector<std::tuple<uint64_t, size_t>> deserialize();
    };

    class FM_Wrapper {
        public:
            size_t rank0;
            size_t rank1;

            FM_Wrapper();
            ~FM_Wrapper();
            std::tuple<size_t, size_t> backward_match(sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index, 
                                                    std::tuple<size_t, size_t>& prev_backward_range,
                                                    char next_char);
            
            size_t get_suffix_array_value(sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index, size_t location);
    };
}

#endif  // RLZ_ALGO_H
