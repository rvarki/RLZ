#ifndef RLZ_ALGO_H
#define RLZ_ALGO_H

#include <sdsl/bit_vectors.hpp>
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
            //void save_as_binary_file(const std::string& infile);
    
        private:
            void load_file_to_bit_vector(const std::string& input_file, sdsl::bit_vector& bit_array);
            void bit_vector_resize(sdsl::bit_vector& bit_array, std::streamsize file_size);
            void serialize(const std::vector<std::tuple<uint64_t, int>>& seq_parse);
            std::vector<std::tuple<uint64_t, int>> deserialize();
    };
}

#endif  // RLZ_ALGO_H
