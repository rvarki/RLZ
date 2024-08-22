#ifndef FM_WRAPPER_H
#define FM_WRAPPER_H

#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <fstream>
#include <vector>
#include <tuple>

class FM_Wrapper 
{
    public:
        size_t rank0;
        size_t rank1;

        FM_Wrapper();
        ~FM_Wrapper();
        std::tuple<size_t, size_t> backward_match(const sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index, 
                                                const std::tuple<size_t, size_t>& prev_backward_range,
                                                const char next_char);
        
        size_t get_suffix_array_value(const sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>& fm_index, 
                                        const size_t location);
};


#endif // FM_WRAPPER_H