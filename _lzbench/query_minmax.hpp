#ifndef QUERY_MIN_MAX_HPP
#define QUERY_MIN_MAX_HPP

#include "eigen/Core"

#include "immintrin.h"
#include <type_traits>

#include "query_common.h"

#include <iostream> // TODO rm

template<class DataT, int OpE>
class OnlineBinaryOpRowmajor {
public:
    using dist_t = typename DataTypeTraits<DataT>::AccumulatorT;
    using Op = BinaryOp<DataT, OpE, DataT>;

    OnlineBinaryOpRowmajor(uint32_t nrows, uint32_t ncols):
        _nrows(nrows), _ncols(ncols), _is_dense(true)
    {
        reset();
    }

    OnlineBinaryOpRowmajor(uint32_t nrows, uint32_t ncols,
        const std::vector<uint16_t>& which_dims):
        _nrows(nrows), _ncols(ncols), _which_dims(which_dims),
        _is_dense(which_dims.size() == 0)
    {
        // if (IsDense) {
        //     printf("ERROR: can't specify subset of dims for Dense OnlineMean!");
        //     exit(1);
        // }
        reset();
    }

    void init(const DataT* window_start) {
        if (_is_dense) {
            for (uint32_t i = 0; i < _nrows; i++) {
                for (uint32_t j = 0; j < _ncols; j++) {
                    _stats[j] = Op{}(_stats[j], window_start[i * _ncols + j]);
                }
            }
        } else {
            for (uint32_t i = 0; i < _nrows; i++) {
                for (uint32_t j_idx = 0; j_idx < _which_dims.size(); j_idx++) {
                    auto j = _which_dims[j_idx];
                    _stats[j_idx] =
                        Op{}(_stats[j_idx], window_start[i * _ncols + j]);
                }
            }
        }
    }

    void update(const DataT* old_window_row, const DataT* new_window_row) {
        if (_is_dense) {
            for (uint32_t j = 0; j < _ncols; j++) {
                _stats[j] = Op{}(_stats[j], new_window_row[j]);
            }
        } else {
            for (uint32_t j_idx = 0; j_idx < _which_dims.size(); j_idx++) {
                auto j = _which_dims[j_idx];
                _stats[j_idx] = Op{}(_stats[j_idx], new_window_row[j]);
            }
        }
    }

    void write_stats(DataT* out) const {
        for (uint32_t j = 0; j < _stats.size(); j++) {
            out[j] = _stats[j];
        }
    }

    void reset() {
        if (_stats.size() == 0) {
            auto sums_size = _is_dense ? _ncols : _which_dims.size();
            for (size_t i = 0; i < sums_size; i++) {
                _stats.push_back(0);
            }
        } else {
            for (size_t i = 0; i < _stats.size(); i++) {
                _stats[i] = 0;
            }
        }
    }

    uint32_t nrows() const { return _nrows(); }
    uint16_t ncols() const { return _ncols(); }

private:

    std::vector<uint16_t> _which_dims;
    std::vector<dist_t> _stats;
    uint32_t _nrows;
    uint16_t _ncols;
    bool _is_dense;
};

template<class DataT, int OpE>
QueryResult sliding_binary_op(const QueryParams& q,
    const DataInfo& di, const DataT* buff)
{
    auto window_nrows = q.window_nrows > 0 ? q.window_nrows : di.nrows;
    // printf("actually running sliding max! window nrows, ncols, stride "
    //     " = %lld, %lld, %lld\n", window_nrows, q.window_ncols, q.window_stride);

    // figure out how long data is, and how many window positions we have
    auto nrows = di.nrows;
    int64_t last_window_start_row = nrows - window_nrows;
    int64_t nwindows = nrows - window_nrows + 1;

    // auto ret_size = nrows * di.ncols;
    size_t sparse_ncols = q.which_cols.size();
    bool sparse = sparse_ncols > 0;
    auto ret_ncols = sparse ? sparse_ncols : di.ncols;
    auto ret_size = nwindows * ret_ncols;

    QueryResult ret;
    auto& ret_vals = QueryResultValsRef<DataT>{}(ret);
    ret_vals.resize(ret_size);

    if (nwindows < 1) { return ret; }

    if (di.storage_order == ROWMAJOR) {
        OnlineBinaryOpRowmajor<DataT, OpE> stat(window_nrows, di.ncols, q.which_cols);
        stat.init(buff);
        auto ret_ptr = ret_vals.data();
        stat.write_stats(ret_ptr);
        for (size_t row = window_nrows; row < last_window_start_row; row++) {
            auto old_row = row - window_nrows;
            auto old_ptr = buff + di.ncols * old_row;
            auto new_ptr = buff + di.ncols * row;
            auto ret_row_ptr = ret_ptr + (old_row + 1) * ret_ncols;
            stat.update(old_ptr, new_ptr);
            stat.write_stats(ret_row_ptr);
        }
        return ret;
    }

    // column-major; treat each col as 1D rowmajor, and also write out results
    // in column-major order
    auto which_cols = q.which_cols;
    if (!sparse) {
        for (int i = 0; i < di.ncols; i++) {
            which_cols.push_back(i);
        }
    }
    for (int j_idx = 0; j_idx < which_cols.size(); j_idx++) {
        OnlineBinaryOpRowmajor<DataT, OpE> stat(window_nrows, 1);
        auto buff_ptr = buff + di.nrows;
        auto ret_ptr = ret_vals.data() + di.nrows; // write to ret in colmajor order
        stat.init(buff_ptr);
        stat.write_stats(ret_ptr);
        for (size_t row = window_nrows; row < last_window_start_row; row++) {
            auto old_row = row - window_nrows;
            auto ret_row_ptr = ret_ptr + (old_row + 1);
            stat.update(buff_ptr + old_row, buff_ptr + row);
            stat.write_stats(ret_row_ptr);
        }
    }
    return ret;
}

// XXX these don't actually do a sliding window---just a (slow) reduction; also,
// they cause segfaults sometimes...

template<class DataT>
QueryResult sliding_min(const QueryParams& q,
    const DataInfo& di, const DataT* buff)
{
    // printf("running sliding min query!\n");
    return sliding_binary_op<DataT, OpE::MIN>(q, di, buff);
}
template<class DataT>
QueryResult sliding_max(const QueryParams& q,
    const DataInfo& di, const DataT* buff)
{
    // printf("running sliding max query!\n");
    return sliding_binary_op<DataT, OpE::MAX>(q, di, buff);
}



// TODO eigen to impl these


template<typename T, typename T2>
static inline auto div_round_up(T x, T2 y) -> decltype(x + y) {
    return (x / y) + ((x % y) > 0);
}


template <typename DataT, int elem_sz=sizeof(DataT)>
void reduce_sum_avx2_rowmajor_ax0(const DataT* buff, int64_t nrows,
    int64_t ncols, int32_t* out)
{
    // using DataT = uint8_t;
    // static const int elem_sz = 1;
    static const int vector_sz_bytes = 32;
    static const int vector_sz_elems = vector_sz_bytes / elem_sz;
    static const int out_stripes_per_in_stripe = 4 / elem_sz;

    int nstripes_in = div_round_up(ncols, vector_sz_elems);
    int nstripes_out = nstripes_in * out_stripes_per_in_stripe;
    // int final_stripe_offset = ncols - vector_sz_elems;
    // int last_stripe_idx = nstripes - 1;

    __m256i sums[nstripes_out];
    for (size_t v = 0; v < nstripes_out; v++) {
        sums[v] = _mm256_setzero_si256();
    }

    // SIMD sum until the last few rows; we leave enough rows at the end
    // so that this loop can read past the row end
    int padding_nrows = div_round_up(ncols, vector_sz_elems);
    size_t i = 0;
    const DataT* read_ptr = buff;
    for ( ; i < nrows - padding_nrows; i++) {
        auto row_start_ptr = buff + i * ncols;
        for (size_t v = 0; v < nstripes_in; v++) {
            auto load_addr = row_start_ptr + v * vector_sz_elems;
            __m256i x = _mm256_load_si256((const __m256i*)load_addr);
            auto x_low = _mm256_extracti128_si256(x, 0);
            auto x_high = _mm256_extracti128_si256(x, 1);
            auto out_idx = v * out_stripes_per_in_stripe;
            if (elem_sz == 1) {
                __m256i x0, x1, x2, x3;
                if (std::is_signed<DataT>::value) {     // int8
                    x0 = _mm256_cvtepi8_epi32(x_low);
                    x1 = _mm256_cvtepi8_epi32(_mm_slli_si128(x_low, 8));
                    x2 = _mm256_cvtepi8_epi32(x_high);
                    x3 = _mm256_cvtepi8_epi32(_mm_slli_si128(x_high, 8));
                } else {                                // uint8
                    x0 = _mm256_cvtepu8_epi32(x_low);
                    x1 = _mm256_cvtepu8_epi32(_mm_slli_si128(x_low, 8));
                    x2 = _mm256_cvtepu8_epi32(x_high);
                    x3 = _mm256_cvtepu8_epi32(_mm_slli_si128(x_high, 8));
                }

                sums[out_idx + 0] = _mm256_add_epi32(sums[out_idx + 0], x0);
                sums[out_idx + 1] = _mm256_add_epi32(sums[out_idx + 1], x1);
                sums[out_idx + 2] = _mm256_add_epi32(sums[out_idx + 2], x2);
                sums[out_idx + 3] = _mm256_add_epi32(sums[out_idx + 3], x3);
            } else if (elem_sz == 2) {
                __m256i x0, x1;
                if (std::is_signed<DataT>::value) {     // int16
                    x0 = _mm256_cvtepi16_epi32(x_low);
                    x1 = _mm256_cvtepi16_epi32(x_high);
                } else {                                // uint16
                    x0 = _mm256_cvtepu16_epi32(x_low);
                    x1 = _mm256_cvtepu16_epi32(x_high);
                }
                sums[out_idx + 0] = _mm256_add_epi32(sums[out_idx + 0], x0);
                sums[out_idx + 1] = _mm256_add_epi32(sums[out_idx + 1], x1);
            } else {
                static_assert(elem_sz == 1 || elem_sz == 2,
                    "only element sizes of 1 and 2 bytes are supported!");
            }
        }
    }

    // dump simd vecs to normal arrays so we can write to them serially
    DataT sums_ar[nstripes_out * vector_sz_elems];
    for (size_t v = 0; v < nstripes_out; v++) {
        auto write_addr = &sums_ar[v * vector_sz_elems];
        _mm256_store_si256((__m256i*)write_addr, sums[v]);
    }

    // add in sums from tail rows we didn't operate on via SIMD stuff
    const DataT* tail_start_ptr = buff + ncols * i;
    for (; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            sums_ar[j] += buff[i * ncols] + j;
        }
    }
    // copy sums to output; didn't just write to output directly since
    // dumping the SIMD vects would write past the end
    for (size_t j = 0; j < ncols; j++) {
        out[j] = sums_ar[j];
    }
}


// template<class DataT>
// // static inline QueryResult poopinize(const QueryParams& q,
// static inline QueryResult frobnicate(const QueryParams& q,
//     const DataInfo& di, const DataT* buff)
// {
//     printf("called frobnicate!\n");
//     exit(1);
//     return QueryResult{};
// }

template<int ElemSz, class DataT>
// static inline QueryResult poopinize(const QueryParams& q,
static inline QueryResult reduce_contiguous(const QueryParams& q,
    const DataInfo& di, const DataT* buff)
{
    // // fprintf(stderr, "running reduce_contiguous\n");
    // printf("******* running reduce_contiguous using query type %d\n", q.type);
    // exit(1);

    // compute actual data type based on elemsz; buff is likely to be a byte*
    // regardless of what those bytes represent
    using RealDataType = typename ElemSizeTraits<ElemSz>::DataT;
    using RowmajorMat = Eigen::Map<const Eigen::Matrix<
        RealDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
    using ColmajorMat = Eigen::Map<const Eigen::Matrix<
        RealDataType, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> >;
    using RowVector = Eigen::Array<RealDataType, 1, Eigen::Dynamic, Eigen::RowMajor>;
    using RowVectorMap = Eigen::Map<RowVector>;
    using RowVector_I32 = Eigen::Array<int32_t, 1, Eigen::Dynamic, Eigen::RowMajor>;

    // using MutVector = Eigen::Map<Eigen::Vector<DataT, Eigen::Dynamic> >;
    // printf("running sliding min query!\n");
    // return sliding_binary_op<DataT, OpE::MIN>(q, di, buff);

    // printf("******* running reduce_contiguous\n");
    // std::cout << "RUN THIS CODE YOU RETARDED COMPILER" << std::endl;
    // exit(1);

    QueryResult ret;
    auto& ret_vals = QueryResultValsRef<RealDataType>{}(ret);
    auto& ret_vals_i32 = ret.vals_i32;

    bool use_i32_output = q.type != QUERY_MIN and q.type != QUERY_MAX;

    // ret_vals.resize(di.ncols);
    ret_vals.reserve(di.ncols);
    if (use_i32_output) { ret_vals_i32.reserve(di.ncols); }
    // ret.idxs.push_back(1);
    // ret.vals_u8.push_back(1);

    // return ret;

    const RealDataType* data_ptr = (const RealDataType*)buff;

    // printf("retvals elem_sz: %lu\n", sizeof(ret_vals[0]));
    // printf("retvals size: %lu\n", ret_vals.size());
    // printf("retvals capacity: %lu\n", ret_vals.capacity());

    // MutRowmajorMat ret_vec(ret_vals.data(), 1, di.ncols);
    RowVectorMap ret_vec(ret_vals.data(), 1, di.ncols);
    auto ret_buff_i32 = ret_vals_i32.data();

    // printf("nrows, ncols: %d, %d\n", (int)di.nrows, (int)di.ncols);
    // exit(1);

    // TODO uncomment below
    //
    if (di.storage_order == ROWMAJOR) {
        RowmajorMat mat(data_ptr, di.nrows, di.ncols);
        RowVector tmp(di.ncols);

        // printf("*** using rowmajor storage order\n");

        // printf("rowmajor size of reduce_sum: %d\n", (int)mat.colwise().sum().size());

        // NOTE: using same eigen cols as below is insanely slow, so we
        // have to implement these reductions ourselves
        switch (q.type) {
        case QUERY_MEAN:
            // printf("running rowmajor query_mean!\n");
            reduce_sum_avx2_rowmajor_ax0(mat.data(), di.nrows, di.ncols, ret_buff_i32);
            // for (size_t j = 0; j < di.ncols; j++) {
            //     ret_buff_i32[j] /= di.nrows;
            // }
        case QUERY_SUM:
            reduce_sum_avx2_rowmajor_ax0(mat.data(), di.nrows, di.ncols, ret_buff_i32);
        case QUERY_MIN:
            ret_vec = mat.row(0);
            for (size_t i = 1; i < di.nrows; i++) {
                ret_vec = ret_vec.min(mat.row(i).array());
            }
            break;
        case QUERY_MAX:
            ret_vec = mat.row(0);
            for (size_t i = 1; i < di.nrows; i++) {
                ret_vec = ret_vec.max(mat.row(i).array());
            }
            break;
        case QUERY_NORM:
            ret_vec = mat.colwise().squaredNorm(); break;
        default:
            printf("Unsupported query type for contiguous data: %d!\n",
                (int)q.type); exit(1);
        }
    } else { // XXX pretty sure all of these (except min/max) ignore overflows
        ColmajorMat mat(data_ptr, di.nrows, di.ncols);

        // printf("colmajor size of reduce_sum: %d\n", (int)mat.colwise().sum().size());
        // printf("colmajor number of rows, cols: %d, %d\n", (int)mat.rows(), (int)mat.cols());
        // printf("colmajor size of ret_vec: %d\n", (int)ret_vec.size());
        // printf("colmajor last val in sum: %d\n", (int)mat.colwise().sum()(0, 17));

        switch (q.type) {
        case QUERY_MEAN:
            ret_vec = mat.colwise().mean(); break;
        case QUERY_SUM:
            /// XXX this will overflow for small int types
            ret_vec = mat.colwise().sum(); break;
        case QUERY_MIN:
            ret_vec = mat.colwise().minCoeff(); break;
        case QUERY_MAX:
            ret_vec = mat.colwise().maxCoeff(); break;
        case QUERY_NORM:
            ret_vec = mat.colwise().squaredNorm(); break;
        default:
            printf("Unsupported query type for contiguous data: %d!\n",
                (int)q.type); exit(1);
        }
    }
    return ret;
}

#endif // QUERY_MIN_MAX_HPP
