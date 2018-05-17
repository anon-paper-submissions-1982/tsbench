
#ifndef QUERY_HPP
#define QUERY_HPP

#include "query_common.h"
#include "query_mean.hpp"
#include "query_minmax.hpp"
#include "query_reduce_row.hpp"

// if we weren't just benchmarking, would need to return something in all
// of these

// template<class DataT>
// void sliding_mean(const QueryParams& q, const DataInfo& di, const DataT* buff) {

// }

// template<class DataT>
// QueryResult sliding_min(const QueryParams& q, const DataInfo& di,
//     const DataT* buff)
// {
//     return QueryResult{}; // TODO
// }

// template<class DataT>
// QueryResult sliding_max(const QueryParams& q, const DataInfo& di,
//     const DataT* buff)
// {
//     return QueryResult{}; // TODO
// }

// template<class DataT>
// QueryResult sliding_l2(const QueryParams& q, const DataInfo& di,
//     const DataT* buff)
// {
//     if (q.reduction == REDUCE_NONE) {

//     } else if (q.reduction == REDUCE_THRESH) {

//     } else if (q.reduction == REDUCE_TOP_K) {

//     } else {
//         printf("Invalid reduction %d for L2 query!\n", (int)q.reduction);
//         exit(1);
//     }
//     return QueryResult{}; // TODO
// }

// template<class DataT>
// QueryResult sliding_dot(const QueryParams& q, const DataInfo& di,
//     const DataT* buff)
// {
//     if (q.reduction == REDUCE_NONE) {

//     } else if (q.reduction == REDUCE_THRESH) {

//     } else if (q.reduction == REDUCE_TOP_K) {

//     } else {
//         printf("Invalid reduction %d for dot product query!\n",
//             (int)q.reduction);
//         exit(1);
//     }
//     return QueryResult{}; // TODO
// }

template<class DataT>
QueryResult corr(const QueryParams& q, const DataInfo& di,
    const DataT* buff)
{

}


template<class DataT>
QueryResult _frobnicate(const QueryParams& q, const DataInfo& di, const DataT* buff) {
    printf("actually running frobnicate; query_type=%d!\n", (int)q.type);
    return QueryResult{};
}

template<class DataT>
QueryResult run_query(const QueryParams& q, const DataInfo& di, const DataT* buff) {

    // printf("actually running run_query; query_type=%d!\n", (int)q.type);

    // QueryResult ret;
    switch (di.element_sz) {
    case 1: return reduce_contiguous<1>(q, di, buff);
    case 2: return reduce_contiguous<2>(q, di, buff);
    // case 1: return frobnicate(q, di, buff);
    // case 2: return frobnicate(q, di, buff);
    default:
        printf("Invalid element size %d!\n", (int)di.element_sz); exit(1);
        exit(1);
    }
    // switch (q.type) {
    //     // case QUERY_MEAN: ret = sliding_mean(q, di, buff); break;
    //     case QUERY_MEAN: ret = reduce_contiguous(q, di, buff); break;
    //     // XXX "sliding" min and max actually just write out min/max seen
    //     // so far, which is a weird thing to do
    //     case QUERY_MIN: ret = reduce_contiguous(q, di, buff); break;
    //     case QUERY_MAX: ret = reduce_contiguous(q, di, buff); break;
    //     // case QUERY_L2: ret = sliding_l2(q, di, buff); break;
    //     // case QUERY_DOT: ret = sliding_dot(q, di, buff); break;
    //     case QUERY_SUM: ret = reduce_contiguous(q, di, buff); break;
    //     default:
    //         printf("Invalid query type %d!\n", (int)q.type); exit(1);
    // }

    // TODO check if query has a reduction here and do it in this one place
    // if so

    // if (q.reduction == REDUCE_NONE) {

    // } else if (q.reduction == REDUCE_THRESH) {

    // } else if (q.reduction == REDUCE_TOP_K) {

    // } else {
    //     printf("Unsupported reduction %d!\n", (int)q.reduction);
    //     exit(1);
    // }
    return QueryResult{}; // can't happen, in theory
}

#endif // QUERY_HPP
