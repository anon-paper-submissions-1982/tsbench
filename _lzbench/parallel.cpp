//
// parallel.cpp
// Created by D Blalock 2018-2-28
//

#include "parallel.h"

#include <chrono>  // TODO rm
#include <iostream>  // TODO rm
#include <future>
#include <thread>

#include "output.h"
#include "preprocessing.h"
#include "query.hpp"
#include "util.h"


size_t _decomp_and_query(lzbench_params_t *params, const compressor_desc_t* desc,
    const uint8_t* comprbuff, size_t comprsize, uint8_t* outbuf, size_t outsize,
    bool already_materialized,
    size_t param1, size_t param2, char* workmem)
{
    // printf("decomp_and_query: running '%s' with insize %lu, outsize %u!\n", desc->name, comprsize, outsize);
    // bool already_materialized = strings_equal(desc->name, "materialized");

    compress_func decompress = desc->decompress;

    size_t dlen = -1;
    if (!already_materialized) {
        // if (comprsize == outsize || true) { // TODO rm
        if (comprsize == outsize) { // uncompressed
            memcpy(outbuf, comprbuff, comprsize);
            dlen = comprsize;
        } else {
            // printf("about to decomp using compressor: %s\n", desc->name);
            dlen = decompress((char*)comprbuff, comprsize, (char*)outbuf,
                              outsize, param1, param2, workmem);
        }
        undo_preprocessors(params->preprocessors, outbuf, dlen,
            params->data_info.element_sz);

        // prevent compiler from not running above command (hopefully...)
        if (params->verbose >= 999 || dlen >= ((int64_t)1) << 31) {
            size_t cmn = common(comprbuff, outbuf, outsize);
            LZBENCH_PRINT(999, "ERROR in %s: only first %d / %d decompressed bytes were correct\n",
                desc->name, (int32_t)cmn, (int32_t)outsize);
        }
    } else {
        dlen = outsize;
    }

    // run query if one is specified
    auto qparams = params->query_params;
    if (qparams.type != QUERY_NONE) {
        // printf("got query type: %d; about to run a query...\n", qparams.type);
        // auto& dinfo = params->data_info;
        DataInfo dinfo = params->data_info;
        if (dinfo.ncols < 1) {
            printf("ERROR: Must specify number of columns in data to run query!\n");
            exit(1);
        }
        dinfo.nrows = dlen / (dinfo.ncols * dinfo.element_sz);
        // printf("dlen: %lld\n", (int64_t)dlen);
        // printf("dinfo nrows, ncols, size: %lu, %lu, %lu\n",
        //     dinfo.nrows, dinfo.ncols, dinfo.nrows * dinfo.ncols);
        QueryResult result = run_query(
            params->query_params, dinfo, outbuf);
        // QueryResult result = frobnicate(                         // TODO rm
        //     params->query_params, dinfo, outbuf);
        // printf("ran query type: %d\n", qparams.type);
        // printf("number of idxs in result: %lu\n", result.idxs.size());


        // hack so it can't pull the above check out of the loop; dummy
        // can be any u8 but next line will always add 0, although compiler
        // doesn't know this (it's 0 because element_sz is in {1,2})
        auto dummy = result.vals_u8.size() > 0 ? result.vals_u8[0] : 0;
        params->verbose += result.idxs.size() > ((int64_t)1e9) ? dummy : 0;

        // prevent compiler from optiming away query
        // XXX does it actually have this effect? could pull this check
        // out of the loop and do nothing if condition is false
        if (params->verbose > 999) {
            printf("query u8 result: ");
            for (auto val : result.vals_u8) { printf("%d ", (int)val); }
            printf("\n");
            printf("query u16 result: ");
            for (auto val : result.vals_u16) { printf("%d ", (int)val); }
            printf("\n");
        }
    }

    return dlen;
}



void parallel_decomp(lzbench_params_t *params,
    std::vector<size_t>& chunk_sizes, const compressor_desc_t* desc,
    std::vector<size_t> &compr_sizes, const uint8_t *inbuf, uint8_t *outbuf,
    uint8_t* tmpbuf, bench_rate_t rate, std::vector<uint64_t> comp_times,
    size_t param1, size_t param2, char* workmem)
{
    // printf("calling parallel decomp for algorithm (T=%d): %s!\n", params->nthreads, desc->name);
    // printf("calling parallel decomp for algorithm %s!\n", desc->name);
    // if (params) {
    //     printf("using nthreads: %d\n", params->nthreads);
    // }

    std::vector<uint64_t> compressed_chunk_starts;
    compressed_chunk_starts.push_back(0);
    for (auto sz : compr_sizes) {
        compressed_chunk_starts.push_back(compressed_chunk_starts.back() + sz);
    }
    compressed_chunk_starts.pop_back(); // last one is just an end idx

    // printf("compr start idxs: ");
    // for (auto start_idx : compressed_chunk_starts) {
    //     printf("%lld, ", start_idx);
    // }
    // printf("\n");

    // printf("param1, param2 = %lu, %lu\n", param1, param2);
    // if (param1 != 80) {
    //     printf("param1 is %lu, not 80!\n", param1);
    // }

    uint64_t run_for_nanosecs = (uint64_t)params->dmintime*1000*1000;

    using result_t = std::tuple<int64_t, int64_t>;

    int nthreads = params->nthreads;
    // std::vector<int64_t> total_scanned_sizes(nthreads);
    // std::vector<std::future<int64_t>> total_scanned_sizes(nthreads);
    std::vector<std::future<result_t>> thread_results_futures;
    std::vector<result_t> thread_results(nthreads);
    // std::vector<result_t> total_scanned_sizes(nthreads);
    // std::vector<std::thread> threads(nthreads);

    auto max_chunk_sz = chunk_sizes[0];
    auto total_raw_sz = 0;
    for (auto sz : chunk_sizes) {
        if (sz > max_chunk_sz) { max_chunk_sz = sz; }
        total_raw_sz += sz;
    }

    bool already_materialized = strings_equal(desc->name, "materialized");

    // for (int i = 0; i < nthreads; i++) {
        // auto& this_total = total_scanned_sizes[i];
        // size_t* this_total = total_scanned_sizes[i];
        // threads[i] = std::thread([&total_scanned_sizes[i]] {
        auto run_in_thread =
            // [i, run_for_nanosecs, max_chunk_sz, total_raw_sz, inbuf, nthreads,
            [run_for_nanosecs, max_chunk_sz, total_raw_sz, inbuf, nthreads,
                params, desc, compr_sizes, chunk_sizes, rate,
                compressed_chunk_starts,
                already_materialized,
                // &total_scanned_sizes,
                param1, param2, workmem](int i) {

            // std::this_thread::sleep_for(std::chrono::duration<double, std::milli>(100));
            // std::this_thread::sleep_for(std::chrono::duration<double>(1));
            // std::this_thread::sleep_for(std::chrono::duration<double>(2));
            // return (int64_t)total_raw_sz;

            //
            // TODO uncomment below here
            // EDIT: this scales linearly, so issue is something below here...
            //

            int64_t decomp_sz = 0;

            bench_timer_t t_end;
            int64_t max_iters = run_for_nanosecs > 0 ? 1000*1000*1000 : -1;
            int64_t niters = 0;
            auto num_chunks = compressed_chunk_starts.size();
            // XXX this is an ugly way to check this

            // printf("using num_chunks: %lld\n", (int64_t)num_chunks);
            // printf("using chunk sizes:"); for (auto sz : chunk_sizes) { printf("%lld, ", (int64_t)sz); } printf("\n");

            // printf("max chunk sz: %lu\n", max_chunk_sz);
            uint8_t* decomp_buff = alloc_data_buffer(max_chunk_sz + 4096);

            int64_t elapsed_nanos = 0;

            bench_timer_t t_start;
            GetTime(t_start);

            do {
                // run multiple iters betwen rtsc calls to avoid sync overhead
                // use nthreads iters as a heuristic so syncs/sec is constant
                // for (int it = 0; it < num_chunks*nthreads; it++) {
                for (int it = 0; it < nthreads; it++) {
                // for (int it = 0; it < 1; it++) { // TODO uncomment above
                    auto chunk_idx = rand() % num_chunks;
                    // auto chunk_idx = it % num_chunks;
                    auto inptr = inbuf + compressed_chunk_starts[chunk_idx];
                    // auto inptr = inbuf; // TODO uncomment above after debug
                    auto insize = compr_sizes[chunk_idx];
                    auto rawsize = chunk_sizes[chunk_idx];

                    _decomp_and_query(params, desc, inptr, insize,
                        decomp_buff, rawsize, already_materialized,
                        param1, param2, workmem);

                    // std::this_thread::sleep_for(std::chrono::duration<double, std::milli>(100));
                    // std::this_thread::sleep_for(std::chrono::duration<double>(1));

                    // this_total += rawsize;
                    decomp_sz += rawsize;
                    // total_scanned_sizes[i] += rawsize;
                    // *(this_total) = *(this_total) + rawsize;
                    niters++;
                }

                // check whether we're done
                GetTime(t_end);
                elapsed_nanos = GetDiffTime(rate, t_start, t_end);

                bool done = elapsed_nanos >= run_for_nanosecs;
                done = done || ((max_iters >= 0) && (niters >= max_iters));
                if (done) {
                // if (true) {
                    LZBENCH_PRINT(8, "%d) elapsed iters, time: %lld, %lld/%lldns\n",
                        i, niters, elapsed_nanos, run_for_nanosecs);
                }
                if (done) { break; }
            } while (true);

            // auto total_comp_size = 0;
            // for (auto sz: compr_sizes) {
            //     total_comp_size += sz;
            // }
            // size_t cmn = common(inbuf, decomp_buff, total_raw_sz);
            // if (cmn < insize) {
            // printf("about to check whether decomp is correct...\n");
            // size_t cmn = common(inbuf, decomp_buff, max_chunk_sz);
            // if (cmn < max_chunk_sz) {
            //     LZBENCH_PRINT(999, "ERROR in %s: only first %d / %d decompressed bytes were correct\n",
            //     desc->name, (int32_t)cmn, (int32_t)max_chunk_sz);
            // }

            free_data_buffer(decomp_buff);

            // return decomp_sz;
            return std::make_tuple(decomp_sz, elapsed_nanos);
        // });
        };
    // }

    // TODO uncomment below
    //
    // auto debug_lambda = [](int64_t i) { std::cout << i << "\n"; return i; };
    for (int i = 0; i < nthreads; i++) {
        thread_results_futures.push_back(std::async(run_in_thread, i));
        // thread_results_futures.push_back(std::async(debug_lambda, i));
    }
    // printf("about to try get()ing all the futures...\n");
    for (int i = 0; i < nthreads; i++) {
        thread_results[i] = thread_results_futures[i].get();
        // thread_results_futures[i].get();
    }


    // // Single threaded version (for debugging)
    // for (int i = 0; i < nthreads; i++) {
    //     thread_results[i] = run_in_thread(i);
    // }

    // for (auto& t : threads) {
    //     t.join();
    // }

    // printf("total sizes: ");
    // for (auto res : thread_results) {
    //     printf("(%lldB, %lldns), ", std::get<0>(res), std::get<1>(res));
    // }
    // printf("\n");

    // compute total amount of data all the threads got through
    int64_t total_scanned_bytes = 0;
    int64_t total_cpu_time = 0;
    for (auto res : thread_results) {
        total_scanned_bytes += std::get<0>(res);
        total_cpu_time += std::get<1>(res);
    }
    double thruput_bytes_per_cpu_ns = total_scanned_bytes / (double)total_cpu_time;
    int64_t thruput_MB_per_sec = (int64_t)((thruput_bytes_per_cpu_ns * 1000) * nthreads);

    // if (!run_for_nanosecs) { // this case shouldn't be used for real results
    //     bench_timer_t t_end;
    //     GetTime(t_end);
    //     // printf("WARNING: minimum run time not specified\n");
    //     run_for_nanosecs = GetDiffTime(rate, t_start, t_end);
    // }
    // auto run_for_usecs = run_for_nanosecs / 1000;
    // auto thruput_MB_per_sec = total_scanned_bytes / run_for_usecs;
    // printf(">> \1%s avg thruput: %lld(MB/s)\n", desc->name, thruput_MB_per_sec);
    // printf(">> \1%s avg thruput: %lld(MB/s)\n", desc->name, thruput_MB_per_sec);

    size_t complen = 0;
    for (auto sz : compr_sizes) { complen += sz; }

    bool decomp_error = false;
    std::vector<uint64_t> decomp_times {total_cpu_time};
    size_t insize = total_scanned_bytes;
    print_stats(params, desc, param1, comp_times, decomp_times, insize,
        complen, decomp_error);

    // printf("------------------------");
}



