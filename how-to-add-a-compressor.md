
Steps

 1. Add all the relevant code in some subdirectory.
 1. Go to `_lzbench/lzbench.h`
    1. Increment `LZBENCH_COMPRESSOR_COUNT` (around line 140)
    1. Add the relevant info in the array `comp_desc`. The fields are described in the struct immediately above it (`compressor_desc_t`).
        1. name and version are self explanatory
        1. First and last level refer to levels of compression. By convention, higher levels mean more compression.
        1. TODO I don't know what `additional_param` means.
        1. `max_block_size` is self-explanatory if your compressor splits the data into blocks and compresses each block (which it probably does).
        1. `compress` and `decompress` need to be pointers to functions matching the signatures immediately above the struct definition. Note that these should be wrapper functions (which we'll define in a moment), not your actual compression and decompression functions.
        1. The same is true for `init` and `deinit`, though these can be `NULL`.
    1. Optionally, add an alias to the `alias_desc` array below and increment `LZBENCH_ALIASES_COUNT`.
 1. Go to `_lzbench/compressors.h` and add in an `ifdef` block declaring the functions you specified in the `comp_desc` array.
 1. Go to `_lzbench/compressors.cpp` and add in an `ifdef` block defining these functions, presumably by calling your own compressor's top-level compress and decompress functions. Note that you should only `#include` your header within this `ifdef`.

 1. Go to the Makefile. Add in a variable containing all the .o files your compressor needs. Then add this variable to the list under the `lzbench` target.
    1. If you need special flags, add a custom rule above this build rule, copying the (simple) structure used by other such rules.


As an example, consider the contents of the `simple_example` directory, and the associated `example_compress`, `example_decompress`, and other functions.

