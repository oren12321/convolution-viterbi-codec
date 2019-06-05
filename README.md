# Convolution Viterbi Codec

This software is a highly flexible parallelized codec that composed from a convolution encoder and a Viterbi decoder.

The encoder uses convolution that “inflates” a given data, such that every bit is expanded by a shift register with memory K and n degrees polynomial generator. The encoding process contains the following steps:
* Lookup tables creation (for fast running of the shift register and the decoder) that contains, according to the current state of the shift register, the next state, the previous state and the possible encoded symbols (by convolution), and all that by pushing 0 or 1 to its memory.
* Scanning of the data bits. Each bit is pushed to the shift register and a code for this bit will be created according to the mentioned lookup tables.

The decoder uses a custom Viterbi algorithm that adjusted for this purpose. The decoding process contains the following steps:
* Trellis diagram creation according to the states and the encoded data. Since the data is binary, the transitions in the table are increase every stage by a factor of 2.
* Backtrace according to the table minimal weights. The minimal weight in the end of the table represents the number of error bits between the encoder to the decoder.

Comment: A bit can be encoded for symbol bigger than 64 bits. In this case, each symbol will be a vector of 64 bit slices.

The program written in C++.and it can deal with encoding of many primitive types by templating. The program is highly parallelized and can be compiled with every combination (or none for serial running) of the following flags: OpenMP, MPI, CUDA.

The parallelism is as follows:
* MPI – Each process receives an equal (as possible) slice of the encoded data. It can be done since big slices of the encoded data are independent (there is a negligible cost in the errors correction but K and n can be increased for compensation).
* OpenMP – The Trellis diagram rows will be divided equally between the CPUs.
* CUDA – There are two possibilities:
  * Lookup tables creation for short time initialization of the codec.
  * Replacement of OpenMP (When CUDA selected, it will replace OpenMP).