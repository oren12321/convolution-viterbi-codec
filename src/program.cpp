/*
 * program.cpp
 *
 *  Created on: Oct 31, 2014
 *      Author: Oren
 */

#include <iostream>
#include <sys/time.h>
#include <stddef.h>

#include <vector>
#include <list>

#include <unistd.h>

#include "api/CViterbiEncoderDecoder.h"

#include <cmath>

using namespace std;

#define MEASURE_CODE(code, start_time_var, end_time_var, description_str) \
	gettimeofday(&start_time_var, NULL); \
	code \
	gettimeofday(&end_time_var, NULL); \
	cout << description_str << " [us] = " << (end_time_var.tv_sec - start_time_var.tv_sec) * 1000000 + (end_time_var.tv_usec - start_time_var.tv_usec) << endl;

template<typename T>
void print_array(vector<T> const& data);

unsigned long estimated_heap_size_in_bytes(unsigned long k, unsigned long n, unsigned long data_size, unsigned long type_size_in_bytes);

template<typename T>
void encoding_decoding_test(unsigned int const k, unsigned int const n, vector<T> const& data, bool const debug);

void tests(void);

int main(int argc, char* argv[])
{
	tests();

	return 0;
}

void tests(void)
{
#define SIZE 100
#define N 20
#define K 12
#define TYPE unsigned long

	cout << "estimated used memory [mb] = " << ((double)estimated_heap_size_in_bytes(K,N,SIZE,sizeof(TYPE)) / 1024.0) / 1024.0 << endl;

	vector<TYPE> data(SIZE);
	for(int i = 0; i < SIZE; i++)
	{
		data[i] = i + 1;
	}

	encoding_decoding_test<TYPE>(K, N, data, false);
}

unsigned long estimated_heap_size_in_bytes(unsigned long k, unsigned long n, unsigned long data_size, unsigned long type_size_in_bytes)
{
	unsigned long result = 0;

	unsigned long k_over_2 = pow(2, k);
	unsigned long number_of_symbols = (data_size + k) * type_size_in_bytes * 8;

	// polynomials
	result += n * sizeof(unsigned int);
	// states lookup tables
	result += k_over_2 * 2 * sizeof(unsigned int) * 2;
	// symbols lookup table
	result += k_over_2 * sizeof(std::list<unsigned int>) * 2 + k_over_2 * 2 * ceil((double)n / (sizeof(unsigned int) * 8)) * sizeof(unsigned int);
	// convolution
	result += number_of_symbols * sizeof(std::list<unsigned int>*);
	// trellis
	result += k_over_2 * sizeof(unsigned int) + k_over_2 * number_of_symbols * sizeof(unsigned int);
	// viterbi
	result += number_of_symbols * sizeof(unsigned int) + data_size * type_size_in_bytes;

	return result;
}

template<typename T>
void encoding_decoding_test(unsigned int const k, unsigned int const n, vector<T> const& data, bool const debug)
{
	timeval start;
	timeval end;

	cout << "provided data :" << endl << "\t";
	print_array<T>(data);

	MEASURE_CODE(
			CViterbiEncoderDecoder<T> vEnDe(k, n, data.size());,
			start,
			end,
			"initialization"
	)

	MEASURE_CODE(
			vector<list<unsigned int>* > const& encoded = vEnDe.Encode(data);,
			start,
			end,
			"encoding"
	)

	MEASURE_CODE(
			vector<T> const& decoded = vEnDe.Decode(encoded);,
			start,
			end,
			"decoding"
	)

	cout << "decoded data :" << endl << "\t";
	print_array<T>(decoded);

	if(debug)
	{
		cout << "debug data :" << endl;
		cout << vEnDe << endl;
	}
}

template<typename T>
void print_array(vector<T> const& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		cout << data[i] << " ";
	}
	cout << endl;
}
