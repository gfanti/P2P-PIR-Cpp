//  Percy++ Copyright 2007,2012,2013 Ian Goldberg <iang@cs.uwaterloo.ca>,
//  Casey Devet <cjdevet@cs.uwaterloo.ca>,
//  Paul Hendry <pshdenry@uwaterloo.ca>,
//  Ryan Henry <rhenry@cs.uwaterloo.ca>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of version 2 of the GNU General Public License as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  There is a copy of the GNU General Public License in the COPYING file
//  packaged with this plugin; if you cannot find it, write to the Free
//  Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//  02110-1301 USA

#ifndef __PERCYSERVER_IMPL_H__
#define __PERCYSERVER_IMPL_H__

#include <iostream>
#include <array>
#include <algorithm>
#include <typeinfo>
#include "percyparams.h"
#include "gf2e.h"
#include "pulse.h"

// this is also defined in percyserver.h!
// constexpr dbsize_t WORDS_PER_BLOCK = 1024;

// Compute single-query outputs
template <typename GF2E_Element>
inline void compute_outputvec_single(
        const GF2E_Element *data, GF2E_Element *inputvec, GF2E_Element *outputvec,
        dbsize_t num_blocks, dbsize_t words_per_block, const vector<dbsize_t> &unsynchronized_files);

template <>
inline void compute_outputvec_single<GF28_Element>(
        const GF28_Element *data, GF28_Element *inputvec, GF28_Element *outputvec,
        dbsize_t num_blocks, dbsize_t words_per_block, const vector<dbsize_t> &unsynchronized_files) {
    const GF28_Element *block = data;
    for (unsigned int j = 0; j < num_blocks; ++j) {
        const GF28_Element *multrow = GF28_mult_table[inputvec[j]];
        const GF28_Element *blockc = block;
        GF28_Element *oc = outputvec;
        GF28_Element *oc_end = oc + (words_per_block & ~7);
        while (oc < oc_end) {
            uint64_t accum = (uint64_t) multrow[*(blockc++)];
            accum |= (uint64_t) multrow[*(blockc++)] << 8;
            accum |= (uint64_t) multrow[*(blockc++)] << 16;
            accum |= (uint64_t) multrow[*(blockc++)] << 24;
            accum |= (uint64_t) multrow[*(blockc++)] << 32;
            accum |= (uint64_t) multrow[*(blockc++)] << 40;
            accum |= (uint64_t) multrow[*(blockc++)] << 48;
            accum |= (uint64_t) multrow[*(blockc++)] << 56;
            *((uint64_t *) oc) ^= accum;
            oc+=8;
        }
        for (unsigned int c = 0; c < (words_per_block & 7); ++c) {
            *(oc++) ^= multrow[*(blockc++)];
        }
        block += words_per_block;
    }
}

template <>
inline void compute_outputvec_single<GF216_Element>(
        const GF216_Element *data, GF216_Element *inputvec, GF216_Element *outputvec,
        dbsize_t num_blocks, dbsize_t words_per_block, const vector<dbsize_t> &unsynchronized_files) {
    const GF216_Element *block = data;
    for (dbsize_t j = 0; j < num_blocks; ++j) {
        if (std::find(unsynchronized_files.begin(), unsynchronized_files.end(), j)==unsynchronized_files.end()){
            GF216_Element inpv_j = inputvec[j];
            if (inpv_j != 0) {
                const GF216_Element *blockc = block;
                GF216_Element log_j = GF216_log_table[inpv_j];
                const GF216_Element *start = GF216_exp_table + log_j;
                GF216_Element *oc = outputvec;
                GF216_Element *oc_end = oc + (words_per_block & ~3);
                GF216_Element block_c;
                while(oc < oc_end) {
                    uint64_t accum = 0;
                    block_c = *(blockc++);
                    if (block_c != 0) {
                        GF216_Element log_c = GF216_log_table[block_c];
                        accum |= (uint64_t) start[log_c];
                    }
                    block_c = *(blockc++);
                    if (block_c != 0) {
                        GF216_Element log_c = GF216_log_table[block_c];
                        accum |= (uint64_t) start[log_c] << 16;
                    }
                    block_c = *(blockc++);
                    if (block_c != 0) {
                        GF216_Element log_c = GF216_log_table[block_c];
                        accum |= (uint64_t) start[log_c] << 32;
                    }
                    block_c = *(blockc++);
                    if (block_c != 0) {
                        GF216_Element log_c = GF216_log_table[block_c];
                        accum |= (uint64_t) start[log_c] << 48;
                    }
                    *((uint64_t *) oc) ^= accum;
                    oc+=4;
                }
                for (dbsize_t c = 0; c < (words_per_block & 3); ++c, ++oc) {
                    block_c = *(blockc++);
                    if (block_c != 0) {
                        GF216_Element log_c = GF216_log_table[block_c];
                        *oc ^= start[log_c];
                    }
                }
            }
        }
        block += words_per_block;
    }
}

// Compute multi-query outputs
template <typename GF2E_Element>
inline void compute_outputvec_multi(
        const GF2E_Element *data, GF2E_Element *inputvec, GF2E_Element *outputvec,
        nqueries_t num_queries, dbsize_t num_blocks, dbsize_t words_per_block,
        const vector<dbsize_t> &unsynchronized_files);

template <>
inline void compute_outputvec_multi<GF28_Element>(
        const GF28_Element *data, GF28_Element *inputvec, GF28_Element *outputvec,
        nqueries_t num_queries, dbsize_t num_blocks, dbsize_t words_per_block, 
        const vector<dbsize_t> &unsynchronized_files) {
    if (num_queries == 2) {
	const GF28_Element *block = data;
	const GF28_Element *inp1 = inputvec;
	const GF28_Element *inp2 = inputvec + num_blocks;
	for (unsigned int j = 0; j < num_blocks; ++j) {
	    const GF28_Element *multrow1 = GF28_mult_table[*(inp1++)];
	    const GF28_Element *multrow2 = GF28_mult_table[*(inp2++)];
	    const GF28_Element *blockc = block;
	    GF28_Element *oc1 = outputvec;
	    GF28_Element *oc2 = outputvec + words_per_block;
	    GF28_Element *oc1_end = oc1 + (words_per_block & ~7);
	    while (oc1 < oc1_end) {
		GF28_Element v1 = *(blockc++);
		uint64_t accum1 = (uint64_t) multrow1[v1];
		uint64_t accum2 = (uint64_t) multrow2[v1];
		GF28_Element v2 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v2] << 8;
		accum2 |= (uint64_t) multrow2[v2] << 8;
		GF28_Element v3 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v3] << 16;
		accum2 |= (uint64_t) multrow2[v3] << 16;
		GF28_Element v4 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v4] << 24;
		accum2 |= (uint64_t) multrow2[v4] << 24;
		GF28_Element v5 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v5] << 32;
		accum2 |= (uint64_t) multrow2[v5] << 32;
		GF28_Element v6 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v6] << 40;
		accum2 |= (uint64_t) multrow2[v6] << 40;
		GF28_Element v7 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v7] << 48;
		accum2 |= (uint64_t) multrow2[v7] << 48;
		GF28_Element v8 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v8] << 56;
		accum2 |= (uint64_t) multrow2[v8] << 56;
		*((uint64_t *) oc1) ^= accum1;
		*((uint64_t *) oc2) ^= accum2;
		oc1+=8;
		oc2+=8;
	    }
	    for (unsigned int c = 0; c < (words_per_block & 7); ++c) {
		GF28_Element v = *(blockc++);
		*(oc1++) ^= multrow1[v];
		*(oc2++) ^= multrow2[v];
	    }
	    block += words_per_block;
	}
    } else {
	const GF28_Element *block = data;
	const GF28_Element **inpv = new const GF28_Element*[num_queries];
	const GF28_Element **multrowv = new const GF28_Element*[num_queries];
	GF28_Element **ocv = new GF28_Element*[num_queries];
	for (unsigned int q=0; q<num_queries; ++q) {
	    inpv[q] = inputvec + q*num_blocks;
	}
	for (unsigned int j = 0; j < num_blocks; ++j) {
	    for (unsigned int q=0; q<num_queries; ++q) {
		multrowv[q] = GF28_mult_table[*(inpv[q]++)];
	    }
	    const GF28_Element *blockc = block;
	    GF28_Element *ocv0 = outputvec;
	    for (unsigned int q=1; q<num_queries; ++q) {
		ocv[q] = outputvec + q*words_per_block;
	    }
	    GF28_Element *oc_end = ocv0 + (words_per_block & ~7);
	    while (ocv0 < oc_end) {
		GF28_Element v1 = *(blockc++);
		GF28_Element v2 = *(blockc++);
		GF28_Element v3 = *(blockc++);
		GF28_Element v4 = *(blockc++);
		GF28_Element v5 = *(blockc++);
		GF28_Element v6 = *(blockc++);
		GF28_Element v7 = *(blockc++);
		GF28_Element v8 = *(blockc++);
		uint64_t accum = (uint64_t) multrowv[0][v1];
		accum |= (uint64_t) multrowv[0][v2] << 8;
		accum |= (uint64_t) multrowv[0][v3] << 16;
		accum |= (uint64_t) multrowv[0][v4] << 24;
		accum |= (uint64_t) multrowv[0][v5] << 32;
		accum |= (uint64_t) multrowv[0][v6] << 40;
		accum |= (uint64_t) multrowv[0][v7] << 48;
		accum |= (uint64_t) multrowv[0][v8] << 56;
		*((uint64_t *) ocv0) ^= accum;
		ocv0 += 8;
		for (unsigned int q=1; q<num_queries; ++q) {
		    uint64_t accum = (uint64_t) multrowv[q][v1];
		    accum |= (uint64_t) multrowv[q][v2] << 8;
		    accum |= (uint64_t) multrowv[q][v3] << 16;
		    accum |= (uint64_t) multrowv[q][v4] << 24;
		    accum |= (uint64_t) multrowv[q][v5] << 32;
		    accum |= (uint64_t) multrowv[q][v6] << 40;
		    accum |= (uint64_t) multrowv[q][v7] << 48;
		    accum |= (uint64_t) multrowv[q][v8] << 56;
		    *((uint64_t *) ocv[q]) ^= accum;
		    ocv[q] += 8;
		}
	    }
	    for (unsigned int c = 0; c < (words_per_block & 7); ++c) {
		GF28_Element v = *(blockc++);
		*(ocv0++) ^= multrowv[0][v];
		for (unsigned int q=1; q<num_queries; ++q) {
		    *(ocv[q]++) ^= multrowv[q][v];
		}
	    }
	    block += words_per_block;
	}
	delete[] inpv;
	delete[] multrowv;
	delete[] ocv;
    }
}

template <>
inline void compute_outputvec_multi<GF216_Element>(
        const GF216_Element *data, GF216_Element *inputvec, GF216_Element *outputvec,
        nqueries_t num_queries, dbsize_t num_blocks, dbsize_t words_per_block, const vector<dbsize_t> &unsynchronized_files) {
    for (nqueries_t q = 0; q < num_queries; q++) {
        compute_outputvec_single<GF216_Element>(data, inputvec + q*num_blocks, 
                outputvec + q*words_per_block, num_blocks, words_per_block, unsynchronized_files);
    }
}


template<dbsize_t arrSize>
void PercyServer::compute_outputvec_sync( 
        const GF216_Element *data, 
        std::vector<std::array<GF216_Element, arrSize> > &outputvec,
        dbsize_t num_blocks, 
        dbsize_t words_per_block,
        dbsize_t max_unsynchronized,
        dbsize_t num_rows,
        GF216_Element q_x ) {
    
    const GF216_Element *block = data;
    
    
    for (dbsize_t j = 0; j < num_blocks; ++j) {
        // std::cerr << "Multiplying everything by " << q_x << std::endl;
        if (std::find(unsynchronized_files.begin(), unsynchronized_files.end(), j)==unsynchronized_files.end()){
            for (GF216_Element k = 0; k < DEGREE; ++k) {
                GF216_Element bin;
                // Which bin should we put this in?
                switch (num_rows/NUM_RATIOS) {
                    case 4:
                        bin = GF216_pulse_mtx_4bins[j][k];
                        break;
                    case 5:
                        bin = GF216_pulse_mtx_5bins[j][k];
                        break;
                    case 6: 
                        bin = GF216_pulse_mtx_6bins[j][k];
                        break;
                    case 7:
                        bin = GF216_pulse_mtx_7bins[j][k];
                        break;
                    case 8: 
                        bin = GF216_pulse_mtx_8bins[j][k];
                        break;
                    case 9:
                        bin = GF216_pulse_mtx_9bins[j][k];
                        break;
                    case 10: 
                        bin = GF216_pulse_mtx_10bins[j][k];
                        break;
                    case 11:
                        bin = GF216_pulse_mtx_11bins[j][k];
                        break;
                    case 12: 
                        bin = GF216_pulse_mtx_12bins[j][k];
                        break;
                    default:
                        bin = GF216_pulse_mtx_1000bins[j][k] % (num_rows/NUM_RATIOS);
                        break;
                }
                
                
                for (GF216_Element m = 0; m < NUM_RATIOS; ++m) {
                    // Multiply the element by the product of alpha^(m*n)*q[j], where
                    //      alpha is a primitive element of GF2^16
                    GF216_Element inpv_j = multiply_GF2E( q_x,
                                            GF216_exp_table[multiply_GF2E(m,(GF216_Element) j)]);
                    if (inpv_j != 0) {
                        // Multiply the next entire block by j, and then xor it with the current results
                        // we simply need to choose which outputvec we're going to add it to.
                        const GF216_Element *blockc = block;
                        GF216_Element log_j = GF216_log_table[inpv_j];
                        const GF216_Element *start = GF216_exp_table + log_j;
                        GF216_Element *oc = &outputvec[bin*NUM_RATIOS+m][0];
                        GF216_Element *oc_end = oc + (words_per_block & ~3);
                        GF216_Element block_c;
                        while(oc < oc_end) {
                            uint64_t accum = 0;
                            block_c = *(blockc++);
                            if (block_c != 0) {
                                GF216_Element log_c = GF216_log_table[block_c];
                                accum |= (uint64_t) start[log_c];
                            }
                            block_c = *(blockc++);
                            if (block_c != 0) {
                                GF216_Element log_c = GF216_log_table[block_c];
                                accum |= (uint64_t) start[log_c] << 16;
                            }
                            block_c = *(blockc++);
                            if (block_c != 0) {
                                GF216_Element log_c = GF216_log_table[block_c];
                                accum |= (uint64_t) start[log_c] << 32;
                            }
                            block_c = *(blockc++);
                            if (block_c != 0) {
                                GF216_Element log_c = GF216_log_table[block_c];
                                accum |= (uint64_t) start[log_c] << 48;
                            }
                            *((uint64_t *) oc) ^= accum;
                            oc+=4;
                        }
                        for (dbsize_t c = 0; c < (words_per_block & 3); ++c, ++oc) {
                            block_c = *(blockc++);
                            if (block_c != 0) {
                                GF216_Element log_c = GF216_log_table[block_c];
                                *oc ^= start[log_c];
                            }
                        }
                    }
                
                }
            }
        }
        // std::cerr << "sync output contains " << outputvec[j][20] << outputvec[j][21] << std::endl;
        block += words_per_block;
    }
}
    
template <typename GF2E_Element>
bool PercyServer::handle_request_GF2E(PercyServerParams &params
	, std::istream &is, std::ostream &os) 
{
    if (is.eof()) {
        return false;
    }

    // Read some values from the params
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t num_blocks = params.num_blocks();

    // Read the number of queries
    unsigned char nq[2];
    is.read((char *)nq, 2);
    if (is.eof()) {
        return false;
    }

    nqueries_t num_queries = (nq[0] << 8) | nq[1];
    if (num_queries == 0) {
        std::cerr << "No more queries\n";
        return false;
    }

    // For each query, read the input vector, which is a sequence of
    // num_blocks entries, each of length sizeof(GF2E_Element) bytes
    GF2E_Element *input = new GF2E_Element[num_queries*num_blocks];
    GF2E_Element *output = new GF2E_Element[num_queries*words_per_block];
    memset(output, '\0', num_queries*words_per_block*sizeof(GF2E_Element));
    is.read((char *)input, num_queries*num_blocks*sizeof(GF2E_Element));
    if (is.eof()) {
        delete[] input;
        delete[] output;
        return false;
    }

    const MemoryDataStore* ds = static_cast<MemoryDataStore*>(datastore);
    const GF2E_Element *data = (const GF2E_Element*)(ds->get_data());

    // Compute the output vector and send it back to the client

    //struct timeval ts, te;
    //gettimeofday(&ts, NULL);
    if (num_queries > 1) {
        compute_outputvec_multi<GF2E_Element>(data, input, output,
            num_queries, num_blocks, words_per_block, unsynchronized_files);
    } else {
        compute_outputvec_single<GF2E_Element>(data, input, output,
            num_blocks, words_per_block, unsynchronized_files);
    }

    //gettimeofday(&te, NULL);
    //int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
    //fprintf(stderr, "%d.%3d msec computation\n", td/1000, td%1000);

    // If the server is Byzantine, give wrong output
    if (byzantine) {
        for (nqueries_t q = 0; q < num_queries; q++) {
            for (dbsize_t i = 0; i < words_per_block; i++) {
                output[q*words_per_block + i]++;
            }
        }
    }

    // Send the output vector to the client via the ostream
    os.write((char *)output,
		num_queries*words_per_block*sizeof(GF2E_Element));
    os.flush();
    delete[] input;
    delete[] output;

    return true;
}

template <typename GF2E_Element>
bool PercyServer::handle_request_RS_Sync(PercyServerParams &params, std::istream &is,
        std::ostream &os)
{
    if(is.eof()) {
        return false;
    }
    // Read some values from the params
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t num_blocks = params.num_blocks();
    // dbsize_t num_bytes = params.num_blocks() / 8;
    dbsize_t max_unsynchronized = params.max_unsynchronized();

    // Read the number of queries
    unsigned char nq[2];
    is.read((char *)nq, 2);
    if (is.eof()) {
        return false;
    }
    nqueries_t num_queries = (nq[0] << 8) | nq[1];
    if (num_queries == 0) {
        std::cerr << "No queries were requested! \n";
        return false;
    }

    // For each query, read the input vector, which is a sequence of
    // num_blocks entries, each of length sizeof(GF2E_Element) bytes
    GF2E_Element *input = new GF2E_Element[num_queries*num_blocks];
    GF2E_Element *output = new GF2E_Element[num_queries*words_per_block];
    memset(output, '\0', num_queries*words_per_block*sizeof(GF2E_Element));
    is.read((char *)input, num_queries*num_blocks*sizeof(GF2E_Element));
    if (is.eof()) {
        delete[] input;
        delete[] output;
        std::cerr << "ERROR: Not enough content in the input stream!\n";
        return false;
    }

    const MemoryDataStore* ds = static_cast<MemoryDataStore*>(datastore);
    const GF2E_Element *data = (const GF2E_Element*)(ds->get_data());

    // Compute the output vector and send it back to the client

    //struct timeval ts, te;
    //gettimeofday(&ts, NULL);
    
    dbsize_t starting_point = 0;
    if (server_unsynchronized) {
        starting_point = max_unsynchronized;
    }
    if (num_queries > 1) {
        compute_outputvec_multi<GF2E_Element>(data, input, output,
            num_queries, num_blocks, words_per_block, unsynchronized_files);
    } else {
        compute_outputvec_single<GF2E_Element>(data, input, output,
            num_blocks, words_per_block, unsynchronized_files);
    }

    //gettimeofday(&te, NULL);
    //int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
    //fprintf(stderr, "%d.%3d msec computation\n", td/1000, td%1000);

    // If the server is Byzantine, give wrong output
    if (byzantine) {
        for (nqueries_t q = 0; q < num_queries; q++) {
            for (dbsize_t i = 0; i < words_per_block; i++) {
            output[q*words_per_block + i]++;
            }
        }
    }

    // Send the output vector to the client via the ostream
    os.write((char *)output,
        num_queries*words_per_block*sizeof(GF2E_Element));
    os.flush();
    delete[] input;
    delete[] output;

    return true;

}

template <typename GF2E_Element>
bool PercyServer::handle_sync_request_RS_Sync(PercyServerParams &params, std::istream &is,
        std::ostream &os)
{
    // Read some values from the params
    // dbsize_t words_per_block = params.words_per_block();
    
    if ( WORDS_PER_BLOCK != params.words_per_block() ) {
        std::cerr << "The words per block are incorrect!\n";
        return false;
    }
    // dbsize_t num_bytes = params.num_blocks() / 8;
    dbsize_t num_blocks = params.num_blocks();
    dbsize_t max_unsynchronized = params.max_unsynchronized();
    dbsize_t num_bins = params.num_bins();
    
    // How many rows do you want to consider?
    dbsize_t num_rows = num_bins * NUM_RATIOS;
    // For each query, read the input vector, which is a sequence of
    // num_blocks entries, each of length sizeof(GF2E_Element) bytes
    // output will have multiple results (one for each bin) 
    std::vector<std::array<GF2E_Element,WORDS_PER_BLOCK> > output(num_rows);
    GF2E_Element filler = 0;
    for (dbsize_t i=0; i<num_rows; i++) {
        output[i].fill(filler);
    }
    // GF2E_Element *output = new GF2E_Element[words_per_block];
    
    // Read in the value of X to evaluate the input polynomial at!
    unsigned char input[sizeof(GF2E_Element)];
    is.read((char *)input, sizeof(GF2E_Element));
    GF2E_Element q_x = (input[0] << 8) | input[1];

    const MemoryDataStore* ds = static_cast<MemoryDataStore*>(datastore);
    const GF2E_Element *data = (const GF2E_Element*)(ds->get_data());

    // Compute the necessary hashes, and return them to the client
    compute_outputvec_sync(data, output, num_blocks, WORDS_PER_BLOCK, 
        max_unsynchronized, num_rows, q_x);

    // If the server is Byzantine, give wrong output
    if (byzantine) {
        for (dbsize_t k = 0; k < num_rows; k++) {
            for (dbsize_t i = 0; i < WORDS_PER_BLOCK; i++) {
                output[k][i]++;
            }
        }
    }
    
    // Send the output vector to the client via the ostream
    for (dbsize_t k = 0; k < num_rows; k++) {
        os.write((char *) &output[k], WORDS_PER_BLOCK * sizeof(GF2E_Element));
    }
    os.flush();

    return true;
}

#endif
