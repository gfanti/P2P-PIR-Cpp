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

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>  
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <vec_ZZ_p.h>
#include <vector>
#include "percytypes.h"
#include "percyserver.h"
#include "percyserver_impl.h"
#include "xor.h"

#ifdef SPIR_SUPPORT
/*
#include "../PolyCommit/PolyCommitCommon.h"
#include "pspir_crypt.h"
*/
#include "spirserver.h"
#endif

// Handle a request.
bool PercyServer::handle_request(PercyServerParams &params, std::istream &is,
        std::ostream &os)
{
    if (datastore == 0) {
        std::cerr << "No datastore!\n";
        return false;
    }
    if (params.get_mode() == MODE_GF28) {
        //return handle_request_GF28(params, is, os);
        std::cerr << "mode is gf2^8\n";
        return handle_request_GF2E<GF28_Element>(params, is, os);
    }
    else if (params.get_mode() == MODE_GF216) {
        std::cerr << "mode is gf2^16\n";
        return handle_request_GF2E<GF216_Element>(params, is, os);
    }
    else if (params.get_mode() == MODE_ZZ_P) {
        std::cerr << "mode is zz_p\n";
        return handle_request_ZZ_p(params, is, os);
    }
    else if (params.get_mode() == MODE_CHOR) {
        std::cerr << "mode is chor\n";
        return handle_request_Chor(params, is, os);
    }
    else if (params.get_mode() == MODE_RS_SYNC) {
        // First it should handle the synchronization errors
        if (!synchronized) {
            set_synchronized();
            return handle_sync_request_RS_Sync<GF216_Element>(params, is, os);
        }
        // Then it can handle the true PIR request
        return handle_request_RS_Sync<GF216_Element>(params, is, os);
    }
    // else if (params.get_mode() == MODE_PULSE_SYNC) {
        // return handle_request_PULSE_Sync(params, is, os);
    // }
    else {
        std::cerr << "Unsupported mode selected in PercyParams." << std::endl;
        exit(1);
    }

    return false;
}

// Handle a request.
bool PercyServer::handle_request_ZZ_p(PercyServerParams &params, std::istream &is,
        std::ostream &os)
{
    // Set the appropriate modulus
    unsigned int modulus_bytes;
    if (params.hybrid()) {
        params.mod_modulussq();
        modulus_bytes = params.modulussq_bytes();
    } else {
        params.mod_modulus();
        modulus_bytes = params.modulus_bytes();
    }

    // Read some values from the params
    bool hybrid_protection = params.hybrid();
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t num_blocks = params.num_blocks();

    // Read the number of queries
    //std::cerr << "Receiving number of queries from client...";
    unsigned char nq[2];
    is.read((char *)nq, 2);
    if (is.eof()) {
        return false;
    }
    nqueries_t num_queries = (nq[0] << 8) | nq[1];
    //std::cerr << "done" << std::endl;

#ifdef USE_W160_OPT
    dbsize_t bytes_per_word = params.bytes_per_word();
    // Use a faster method if w160 is selected and there is only one query
    bool use_w160_method = bytes_per_word == 20 && modulus_bytes == 21
        && num_queries == 1 && num_blocks < (1ULL<<32);
#endif

#ifdef SPIR_SUPPORT
    PolyCommitParams * pcparamsp = params.get_pcparamsp();

    //Save ZZ_p context
    ZZ_pContext savectx;
    if (params.spir()) {
        savectx.save();
        ZZ_p::init(pcparamsp->get_order());
    }

    // Read in the servers' shared secret key from a file. This too should probably 
    // be read in from somewhere else. Later we combine this with a hash to seed 
    // NTL's PRNG.
    unsigned char prng_seed[24];
    ifstream infile2("server_key.txt");
    if(!infile2.is_open()) {
        std::cerr << "Error: Cannot open server key file." << std::endl;
        exit(1);
    }
    infile2.read((char*)prng_seed, 16);
    infile2.close();

    SPIRServerQuery * spir_queries = new SPIRServerQuery[num_queries];
#endif

    //Read in the input
    vec_ZZ_p * inputvector = new vec_ZZ_p[num_queries]; //Used for non-w160, and w160 with SPIR
#ifdef USE_W160_OPT
    unsigned char *indata = NULL; //Used for w160 without SPIR
    if (use_w160_method) {
        indata = new unsigned char[num_blocks*modulus_bytes];
    }
#endif
    for(nqueries_t q = 0; q < num_queries; ++q) {
        //std::cerr << "Receiving query " << (q+1) << " of " << num_queries << " from client...";
#ifdef SPIR_SUPPORT
        if (params.spir()) {
	    bool ret = spir_queries[q].read_spir_input(params, is);
	    if (!ret) {
		delete[] spir_queries;
		return false;
	    }
        }
#endif
#ifdef USE_W160_OPT
        //If the w160 method is used, read into the indata array, and then if SPIR is 
        //being used, copy this data into the inputvector as well
        if (use_w160_method) {
            is.read((char *)indata, num_blocks*modulus_bytes);
            if (is.eof()) {
#ifdef SPIR_SUPPORT
		delete[] spir_queries;
#endif
                return false;
            }
            if (params.spir()) {
                inputvector[q].SetLength(num_blocks);
                for(dbsize_t i = 0; i < num_blocks; ++i) {
                    ZZ inputz;
                    ZZFromBytes(inputz, indata + i*modulus_bytes, modulus_bytes);
                    inputvector[q][i] = to_ZZ_p(inputz);
                }
            }
        } else {
#endif
            // Otherwise:
            // For each query, read the input vector, which is a sequence of
            // num_blocks entries, each of length modulus_bytes.
            inputvector[q].SetLength(num_blocks);
            unsigned char * inputdata = new unsigned char[modulus_bytes];
            for (dbsize_t i = 0; i < num_blocks; ++i) {
                memset(inputdata, 0, modulus_bytes);
                ZZ inputz;
                is.read((char *)inputdata, modulus_bytes);
                if (is.eof()) {
#ifdef SPIR_SUPPORT
		    delete[] spir_queries;
#endif
                    return false;
                }
                ZZFromBytes(inputz, inputdata, modulus_bytes);
                inputvector[q][i] = to_ZZ_p(inputz);
            }
            delete[] inputdata;
#ifdef USE_W160_OPT
        }
#endif
    }

#ifdef SPIR_SUPPORT
    if (params.spir()) {
        for(unsigned int q = 0; q < num_queries; ++q) {
	    bool ret = spir_queries[q].query_verification(params, inputvector[q]);
	    if (!ret) {
		delete[] spir_queries;
		return false;
	    }
        }

        for(unsigned int q = 0; q < num_queries; ++q) {
	    spir_queries[q].init_randomization(params, prng_seed);
        }

        savectx.restore();
    }
#endif

#ifdef USE_W160_OPT
    // Compute the output vector and send it back to the client
    if (use_w160_method) {
        // Get a pointer to the database data
        const MemoryDataStore* ds = static_cast<MemoryDataStore*>(datastore);
        const unsigned char *dbdata = ds->get_data();

        // Initialize the output to 0.  Each word of the output will
        // temporarily be stored as 5 128-bit values, and we'll do the
        // modular reduction right at the end.  The 5 128-bit values
        // represent the number
        // a[0] + a[1]<<56 + a[2]<<112 + a[3]<<168 + a[4]<<224.
        __uint128_t out128s[5*words_per_block];
        memset(out128s, 0, sizeof(out128s));

        // Process each block of the database
        unsigned char *multdata = indata;
        for (dbsize_t j = 0; j < num_blocks;
                multdata += modulus_bytes, ++j) {
            __uint128_t multiplier[3];
            multiplier[0] = (*(uint64_t*)multdata) & 0x00ffffffffffffffULL;
            multiplier[1] = (*(uint64_t*)(multdata+7)) & 0x00ffffffffffffffULL;
            multiplier[2] = (*(uint64_t*)(multdata+14)) & 0x00ffffffffffffffULL;

            __uint128_t *outval = out128s;

            // Process each word in the block
            for (dbsize_t c = 0; c < words_per_block-1;
                    dbdata += bytes_per_word, outval += 5, ++c) {
                uint64_t dbval[3];
                dbval[0] = (*(uint64_t*)dbdata) & 0x00ffffffffffffffULL;
                dbval[1] = (*(uint64_t*)(dbdata+7)) & 0x00ffffffffffffffULL;
                dbval[2] = (*(uint64_t*)(dbdata+14)) & 0x0000ffffffffffffULL;

                // outval += multiplier * dbval
                outval[0] += multiplier[0] * dbval[0];
                outval[1] += multiplier[0] * dbval[1]
                    + multiplier[1] * dbval[0];
                outval[2] += multiplier[0] * dbval[2]
                    + multiplier[1] * dbval[1]
                    + multiplier[2] * dbval[0];
                outval[3] += multiplier[1] * dbval[2]
                    + multiplier[2] * dbval[1];
                outval[4] += multiplier[2] * dbval[2];

                // Do the carries
                outval[1] += (outval[0]>>56);
                outval[0] &= 0x00ffffffffffffffULL;
                outval[2] += (outval[1]>>56);
                outval[1] &= 0x00ffffffffffffffULL;
                outval[3] += (outval[2]>>56);
                outval[2] &= 0x00ffffffffffffffULL;
                outval[4] += (outval[3]>>56);
                outval[3] &= 0x00ffffffffffffffULL;
            }
            // Last block: we may not be able to read past the end of
            // the database
            {
                uint64_t dbval[3];
                dbval[0] = (*(uint64_t*)dbdata) & 0x00ffffffffffffffULL;
                dbval[1] = (*(uint64_t*)(dbdata+7)) & 0x00ffffffffffffffULL;
                dbval[2] = 0;
                memmove(dbval+2, dbdata+14, 6);

                // outval += multiplier * dbval
                outval[0] += multiplier[0] * dbval[0];
                outval[1] += multiplier[0] * dbval[1]
                    + multiplier[1] * dbval[0];
                outval[2] += multiplier[0] * dbval[2]
                    + multiplier[1] * dbval[1]
                    + multiplier[2] * dbval[0];
                outval[3] += multiplier[1] * dbval[2]
                    + multiplier[2] * dbval[1];
                outval[4] += multiplier[2] * dbval[2];

                // Do the carries
                outval[1] += (outval[0]>>56);
                outval[0] &= 0x00ffffffffffffffULL;
                outval[2] += (outval[1]>>56);
                outval[1] &= 0x00ffffffffffffffULL;
                outval[3] += (outval[2]>>56);
                outval[2] &= 0x00ffffffffffffffULL;
                outval[4] += (outval[3]>>56);
                outval[3] &= 0x00ffffffffffffffULL;

                dbdata += bytes_per_word;
            }
        }

        // Now convert the 5 128-bit values in each output word into a
        // single ZZ, reduce it mod the modulus, and output the
        // resulting 21-byte value
        __uint128_t *outval = out128s;
        //std::cerr << "Sending response to client...";
	vec_ZZ_p reponse;
	response.SetLength(words_per_block);
        for (dbsize_t c = 0; c < words_per_block; outval += 5, ++c) {
            unsigned char valbuf[4*7+16];
            memmove(valbuf, (unsigned char*)(outval), 7);
            memmove(valbuf+7, (unsigned char*)(outval+1), 7);
            memmove(valbuf+14, (unsigned char*)(outval+2), 7);
            memmove(valbuf+21, (unsigned char*)(outval+3), 7);
            memmove(valbuf+28, (unsigned char*)(outval+4), 16);
            response[c] = to_ZZ_p(ZZFromBytes(valbuf,sizeof(valbuf)));
	}

#ifdef SPIR_SUPPORT
	if (params.spir()) {
	    spir_queries[0].randomize_response(response);
	}
#endif

	unsigned char outmodbuf[modulus_bytes];
	for (dbsize_t c = 0; c < words_per_block; ++c) {
            BytesFromZZ(outmodbuf, rep(response[c]), modulus_bytes);
            os.write((char *)outmodbuf, modulus_bytes);
        }
    } else {
#endif
	
	vec_ZZ_p * responses = new vec_ZZ_p[num_queries];
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    responses[q].SetLength(words_per_block);
	}
	//ZZ_p responses[words_per_block][num_queries];
        for (dbsize_t c = 0; c < words_per_block; ++c) {
            ZZ_p * reswrds = new ZZ_p[num_queries];
            compute_one(reswrds, hybrid_protection, num_blocks,
                    num_queries, inputvector, c);
	    for (unsigned int q = 0; q < num_queries; ++q) {
		responses[q][c] = reswrds[q];
	    }
	    delete[] reswrds;
        }

#ifdef SPIR_SUPPORT
	if (params.spir()) {
	    for (nqueries_t q = 0; q < num_queries; ++q) {
		spir_queries[q].randomize_response(responses[q]);
	    }
	}
#endif

	unsigned char * bytes = new unsigned char[modulus_bytes];
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    for (dbsize_t c = 0; c < words_per_block; ++c) {
		BytesFromZZ(bytes, rep(responses[q][c]), modulus_bytes);
		os.write((char *)bytes, modulus_bytes);
	    }
	}

	delete[] bytes;
	delete[] responses;
#ifdef USE_W160_OPT
    }
#endif
    os.flush();

#ifdef USE_W160_OPT
    if (use_w160_method) {
        delete[] indata;
    }
#endif

#ifdef SPIR_SUPPORT
    delete[] spir_queries;
#endif
    delete[] inputvector;

    return true;
}

bool PercyServer::handle_request_Chor(PercyServerParams &params, std::istream &is,
        std::ostream &os)
{
    if(!is.eof()) {
        // Read some values from the params
        dbsize_t words_per_block = params.words_per_block();
        dbsize_t num_bytes = params.num_blocks() / 8;
        nqueries_t q;

        // Read the number of queries
        unsigned char nq[2];
        is.read((char *)nq, 2);
        if (is.eof()) {
            return false;
        }
        nqueries_t num_queries = (nq[0] << 8) | nq[1];

        // For each query, read the input vector, which is a sequence of
        // num_blocks/8 entries, each of length 1 byte
        unsigned char *inputvector = new unsigned char[num_queries*num_bytes];
        unsigned char *outputvector = new unsigned char[num_queries*words_per_block];
        memset(outputvector, '\0', num_queries*words_per_block);

        is.read((char *)inputvector, num_queries*num_bytes);
        if (is.eof()) {
	    delete[] inputvector;
	    delete[] outputvector;
            return false;
        }

        /*
        std::cerr << "Server " << ": ";
        printBS_debug(inputvector, num_blocks/8);
        std::cerr << std::endl;
        */

        const MemoryDataStore* ds = static_cast<MemoryDataStore*>(datastore);

        const unsigned char *data = (const unsigned char*)(ds->get_data());

        // Compute the output vector and send it back to the client
        struct timeval ts, te;
        gettimeofday(&ts, NULL);
        for (q = 0; q < num_queries; q++) {
            for (unsigned int i = 0; i < num_bytes; i++) {
                unsigned char query_byte = inputvector[q*num_bytes + i];
                if (query_byte & 128)
                    XOR_equal(outputvector + q*words_per_block,
                            data + 8*i*words_per_block, words_per_block);
                if (query_byte & 64)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+1)*words_per_block, words_per_block);
                if (query_byte & 32)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+2)*words_per_block, words_per_block);
                if (query_byte & 16)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+3)*words_per_block, words_per_block);
                if (query_byte & 8)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+4)*words_per_block, words_per_block);
                if (query_byte & 4)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+5)*words_per_block, words_per_block);
                if (query_byte & 2)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+6)*words_per_block, words_per_block);
                if (query_byte & 1)
                    XOR_equal(outputvector + q*words_per_block,
                            data + (8*i+7)*words_per_block, words_per_block);
            }
        }

        gettimeofday(&te, NULL);
        int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
        fprintf(stderr, "Chor: %d.%03d msec computation\n", td/1000, td%1000);

        if (byzantine) {
            for (dbsize_t i=0; i<num_queries*words_per_block; ++i) {
                outputvector[i]++;
            }
        }

        /*
        std::cerr << "Server: ";
        printBS_debug(outputvector+(num_queries-1)*words_per_block, words_per_block);
        std::cerr << std::endl;
        */

        os.write((char *)outputvector, num_queries*words_per_block);
        os.flush();
	delete[] inputvector;
	delete[] outputvector;
        return true;
    } else {
        return false;
    }
}

// Tell the server that it should hold unsynchronized files
void PercyServer::set_server_unsynchronized(PercyServerParams &params) {
    server_unsynchronized = true; 
    srand (time(NULL));
    // decide which files will be unsynchronized, and print them out
    for (dbsize_t i=0; i<params.max_unsynchronized(); i++) {
        // pick a random number in num_blocks
        
        dbsize_t idx = std::rand() % params.num_blocks();
        unsynchronized_files.push_back(idx);
        std::cerr << "Unsynchronized file at " << idx << std::endl;
    }
}

void PercyServer::compute_one(ZZ_p *value, bool hybrid_protection,
        dbsize_t num_blocks, nqueries_t num_queries,
        const vec_ZZ_p *inputvector, dbsize_t c)
{
    for(nqueries_t q = 0; q < num_queries; ++q) {
        value[q] = hybrid_protection ? 1 : 0;
    }

    for (dbsize_t j = 0; j < num_blocks; ++j) {
        // The cth word of the jth block
        ZZ wrd = datastore->get_word(c, j);

        if (hybrid_protection) {
            for(nqueries_t q = 0; q < num_queries; ++q) {
                value[q] = value[q] * power(inputvector[q][j], wrd);
            }
        } else {
            for(nqueries_t q = 0; q < num_queries; ++q) {
                value[q] = value[q] + inputvector[q][j] * to_ZZ_p(wrd);
            }
        }
    }

    if (byzantine) {
        // Produce a *consistent* incorrect value for maximal client
        // confusion.
        for(nqueries_t q = 0; q < num_queries; ++q) {
            value[q] += 1;
        }
    }
}
