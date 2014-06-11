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

#include <string.h>
#include <fstream>
#include "percyio.h"
#include "distserver.h"
#include "xor.h"

#ifdef SPIR_SUPPORT
#include "spirserver.h"
#endif

bool PercyMasterServer::handle_request(PercyServerParams &params, std::istream &is, 
	std::ostream &os)
{
    PercyDistServerParams * dparams = static_cast<PercyDistServerParams*>(&params);
    if (dparams == NULL) {
	std::cerr << "params is not a PercyDistServerParams\n";
	return false;
    }

    // Get parameters
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t num_blocks = params.num_blocks();
    dbsize_t bytes_per_word = params.bytes_per_word();
    dbsize_t vsplit = dparams->get_vsplit();
    dbsize_t hsplit = dparams->get_hsplit();
    dbsize_t num_subdb = vsplit * hsplit;
    dbsize_t subdb_words_per_block = words_per_block / hsplit;
    dbsize_t subdb_num_blocks = num_blocks / vsplit;
    std::vector<std::iostream*>& workerio = dparams->get_workerio();
    PercyMode mode = params.get_mode();

    // Get modulus bytes for ZZ_p
    if (mode == MODE_ZZ_P) {
	if (params.hybrid()) {
	    params.mod_modulussq();
	    bytes_per_word = params.modulussq_bytes();
	} else {
	    params.mod_modulus();
	    bytes_per_word = params.modulus_bytes();
	}
    }

    // Change num_blocks for Chor
    if (mode == MODE_CHOR) {
	num_blocks /= 8;
	subdb_num_blocks /= 8;
    }

    // Read the number of queries
    unsigned char nq[2];
    is.read((char*)nq, 2);
    if (is.eof()) {
	return false;
    }
    nqueries_t num_queries = (nq[0] << 8) | nq[1];

    // Send the number of queries to each worker
    for (dbsize_t i = 0; i < num_subdb; ++i) {
	workerio[i]->write((char*)nq, 2);
    }

    // Check if done
    if (num_queries == 0) {
	return false;
    }

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
    std::ifstream infile2("server_key.txt");
    if(!infile2.is_open()) {
        std::cerr << "Error: Cannot open server key file." << std::endl;
        exit(1);
    }
    infile2.read((char*)prng_seed, 16);
    infile2.close();

    SPIRServerQuery * spir_queries = new SPIRServerQuery[num_queries];
    vec_ZZ_p * inputvector = new vec_ZZ_p[num_queries];
#endif

    // Read the queries
    dbsize_t query_size = num_blocks * bytes_per_word;
    unsigned char * queries = new unsigned char[num_queries * query_size];
    for (nqueries_t q = 0; q < num_queries; ++q) {
#ifdef SPIR_SUPPORT
        if (params.spir()) {
	    bool ret = spir_queries[q].read_spir_input(params, is);
	    if (!ret) {
		delete[] spir_queries;
		delete[] inputvector;
		delete[] queries;
		return false;
	    }
        }
#endif
	is.read((char*)(queries + (q * query_size)), query_size);
	if (is.eof()) {
#ifdef SPIR_SUPPORT
	    delete[] spir_queries;
	    delete[] inputvector;
#endif
	    delete[] queries;
	    return false;
	}
#ifdef SPIR_SUPPORT
	if (params.spir()) {
	    inputvector[q].SetLength(num_blocks);
	    for (dbsize_t i = 0; i < num_blocks; ++i) {
		ZZ inputz;
		ZZFromBytes(inputz, (queries + (query_size * q) + (i * bytes_per_word)), 
			bytes_per_word);
		inputvector[q][i] = to_ZZ_p(inputz);
	    }
	}
#endif
    }

#ifdef SPIR_SUPPORT
    if (params.spir()) {
        for(unsigned int q = 0; q < num_queries; ++q) {
	    bool ret = spir_queries[q].query_verification(params, inputvector[q]);
	    if (!ret) {
		delete[] spir_queries;
		delete[] inputvector;
		delete[] queries;
		return false;
	    }
        }

        for(unsigned int q = 0; q < num_queries; ++q) {
	    spir_queries[q].init_randomization(params, prng_seed);
        }

        savectx.restore();
    }
    delete[] inputvector;
#endif

    // Send the partial queries
    dbsize_t worker_query_size = subdb_num_blocks * bytes_per_word;
    dbsize_t qoffset = 0;
    for (nqueries_t k = 0; k < num_queries; ++k) {
	for (dbsize_t i = 0; i < vsplit; ++i) {
	    for (dbsize_t j = 0; j < hsplit; ++j) {
		std::iostream& wio = *(workerio[(hsplit * i) + j]);
		wio.write((char*)(queries + qoffset), worker_query_size);
		wio.flush();
	    }
	    qoffset += worker_query_size;
	}
    }
    delete[] queries;

    // Receive the replies
    unsigned char * responses = new unsigned char[num_queries * words_per_block * bytes_per_word];
    memset(responses, 0, num_queries * words_per_block * bytes_per_word);

    vec_ZZ_p responses_ZZ_p;
    if (mode == MODE_ZZ_P) {
	responses_ZZ_p.SetLength(num_queries * words_per_block);
	for (dbsize_t i = 0; i < num_queries * words_per_block; ++i) {
	    responses_ZZ_p[i] = ZZ_p::zero();
	}
    }

    dbsize_t worker_response_size = subdb_words_per_block * bytes_per_word;
    unsigned char * worker_response = new unsigned char[worker_response_size];
    //unsigned char worker_response[worker_response_size];
    for (dbsize_t i = 0; i < vsplit; ++i) {
	for (dbsize_t j = 0; j < hsplit; ++j) {
	    std::iostream& wio = *(workerio[(hsplit * i) + j]);
	    for (nqueries_t k = 0; k < num_queries; ++k) {
		wio.read((char *)worker_response, worker_response_size);
		if (wio.eof() || (dbsize_t)(wio.gcount()) != worker_response_size) {
		    std::cerr << "Error: Worker " << (hsplit * i) + j
			    << "did not send a full response.\n";
		    delete[] responses;
		    delete[] worker_response;
#ifdef SPIR_SUPPORT
		    delete[] spir_queries;
#endif
		    return false;
		}
		// Add result to responses
		switch(mode) {
		case MODE_GF28:
		case MODE_GF216:
		case MODE_CHOR:
		    XOR_equal(responses + worker_response_size * (k * hsplit + j), 
			    worker_response, worker_response_size);
		    break;
		case MODE_ZZ_P:
		    for (dbsize_t w = 0; w < subdb_words_per_block; ++w) {
			ZZ_p addword = to_ZZ_p(ZZFromBytes(worker_response + (bytes_per_word * w), 
				bytes_per_word));
			if (IsZero(addword)) {
			    continue;
			}
			responses_ZZ_p[subdb_words_per_block * (k * hsplit + j) + w] += addword;
		    }
		    break;
		}
	    }
	}
    }
    delete[] worker_response;

    if (mode == MODE_ZZ_P) {
	for (nqueries_t k = 0; k < num_queries; ++k) {
	    vec_ZZ_p response;
	    response.SetLength(words_per_block);
	    for (dbsize_t w = 0; w < words_per_block; ++w) {
		response[w] = responses_ZZ_p[k * words_per_block + w];
	    }
#ifdef SPIR_SUPPORT
	    if (params.spir()) {
		spir_queries[k].randomize_response(response);
	    }
#endif
	    for (dbsize_t w = 0; w < words_per_block; ++w) {
		BytesFromZZ(responses + (k * words_per_block + w) * bytes_per_word,
			rep(response[w]), bytes_per_word);
	    }
	}
    }

    // Send the responses
    os.write((char*)responses, num_queries * words_per_block * bytes_per_word);
    os.flush();

    delete[] responses;

    return true;
}

