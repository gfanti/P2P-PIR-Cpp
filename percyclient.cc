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
#include <vec_vec_ZZ_p.h>
#include <ZZ_pX.h>
#include <algorithm>
#include "percyclient.h"
#include "rsdecoder.h"
#include "rsdecoder_impl.h"
#include "gf2e.h"

#ifdef SPIR_SUPPORT
#include "spirclient.h"
#endif

// Simple protected constructor for PercyClient
PercyClient::PercyClient (PercyClientParams &params, nservers_t num_servers, 
	nservers_t t) :
    params(params),
    num_servers(num_servers),
    t(t),
    goodservers(),
    requested_blocks(),
    received_blocks()
{}

PercyClient::PercyClient (const PercyClient &other) :
    params(other.params),
    num_servers(other.num_servers),
    t(other.t),
    goodservers(other.goodservers),
    requested_blocks(other.requested_blocks),
    received_blocks(other.received_blocks)
{}

PercyClient::~PercyClient () {}

// A factory method used to get a PercyClient for the given mode
PercyClient * PercyClient::make_client (PercyClientParams &params, 
	nservers_t num_servers, nservers_t t, sid_t * sids)
{
    PercyClient * retptr = NULL;
    std::cerr << "Mode is " << params.get_mode() << std::endl;
    switch (params.get_mode()) {
        case MODE_ZZ_P:
            retptr = new PercyClient_ZZ_p(params, num_servers, t, sids);
            break;
        case MODE_CHOR:
            retptr = new PercyClient_Chor(params, num_servers, t);
            break;
        case MODE_GF28:
            retptr = new PercyClient_GF2E<GF28_Element>(params, num_servers, t, sids);
            break;
        case MODE_GF216:
            retptr = new PercyClient_GF2E<GF216_Element>(params, num_servers, t, sids);
            break;
        case MODE_RS_SYNC:
            retptr = new PercyClient_RS_Sync<GF216_Element>(params, num_servers, t, sids);
            break;
        // case MODE_PULSE_SYNC:
        // retptr = new PercyClient_PULSE_Sync<GF216_Element>(params, num_servers, t, sids);
        // break;
        default:
            break;
    }
    return retptr;
}

// Do all steps in our method.
// current_results will be set to be the results of decoding the blocks
//     requested with the current call.
// previous_results will be set to be the results of decoding all blocks
//     that were not decoded previously.
nqueries_t PercyClient::fetch_blocks(vector<dbsize_t> block_numbers,
	vector<ostream*> &osvec, vector<istream*> &isvec,
	vector<PercyBlockResults> &current_results,
	vector<PercyBlockResults> &previous_results)
{
    std::cerr << "New call of fetch_blocks \n";
    int res = send_request(block_numbers, osvec);
    std::cerr << "made it back from fetch_blocks \n";
    if (res < 0) {
        return block_numbers.size();
    }
    nservers_t k = receive_replies(isvec);
    std::cerr << "Just received replies \n";
    // Calculate what h to use
    nservers_t h = (nservers_t)(floor(sqrt((t+params.tau())*k)))+1;
    const char *envh = getenv("PIRC_H");
    if (envh) {
        nservers_t override_h = atoi(envh);
        if (override_h > 0) {
            h = override_h;
        }
    }
    nqueries_t num_results = received_blocks.size();
    nqueries_t prev_queries = num_results - block_numbers.size();
    vector<PercyBlockResults> results;
    dbsize_t pr_res = process_replies(h, results);
    for (nqueries_t q = 0; q < prev_queries; ++q) {
        std::cerr << "looping/pushing back results\n";
        previous_results.push_back(results[q]);
    }
    for (nqueries_t q = prev_queries; q < num_results; ++q) {
        std::cerr << "looping/pushing back results, 2nd loop\n";
        previous_results.push_back(results[q]);
    }
    return pr_res;
}

nqueries_t PercyClient::fetch_blocks(vector<dbsize_t> block_numbers,
	vector<ostream*> &osvec, vector<istream*> &isvec,   	
	vector<PercyBlockResults> &current_results)
{
    vector<PercyBlockResults> dummy_previous;
    return fetch_blocks(block_numbers, osvec, isvec, current_results,
	    dummy_previous);
}


// Some helper methods

void XOR_equal(unsigned char *dst, const unsigned char *src, unsigned const int len) {
    for(unsigned int i=0; i<len; i++) {
        *(dst++) ^= *(src++);
    }
}

// Generate t-private (t+1)-of-l shares of a given secret value.
// If a non-NULL pointer *polyp is provided, the random polynomial used is 
// copied to that address.
static void genshares(nservers_t t, nservers_t l,
        const vec_ZZ_p &indices, vec_ZZ_p &values, dbsize_t secret, 
        ZZ_pX *polyp = NULL)
{
    // Pick a random polynomial of degree t
    ZZ_pX randpoly = random_ZZ_pX(t+1);
    // Set the constant term to the secret
    SetCoeff(randpoly, 0, secret);
    // std::cerr << "Poly is (" << randpoly << ")\n";

    // Evaluate the polynomial at each of the indices
    for (nservers_t i=0; i<l; ++i) {
        eval(values[i], randpoly, indices[i]);
        // std::cerr << "(" << indices[i] << ", " << values[i] << ")\n\n";
    }

    // Store this polynomial if requested by the caller
    if (polyp != NULL) {
        *polyp = randpoly;
    }
}


// ZZ_P CLASS

// Constructor and destructor
PercyClient_ZZ_p::PercyClient_ZZ_p (PercyClientParams &params, 
	nservers_t num_servers, nservers_t t, sid_t * sids) :
    PercyClient(params, num_servers, t),
    indices(),
    randmults(),
    answers(),
    unfinished_results(),
    decoded()
{
    choose_indices(sids);
}

PercyClient_ZZ_p::PercyClient_ZZ_p (const PercyClient_ZZ_p &other) :
    PercyClient(other),
    indices(other.indices),
    randmults(other.randmults),
    answers(other.answers),
    unfinished_results(other.unfinished_results),
    decoded(other.decoded)
{}

PercyClient_ZZ_p& PercyClient_ZZ_p::operator= (PercyClient_ZZ_p other)
{
    PercyClient::operator=(other);
    indices = other.indices;
    randmults = other.randmults;
    answers = other.answers;
    unfinished_results = other.unfinished_results;
    decoded = other.decoded;
    return *this;
}

PercyClient_ZZ_p::~PercyClient_ZZ_p () {}

// TODO
//void PercyClient::choose_indices_ZZ_p (vec_ZZ_p * indices, sid_t *sids) {
void PercyClient_ZZ_p::choose_indices (sid_t * sids) 
{
    indices.SetLength(num_servers);
    for (nservers_t j = 0; j < num_servers; ++j) {
        if (params.tau() > 0 || params.spir()) {
            // Use the constant indices 1, 2, ..., num_servers.
            indices[j] = to_ZZ_p(sids[j]);
        } else {
            // Use random indices
            ZZ_p r;
            bool ok = false;
            do {
                r = random_ZZ_p();
                if (r != 0) {
                    ok = true;
                    for (nservers_t k = 0; k < j; ++k) {
                        if (indices[k] == r) {
                            ok = false;
                        }   
                    }   
                }   
            } while (!ok);
            indices[j] = r;
        }
    }
}

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
int PercyClient_ZZ_p::send_request(vector<dbsize_t> block_numbers,
        std::vector<ostream*> &osvec)
{
    nqueries_t num_queries = block_numbers.size();
 
    if (num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << num_servers
            << ".\n";
        return -1;
    }

    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();
    params.mod_modulus();

    dbsize_t num_blocks = params.num_blocks();

#ifdef SPIR_SUPPORT
    SPIRClientQuery * spir_queries = new SPIRClientQuery[num_queries];
#endif

    // Construct the vector of server indices
    nqueries_t q;

    // Generate random multiples (!= 0)
    nqueries_t previous_queries = requested_blocks.size();
    if (randomize && !params.spir()) {
        if (params.spir()) {
            std::cerr << "Cannot randomize queries with SPIR\n";
    #ifdef SPIR_SUPPORT
            delete[] spir_queries;
    #endif
            return -1;
        }
        for (nqueries_t q = 0; q < num_queries; ++q) {
            randmults.push_back(vector<ZZ_p>());
            vector<ZZ_p>& query_randmults = randmults.back();
            for (nservers_t j = 0; j < num_servers; ++j) {
                ZZ_p r;
                while (IsZero(r)) {
                    r = random_ZZ_p();
                }
                query_randmults.push_back(r);
            }
        }
    }

    // Construct the shares of the e_{index} vector
    vec_vec_ZZ_p * shares = new vec_vec_ZZ_p[num_queries];
    vec_ZZ_pX * polyvec = new vec_ZZ_pX[num_queries];

    for (q=0; q<num_queries; ++q) {
	shares[q].SetLength(num_blocks);
	polyvec[q].SetLength(num_blocks);
	shares[q].SetLength(num_blocks);
	for (dbsize_t i = 0; i < num_blocks; ++i) {
	    params.mod_modulussq();
	    shares[q][i].SetLength(num_servers);
	    params.mod_modulus();
	    genshares(t, num_servers, indices, shares[q][i],
		    i == block_numbers[q], &(polyvec[q][i]));
	}
#ifdef SPIR_SUPPORT
        if (params.spir()) {
	    spir_queries[q].init_parameters(params, num_servers, t,
		    block_numbers[q], indices, shares[q], polyvec[q]);
        }
#endif
    }

    // Multiply shares by random multiples
    if (randomize && !params.spir()) {
	for (nservers_t j = 0; j < num_servers; ++j) {
	    for (q = 0; q < num_queries; ++q) {
                for (dbsize_t i = 0; i < num_blocks; ++i) {
                    shares[q][i][j] *= randmults[q+previous_queries][j];
                }
            }
        }
    }

    // Optionally encrypt the shares
    if (params.hybrid()) {
        for (dbsize_t i = 0; i < num_blocks; ++i) {
            for (nservers_t j = 0; j < num_servers; ++j) {
                for (q=0; q<num_queries; ++q) {
                    params.mod_modulus();
                    ZZ share = rep(shares[q][i][j]);
                    params.mod_modulussq();
                    shares[q][i][j] = params.encrypt(share);
                }
                std::cerr << i << " " << j << "\n";
            }
        }
        // If we're encrypting, leave modulus^2 as the active modulus
        params.mod_modulussq();
    }

    // Send the query to each server
    unsigned int modulus_bytes = params.hybrid() ? params.modulussq_bytes() :
        params.modulus_bytes();
    for (nservers_t j = 0; j < num_servers; ++j) {
        //	std::cerr << "Sending number of queries to server " << j << "...";
        unsigned char nq[2];
        nq[0] = (num_queries >> 8) & 0xff;
        nq[1] = (num_queries) & 0xff;
        osvec[j]->write((char *)nq, 2);
        //	cerr << "done" << std::endl;

        for (q=0; q<num_queries; ++q) {
#ifdef SPIR_SUPPORT
            if (params.spir()) {
		spir_queries[q].send_to_server(params, j, *(osvec[j]));
            }
#endif

            //		std::cerr << "Sending query " << q << " to server " << j << "...";
	    unsigned char * bytes = new unsigned char[modulus_bytes];
            for (dbsize_t i = 0; i < num_blocks; ++i) {
                BytesFromZZ(bytes, rep(shares[q][i][j]), modulus_bytes);
                osvec[j]->write((char *)bytes, modulus_bytes);
            }
	    delete[] bytes;
            //		std::cerr << "done" << std::endl;
        }
        osvec[j]->flush();
    }

    // Add block numbers to requests_blocks
    for (q = 0; q < num_queries; ++q) {
	requested_blocks.push_back(block_numbers[q]);
    }

#ifdef SPIR_SUPPORT
    delete[] spir_queries;
#endif
    delete[] shares;
    delete[] polyvec;

    savectx.restore();
    return 0;
}

// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
nservers_t PercyClient_ZZ_p::receive_replies(std::vector<istream*> &isvec)
{
    dbsize_t words_per_block = params.words_per_block();
    nqueries_t num_queries = requested_blocks.size();
    nqueries_t q;

    // Choose the right modulus
    unsigned int modulus_bytes;
    if (params.hybrid()) {
        modulus_bytes = params.modulussq_bytes();
        params.mod_modulussq();
    } else {
        modulus_bytes = params.modulus_bytes();
        params.mod_modulus();
    }

    // The vector of servers that have responded properly
    goodservers.clear();

    // Add to the appropriate vectors
    nqueries_t previous_queries = received_blocks.size();
    for (q=0; q<num_queries; ++q) {
	answers.push_back(vector<vec_ZZ_p>(words_per_block));
        for (dbsize_t c = 0; c < words_per_block; ++c) {
            answers[q+previous_queries][c].SetLength(num_servers);
        }
    }

    // Read the replies
    unsigned char * bytes = new unsigned char[modulus_bytes];
    for (nservers_t j = 0; j < num_servers; ++j) {
        bool isgood = true;
	for (q=0; q<num_queries; ++q) {
	    for (dbsize_t i = 0; isgood && i < words_per_block; ++i) {
                isvec[j]->read((char *)bytes, modulus_bytes);
                if (isvec[j]->eof()) {
                    std::cerr << "Server " << j+1 << " did not send complete reply.\n";
                    std::cerr << "Marking server " << j+1 << " as bad.\n";
                    isgood = false;
                    break;
                }
                if ((dbsize_t)(isvec[j]->gcount()) < modulus_bytes) {
                    // Mark this server as bad
                    std::cerr << "Marking server " << j+1 << " as bad.\n";
                    isgood = false;
                    break;
                }
                ZZ ans;
                ZZFromBytes(ans, bytes, modulus_bytes);
                answers[q+previous_queries][i][j] = to_ZZ_p(ans);
            }
        }
        if (isgood) {
            goodservers.push_back(j);
        }
    }
    delete[] bytes;

    // Add to unfinished_results and decoded
    for (q=0; q<num_queries; ++q) {
	unfinished_results.push_back(
		vector<DecoderResult<ZZ_p> >(1, DecoderResult<ZZ_p>(goodservers, map<dbsize_t, ZZ_p>())));
	decoded.push_back(std::set<dbsize_t>());
    }

    // Moved block numbers from requests_blocks to received_blocks
    for (q = 0; q < num_queries; ++q) {
	received_blocks.push_back(requested_blocks[q]);
    }
    requested_blocks.clear();

    // Optionally decrypt the answers
    if (params.hybrid()) {
        params.mod_modulussq();
        for (dbsize_t i = 0; i < words_per_block; ++i) {
            for (nservers_t j = 0; j < num_servers; ++j) {
                for (q=0; q<num_queries; ++q) {
                    ZZ_p dec = params.decrypt(answers[q+previous_queries][i][j]);
                    clear(answers[q+previous_queries][i][j]);
                    answers[q+previous_queries][i][j] = dec;
                }
            }
        }
    }

    // Remove random multiple
    if (randomize && !params.spir()) {
        for (nservers_t j = 0; j < num_servers; ++j) {
            //ZZ_p randmult_inv = inv(randmults[j]);
            for (q = 0; q < num_queries; ++q) {
                ZZ_p randmult_inv = inv(randmults[q][j]);
                for (dbsize_t i = 0; i < words_per_block; ++i) {
                    answers[q+previous_queries][i][j] *= randmult_inv;
                }
            }
        }
	randmults.clear();
    }

    // Now we're mod modulus for sure
    params.mod_modulus();

    return goodservers.size();
}

nqueries_t PercyClient_ZZ_p::process_replies(nservers_t h, vector<PercyBlockResults> &results)
{
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t bytes_per_word = params.bytes_per_word();
    nservers_t tau = params.tau();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return false;
    }

    std::cerr << goodservers.size() << " of " << num_servers << " servers responded.\n";

    RSDecoder_ZZ_p decoder(params.get_p1(), params.get_p2());

    // Call the decoder's Recover method
    bool res = decoder.Recover(bytes_per_word, t+tau, h, goodservers, answers, 
	    indices, unfinished_results, decoded);

    // Convert the results to PercyBlockResults
    // Remove queries that are completely decoded
    for (nqueries_t q = 0; q < received_blocks.size(); ++q) {
	results.push_back(PercyBlockResults());
	PercyBlockResults &curr_block_res = results.back();
	curr_block_res.block_number = received_blocks[q];

	if (decoded[q].size() == answers[q].size()) {
	    unsigned char * sigma = new unsigned char[words_per_block * bytes_per_word];
	    for (dbsize_t i = 0; i < unfinished_results[q].size(); ++i) {
		for (dbsize_t j = 0; j < words_per_block; ++j) {
		    BytesFromZZ(sigma+(j*bytes_per_word),
			    rep(unfinished_results[q][i].recovered[j]), 
			    bytes_per_word);
		}
		curr_block_res.results.push_back(PercyResult(unfinished_results[q][i].G,
			string((char *)sigma, words_per_block * bytes_per_word)));
	    }
	    delete[] sigma;

	    // Remove query
	    received_blocks.erase(received_blocks.begin() + q);
	    answers.erase(answers.begin() + q);
	    unfinished_results.erase(unfinished_results.begin() + q);
	    decoded.erase(decoded.begin() + q);
	    --q;
	}
    }

    (void)res;
    return received_blocks.size();
}


// CHOR CLASS

// Constructor and destructor
PercyClient_Chor::PercyClient_Chor (PercyClientParams &params, 
	nservers_t num_servers, nservers_t t) :
    PercyClient(params, num_servers, t),
    answers()
{}

PercyClient_Chor::PercyClient_Chor (const PercyClient_Chor &other) :
    PercyClient(other),
    answers()
{
    dbsize_t words_per_block = params.words_per_block();
    for (nqueries_t q = 0; q < other.answers.size(); ++q) {
	answers.push_back(new unsigned char[num_servers * words_per_block]);
	memcpy(answers[q], other.answers[q], num_servers * words_per_block);
    }
}

PercyClient_Chor& PercyClient_Chor::operator= (PercyClient_Chor other)
{
    PercyClient::operator=(other);
    std::swap(answers, other.answers);
    return *this;
}

PercyClient_Chor::~PercyClient_Chor ()
{
    while (answers.size() > 0) {
	if (answers.back() != NULL) {
	    delete[] answers.back();
	}
	answers.pop_back();
    }
}


//int PercyClient::send_request_Chor(vector<dbsize_t> block_numbers,
//        std::vector<iostream*> &iosvec, unsigned char * &shares)
int PercyClient_Chor::send_request(vector<dbsize_t> block_numbers,
        std::vector<ostream*> &osvec)
{
    std::cerr << "send_request chor\n";
    nqueries_t num_queries = block_numbers.size();

    if (num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request_Chor.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << num_servers
            << ".\n";
        return -1;
    }

    dbsize_t num_bytes = params.num_blocks() / 8;

    unsigned char * shares = new unsigned char [num_servers*num_queries*num_bytes];
    for (nqueries_t q = 0; q < num_queries; q++) {
        // Generate random bytestrings for each server but the last
        nservers_t last_server = num_servers - 1;
        for (nservers_t j = 0; j < last_server; j++) {
            for (unsigned int i = 0; i < num_bytes; i++) {
                shares[(j*num_queries + q)*num_bytes + i] = RandomBits_ulong(8);
            }
        }
        // Create a bytestring out of the block_numbers vector, storing it in the location of the 
        // final server's bytestring
        memset(shares + (last_server*num_queries+q)*num_bytes, '\0', num_bytes);
        unsigned int block_byte = (unsigned int) block_numbers[q] / 8;
        unsigned char block_pos = (unsigned char) block_numbers[q] % 8;
        unsigned char byte_val = 1 << (7 - block_pos);
        shares[(last_server*num_queries + q)*num_bytes + block_byte] = byte_val;
        // Compute the final bytestring (such that the XOR of all the 
	// server bytestrings is the bytestring 
        // representing the query vector)
        for (nservers_t i = 0; i < last_server; i++) {
            XOR_equal(
                    shares + (last_server*num_queries + q)*num_bytes,
                    shares + (i*num_queries + q)*num_bytes, num_bytes);
        }
    }

    // Send the query to each server
    for (nservers_t j = 0; j < num_servers; ++j) {
        unsigned char nq[2];
        nq[0] = (num_queries >> 8) & 0xff;
        nq[1] = num_queries & 0xff;
        osvec[j]->write((char *)nq, 2);

        osvec[j]->write((char *)(shares + j*num_queries*num_bytes), num_queries*num_bytes);
        osvec[j]->flush();
    }
    delete[] shares;

    // Add block numbers to requests_blocks
    for (nqueries_t q = 0; q < num_queries; ++q) {
	requested_blocks.push_back(block_numbers[q]);
    }

    return 0;
}

//nservers_t PercyClient::receive_replies_Chor(std::vector<iostream*> &iosvec,
//        unsigned char * &answers)
nservers_t PercyClient_Chor::receive_replies(std::vector<istream*> &isvec)
{
    dbsize_t words_per_block = params.words_per_block();
    nqueries_t num_queries = requested_blocks.size();

    // The vector of servers that have responded properly
    goodservers.clear();

    // Allocate space for answers
    nqueries_t previous_queries = answers.size();
    for (nqueries_t q = 0; q < num_queries; ++q) {
	answers.push_back(new unsigned char[num_servers * words_per_block]);
    }

    // Read the replies
    for (nservers_t j = 0; j < num_servers; ++j) {
        bool isgood = true;
        for (nqueries_t q=0; q<num_queries; ++q) {
            for (dbsize_t i = 0; isgood && i < words_per_block; ++i) {
                isvec[j]->read((char *)(answers[q+previous_queries] 
			+ j*words_per_block + i), 1);
                if ((dbsize_t)(isvec[j]->gcount()) < 1) {
                    // Mark this server as bad
                    std::cerr << "Marking server " << j+1 << " as bad.\n";
                    isgood = false;
                    break;
                }
            }
        }
        if (isgood) {
            goodservers.push_back(j);
        }
    }

    // Moved block numbers from requests_blocks to received_blocks
    for (nqueries_t q = 0; q < num_queries; ++q) {
        received_blocks.push_back(requested_blocks[q]);
    }
    requested_blocks.clear();

    std::cerr << "Got replies...\n";

//    for(nservers_t j = 0; j < num_servers; j++) {
//    std::cerr << "Client: Server " << j << ": ";
//    printBS_debug(answers+(j*num_queries+num_queries-1)*words_per_block, words_per_block);
//    std::cerr << std::endl;

    return goodservers.size();
}

//vector< vector<PercyResult> > PercyClient::process_replies_Chor(
//        nservers_t h, unsigned char * &shares, unsigned char * &answers,
//        vector<dbsize_t> block_numbers, vector<iostream*>& iosvec)
nqueries_t PercyClient_Chor::process_replies (nservers_t h,
	vector<PercyBlockResults> &results)
{
    dbsize_t words_per_block = params.words_per_block();
    nqueries_t num_queries = received_blocks.size();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return num_queries;
    }

    std::cerr << goodservers.size() << " of " << num_servers << " servers responded." << std::endl;

    if (goodservers.size() < num_servers) {
        std::cerr << "Not enough servers responded." << std::endl;
	return num_queries;
    }

    unsigned char *result = new unsigned char[words_per_block];
    for (nqueries_t q = 0; q < num_queries; q++) {
	results.push_back(PercyBlockResults());
	results.back().block_number = received_blocks[q];
	vector<PercyResult>& block_result = results.back().results;

        // Recover the query block by XORing all the replies together
        memset(result, '\0', words_per_block);
        for (nservers_t j = 0; j < num_servers; j++) {
            XOR_equal(result, answers[q] + j*words_per_block, words_per_block);
        }

        vector<nservers_t> empty;
        block_result.push_back(PercyResult(empty, string((char *)result, words_per_block)));

	delete[] answers[q];
    }
    delete[] result;
    received_blocks.clear();
    answers.clear();

    return 0;
}

