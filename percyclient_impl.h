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

#ifndef __PERCYCLIENT_IMPL_H__
#define __PERCYCLIENT_IMPL_H__

#include <iostream>
#include <algorithm>
#include <deque>
#include "gf2e.h"
#include "pulse.h"
#include "percyparams.h"
#include "percyresult.h"

template <typename GF2E_Element>
inline void setCoeffs(GF2X &GF2X_P);

template <>
inline void setCoeffs<GF28_Element>(GF2X &GF2X_P) {
    SetCoeff(GF2X_P, 8, 1);
    SetCoeff(GF2X_P, 4, 1);
    SetCoeff(GF2X_P, 3, 1);
    SetCoeff(GF2X_P, 1, 1);
    SetCoeff(GF2X_P, 0, 1);
}

template <>
inline void setCoeffs<GF216_Element>(GF2X &GF2X_P) {
    SetCoeff(GF2X_P, 16, 1);
    SetCoeff(GF2X_P, 14, 1);
    SetCoeff(GF2X_P, 12, 1);
    SetCoeff(GF2X_P, 7, 1);
    SetCoeff(GF2X_P, 6, 1);
    SetCoeff(GF2X_P, 4, 1);
    SetCoeff(GF2X_P, 2, 1);
    SetCoeff(GF2X_P, 1, 1);
    SetCoeff(GF2X_P, 0, 1);
}

// Constructor and destructor
template <typename GF2E_Element>
PercyClient_GF2E<GF2E_Element>::PercyClient_GF2E (PercyClientParams &params,
	nservers_t num_servers, nservers_t t, sid_t * sids) :
    PercyClient(params, num_servers, t),
    indices(),
    randmults(),
    answers(),
    unfinished_results(),
    decoded()
{
    // Initialize the GF2E modulus
    GF2X MODULUS_P;
    setCoeffs<GF2E_Element>(MODULUS_P);

    GF2E::init(MODULUS_P);
    GF2X::HexOutput = 1;

    choose_indices(sids);
}

template <typename GF2E_Element>
PercyClient_GF2E<GF2E_Element>::PercyClient_GF2E (
	const PercyClient_GF2E<GF2E_Element> &other) :
    PercyClient(other),
    indices(new GF2E_Element[num_servers]),
    indices_ntl(other.indices_ntl),
    answers(vector<GF2E_Element *>(other.answers.size())),
    answers_ntl(other.answers_ntl),
    unfinished_results(other.unfinished_results),
    decoded(other.decoded)
{
    dbsize_t words_per_block = params.words_per_block();
    memcpy(indices, other.indices, num_servers * sizeof(GF2E_Element));
    for (nqueries_t q = 0; q < answers.size(); ++q) {
	answers[q] = new GF2E_Element[words_per_block * num_servers];
	memcpy(answers[q], other.answers[q], words_per_block * num_servers *
		sizeof(GF2E_Element));
    }
}

template <typename GF2E_Element>
PercyClient_GF2E<GF2E_Element>& PercyClient_GF2E<GF2E_Element>::operator= (
	PercyClient_GF2E<GF2E_Element> other)
{
    PercyClient::operator=(other);
    std::swap(indices, other.indices);
    indices_ntl = other.indices_ntl;
    randmults = other.randmults;
    std::swap(answers, other.answers);
    answers_ntl = other.answers_ntl;
    unfinished_results = other.unfinished_results;
    decoded = other.decoded;
    return *this;
}

template <typename GF2E_Element>
PercyClient_GF2E<GF2E_Element>::~PercyClient_GF2E ()
{
    if (indices != NULL) {
        delete[] indices;
    }
    while (answers.size() > 0) {
        if (answers.back() != NULL) {
            delete[] answers.back();
        }
        answers.pop_back();
    }
}


// Generate t-private (t+1)-of-l shares of a given secret value.
template <typename GF2E_Element>
static void genshares_GF2E(nservers_t t, nservers_t l,
    const GF2E_Element *indices, GF2E_Element *values, GF2E_Element secret) {
    // Pick a random polynomial of degree t with the right constant term
    GF2E_Element * coeffs = new GF2E_Element[t+1];
    coeffs[0] = secret;
    for (nservers_t i=1;i<=t;++i) {
        coeffs[i] = RandomBits_ulong(8*sizeof(GF2E_Element));
    }

    // Evaluate the polynomial at each of the indices
    for (nservers_t i=0; i<l; ++i) {
        values[i] = evalpoly_GF2E<GF2E_Element>(coeffs, t, indices[i]);
    }
    delete[] coeffs;
}

// Swap the query polynomials of unsynchronized files (this actually just zeroes out the appropriate entries)
template <typename GF2E_Element>
static void swap_symbols(nservers_t t, nservers_t num_servers, GF2E_Element *values, 
    dbsize_t num_blocks, const std::vector<dbsize_t> &sync_error_locs) {
    
    for (std::vector<dbsize_t>::const_iterator it = sync_error_locs.begin(); it != sync_error_locs.end(); ++it) {
        bool zeroed = true;
        dbsize_t sync_idx = *it;
        for (nservers_t j = 0; j < num_servers; j++) {
            // check if the appropriate values of values is nonzero
            if (values[(sync_idx * num_servers) + j] != 0) {
                zeroed = false;
                break;
            }
        }
        // if the unsynchronized file will be touched by a server, swap its bits with zeros elsewhere
        if (!zeroed) {
            for (nservers_t j = 0; j < num_servers; j++) {
                // // find zeros elsewhere in the vector
                // dbsize_t swap_idx = sync_idx;
                
                //For now, we are just zeroing out the offending entries instead of swapping them!
                values[(sync_idx * num_servers) + j] = 0;
                
                // // how do you check for equality ?
                // while (swap_idx == sync_idx) {
                    // // pick a random new index in the range [0,num_blocks-1]
                    // dbsize_t rand_idx = RandomBits_ulong(8*sizeof(dbsize_t));
                    // // make sure the randomly drawn index is in range and not mis-synchronized
                    // if ((rand_idx >= num_blocks) ||
                        // (std::find(sync_error_locs.begin(), sync_error_locs.end(), rand_idx) == sync_error_locs.end())){
                        // continue;
                    // }
                    
                    // // otherwise, swap the appropriate entries
                    // swap_idx = rand_idx;
                    // // std::swap(values[(num_blocks + sync_idx) * num_servers + j], values[(num_blocks + swap_idx) * num_servers + j]);
                    
                    
                // }
            }
        }
    }
    
    // Sanity check
    // for (std::vector<dbsize_t>::const_iterator it = sync_error_locs.begin(); it != sync_error_locs.end(); ++it) {
        // dbsize_t sync_idx = *it;
        // for (nservers_t j = 0; j < num_servers; j++) {
            // std::cerr << "Printing mis-synched random values: " << values[(num_blocks + sync_idx) * num_servers + j] << std::endl;
        // }
    // }
    
}

/* GF:This method seems to choose which indices X we want to evaluate the polynomials at. 
   I.e. if the query is a random polynomial q(X), we need to evaluate it at some random X. */
template <typename GF2E_Element>
void PercyClient_GF2E<GF2E_Element>::choose_indices(sid_t *sids) {
    indices = new GF2E_Element[num_servers];
    for (nservers_t j=0; j<num_servers; ++j) {
        if (params.tau()) {
            // Use the indices provided
            indices[j] = (GF2E_Element)sids[j];
        } else {
            // Use random indices
            GF2E_Element r;
            bool ok = false;
            do {
                r = RandomLen_long(8*sizeof(GF2E_Element));
                if (r != 0) {
                    ok = true;
                    for (nservers_t k=0;k<j;++k) {
                        if (indices[k] == r) {
                            ok = false;
                        }
                    }
                }
            } while (!ok);
            indices[j] = r;
        }
    }

    // Create NTL version
    indices_ntl.SetLength(num_servers);
    for (nservers_t ix=0; ix<num_servers; ++ix) {
	conv(indices_ntl[ix],
		GF2XFromBytes((unsigned char *)(indices + ix), 1));
    }

}

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
template <typename GF2E_Element>
int PercyClient_GF2E<GF2E_Element>::send_request(vector<dbsize_t> block_numbers,
        std::vector<ostream*> &osvec)
{
    nqueries_t num_queries = block_numbers.size();

    if (num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request_GF2E.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << num_servers
            << ".\n";
        return -1;
    }

    // Construct the vector of server indices
    nqueries_t q;

    // Generate random multiples (!= 0)
    nqueries_t previous_queries = requested_blocks.size();
    if (randomize) {
        for (nqueries_t q = 0; q < num_queries; ++q) {
            randmults.push_back(vector<GF2E_Element>());
            vector<GF2E_Element>& query_randmults = randmults.back();
            for (nservers_t j = 0; j < num_servers; ++j) {
                GF2E_Element r = 0;
                while (r == 0) {
                    r = RandomLen_long(8*sizeof(GF2E_Element));
                }
                query_randmults.push_back(r);
            }
        }
    }

    // Construct the shares of the e_{index} vector
    dbsize_t num_blocks = params.num_blocks();
    GF2E_Element * shares = new GF2E_Element[num_queries * num_blocks * num_servers];

    for (q=0; q<num_queries; ++q) {
        for (dbsize_t i = 0; i < num_blocks; ++i) {
            genshares_GF2E<GF2E_Element>(t, num_servers, indices,
                    shares + (q * num_blocks + i) * num_servers, i == block_numbers[q]);
        }
    }

    // Multiply shares by random multiples
    if (randomize) {
        for (nservers_t p = 0; p < num_servers; ++p) {
            for (nqueries_t i = 0; i < num_queries; ++i) {
                for (dbsize_t j = 0; j < num_blocks; ++j) {
                    dbsize_t index = (i * num_blocks + j) * num_servers + p;
                    shares[index] = 
                                multiply_GF2E<GF2E_Element>(shares[index],
                                        randmults[i+previous_queries][p]);
                }
            }
        }
    }

    // Send the params and query to each server
    for (nservers_t j = 0; j < num_servers; ++j) {
        unsigned char nq[2];
        nq[0] = (q >> 8) & 0xff;
        nq[1] = (q) & 0xff;
        osvec[j]->write((char *)nq, 2);
        for (q=0; q<num_queries; ++q) {
            for (dbsize_t i = 0; i < num_blocks; ++i) {
                dbsize_t index = (q * num_blocks + i) * num_servers + j;
                char *shareptr = (char *)&(shares[index]);
                osvec[j]->write(shareptr, sizeof(GF2E_Element));
            }
        }
        osvec[j]->flush();
    }

    // Add block numbers to requests_blocks
    for (q = 0; q < num_queries; ++q) {
	requested_blocks.push_back(block_numbers[q]);
    }

    delete[] shares;

    return 0;
}


// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
template <typename GF2E_Element>
nservers_t PercyClient_GF2E<GF2E_Element>::receive_replies (
	std::vector<istream*> &isvec)
{
    dbsize_t words_per_block = params.words_per_block();
    nqueries_t num_queries = requested_blocks.size();
    nqueries_t q;

    // The vector of servers that have responded properly
    goodservers.clear();

    // The responses from the servers
    nqueries_t previous_queries = received_blocks.size();
    for (q = 0; q < num_queries; ++q) {
	answers.push_back(new GF2E_Element[words_per_block * num_servers]);
    }
    nqueries_t total_queries = answers.size();

    // Read the replies
    for (nservers_t j = 0; j < num_servers; ++j) {
        bool isgood = true;
        for (q=previous_queries; q<total_queries; ++q) {
            for (dbsize_t i = 0; isgood && i < words_per_block; ++i) {
                isvec[j]->read((char *)(answers[q] + i * num_servers + j),
			sizeof(GF2E_Element));
                if (isvec[j]->eof()) {
                    std::cerr << "Server " << j+1 << " did not send complete reply.\n";
                    std::cerr << "Marking server " << j+1 << " as bad.\n";
                    isgood = false;
                    break;
                }
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

    // Add to unfinished_results and decoded
    for (q=0; q<num_queries; ++q) {
	unfinished_results.push_back(
		vector<DecoderResult<GF2E> >(1, DecoderResult<GF2E>(goodservers, map<dbsize_t, GF2E>())));
	decoded.push_back(std::set<dbsize_t>());
    }

    // Moved block numbers from requests_blocks to received_blocks
    for (q = 0; q < num_queries; ++q) {
	received_blocks.push_back(requested_blocks[q]);
    }
    requested_blocks.clear();

    // Remove random multiple
    if (randomize) {
        for (nservers_t j = 0; j < num_servers; ++j) {
            for (nqueries_t q = 0; q < num_queries; ++q) {
                GF2E_Element randmult_inv =
                        inverse_GF2E<GF2E_Element>(randmults[q][j]);
                for (dbsize_t i = 0; i < words_per_block; ++i) {
                    answers[q+previous_queries][i*num_servers + j] = multiply_GF2E<GF2E_Element>(
			    answers[q+previous_queries][i*num_servers + j], randmult_inv);
                }
            }
        }
	randmults.clear();
    }

    std::cerr << "Got replies...\n";

    return goodservers.size();
}

// Construct the t+1 Lagrange coefficents to interpolate at the point
// alpha from the source points
// indices[goodservers[firstpoint .. firstpoint+numpoints-1]]
template <typename GF2E_Element>
void PercyClient_GF2E<GF2E_Element>::construct_lagrange_coeffs(GF2E_Element *coeffs,
	GF2E_Element alpha, nservers_t firstpoint, nservers_t numpoints)
{
    for (nservers_t i=0; i < numpoints; ++i) {
        GF2E_Element numer = 1, denom = 1;
        for (nservers_t j=0; j < numpoints; ++j) {
            if (j==i) continue;
            GF2E_Element numerdiff = indices[goodservers[firstpoint+j]] ^
	    		             alpha;
            GF2E_Element denomdiff = indices[goodservers[firstpoint+j]] ^
				     indices[goodservers[firstpoint+i]];
            numer = multiply_GF2E<GF2E_Element>(numer, numerdiff);
            denom = multiply_GF2E<GF2E_Element>(denom, denomdiff);
        }
        coeffs[i] = multiply_GF2E<GF2E_Element>(numer,
			    inverse_GF2E<GF2E_Element>(denom));
    }
}

// Do Lagrange interpolation of source values
// answers[goodservers[firstpoint..firstpoint+numpoints-1]] with the
// coefficients coeffs[0..t]
template <typename GF2E_Element>
inline GF2E_Element PercyClient_GF2E<GF2E_Element>::interpolate(
	const GF2E_Element *word_answers, const GF2E_Element *coeffs, 
	nservers_t firstpoint, nservers_t numpoints)
{
    GF2E_Element res = 0;

    for (nservers_t i = 0; i < numpoints; ++i) {
	res ^= multiply_GF2E<GF2E_Element>(coeffs[i],
		word_answers[goodservers[firstpoint+i]]);
    }

    return res;
}

template <typename GF2E_Element>
bool PercyClient_GF2E<GF2E_Element>::try_fast_recover(nservers_t h,
	vector<PercyBlockResults> &results)
{
    dbsize_t words_per_block = params.words_per_block();
    nqueries_t num_queries = answers.size();
    nservers_t k = goodservers.size();

    if (k <= t) return false;

    // Construct the Lagrange coefficients for the target index 0, as
    // well as for target indices[goodservers[0..(k-t-2)]] from the
    // source indices[goodservers[k-t-1 .. k-1]]

    // lagrange_coeffs[0] is the array of t+1 coeffs to interpolate 0
    // largrnge_coeffs[a+1] is the array of t+1 coeffs to interpolate
    //                                        indices[goodservers[a]]
    GF2E_Element * lagrange_coeffs = new GF2E_Element[(k-t) * (t+1)];
    construct_lagrange_coeffs(lagrange_coeffs, 0, k-t-1, t+1);
    for (nservers_t a = 0; a < k-t-1; ++a) {
	construct_lagrange_coeffs(lagrange_coeffs + (a+1) * (t+1),
	    indices[goodservers[a]], k-t-1, t+1);
    }

    results.clear();
    GF2E_Element *resblock = new GF2E_Element[words_per_block];
    for (nqueries_t q = 0; q < num_queries; ++q) {
	results.push_back(PercyBlockResults());
	results.back().block_number = received_blocks[q];
	for (dbsize_t word = 0; word < words_per_block; ++word) {
	    // Check if all the servers agree
	    bool match = true;
	    for (nservers_t a = 0; a < k-t-1; ++a) {
		GF2E_Element expectedans =
			answers[q][word*num_servers + goodservers[a]];
		GF2E_Element actualans = interpolate(
			answers[q] + word*num_servers, lagrange_coeffs + (a+1) * (t+1), 
			k-t-1, t+1);

		if (expectedans != actualans) {
		    match = false;
		    break;
		}
	    }
	    if (!match) {
		delete[] lagrange_coeffs;
		delete[] resblock;
		results.clear();
		return false;
	    }
	    resblock[word] = interpolate(answers[q] + word*num_servers, 
		    lagrange_coeffs, k-t-1, t+1);
	}
	PercyResult thisresult(goodservers, string((char *)resblock, 
		words_per_block*sizeof(GF2E_Element)));
	results.back().results.push_back(thisresult);
    }

    delete[] lagrange_coeffs;
    delete[] resblock;
    return true;
}

//Process the received replies and return the decoded results as a
//PercyResults object.
template <typename GF2E_Element>
nqueries_t PercyClient_GF2E<GF2E_Element>::process_replies (
	nservers_t h, vector<PercyBlockResults> &results)
{
    dbsize_t words_per_block = params.words_per_block();
    nservers_t tau = params.tau();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return false;
    }

    std::cerr << goodservers.size() << " of " << num_servers << " servers responded.\n";

    if (try_fast_recover(h, results)) {
	// Clean up
	received_blocks.clear();
	while (answers.size() > 0) {
	    if (answers.back() != NULL) {
		delete[] answers.back();
	    }
	    answers.pop_back();
	}
	answers_ntl.clear();
	unfinished_results.clear();
	decoded.clear();
	return 0;
    }

    // Create NTL answers not yet created
    for (nqueries_t q = answers_ntl.size(); q < received_blocks.size(); ++q) {
	answers_ntl.push_back(vector<vec_GF2E>(words_per_block));
	for (dbsize_t c = 0; c < words_per_block; ++c) {
	    answers_ntl[q][c].SetLength(num_servers);
	    for (nservers_t j = 0; j < num_servers; ++j) {
		conv(answers_ntl[q][c][j],
			GF2XFromBytes((unsigned char *)(answers[q] +
			    c * num_servers + j), sizeof(GF2E_Element)));
	    }
	}
    }

    RSDecoder_GF2E decoder;

    // Call the decoder's Recover method
    bool res = decoder.Recover(sizeof(GF2E_Element), t+tau, h, goodservers,
	    answers_ntl, indices_ntl, unfinished_results, decoded);

    // Convert the results to PercyBlockResults
    // Remove queries that are completely decoded
    for (nqueries_t q = 0; q < received_blocks.size(); ++q) {
	results.push_back(PercyBlockResults());
	results.back().block_number = received_blocks[q];
	vector<PercyResult> &block_results = results.back().results;

	if (decoded[q].size() == words_per_block) {
	    GF2E_Element *sigma = new GF2E_Element[words_per_block];
	    for (dbsize_t i = 0; i < unfinished_results[q].size(); ++i) {
		for (dbsize_t j = 0; j < words_per_block; ++j) {
		    BytesFromGF2X((unsigned char *)(sigma + j),
			    rep(unfinished_results[q][i].recovered[j]), 
			    sizeof(GF2E_Element));
		}
		block_results.push_back(PercyResult(unfinished_results[q][i].G, 
			string((char *)sigma, words_per_block * sizeof(GF2E_Element))));
	    }
	    delete[] sigma;

	    // Remove query
	    received_blocks.erase(received_blocks.begin() + q);
	    if (answers[q] != NULL) {
		delete[] answers[q];
	    }
	    answers.erase(answers.begin() + q);
	    answers_ntl.erase(answers_ntl.begin() + q);
	    unfinished_results.erase(unfinished_results.begin() + q);
	    decoded.erase(decoded.begin() + q);
	    --q;
	}
    }

    (void)res;
    return received_blocks.size();
}


// Constructor and destructor
template <typename GF2E_Element>
PercyClient_RS_Sync<GF2E_Element>::PercyClient_RS_Sync (PercyClientParams &params,
	nservers_t num_servers, nservers_t t, sid_t * sids) :
    PercyClient_GF2E<GF2E_Element>::PercyClient_GF2E(params, num_servers, t, sids),
    sync_error_locs()
{ }

template <typename GF2E_Element>
PercyClient_RS_Sync<GF2E_Element>::PercyClient_RS_Sync (
	const PercyClient_RS_Sync<GF2E_Element> &other) :
    PercyClient_GF2E<GF2E_Element>::PercyClient_GF2E (&other),
    sync_error_locs(other.sync_error_locs)
{}

template <typename GF2E_Element>
PercyClient_RS_Sync<GF2E_Element>& PercyClient_RS_Sync<GF2E_Element>::operator= (
	PercyClient_RS_Sync<GF2E_Element> other)
{
    PercyClient_GF2E<GF2E_Element>::operator=(other);
    sync_error_locs = other.sync_error_locs;
    return *this;
}

template <typename GF2E_Element>
PercyClient_RS_Sync<GF2E_Element>::~PercyClient_RS_Sync ()
{
    PercyClient_GF2E<GF2E_Element>::~PercyClient_GF2E ();
}

template <typename GF2E_Element>
void PercyClient_RS_Sync<GF2E_Element>::choose_indices(sid_t *sids) {
    
    PercyClient_GF2E<GF2E_Element>::choose_indices(sids);
}


// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
template <typename GF2E_Element>
int PercyClient_RS_Sync<GF2E_Element>::send_request(vector<dbsize_t> block_numbers,
        std::vector<ostream*> &osvec)
{
    nqueries_t num_queries = block_numbers.size();

    if (this->num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request_GF2E.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << this->num_servers
            << ".\n";
        return -1;
    }

    // Construct the vector of server indices
    nqueries_t q;

    // Generate random multiples (!= 0)
    nqueries_t previous_queries = this->requested_blocks.size();
    if (this->randomize) {
        for (nqueries_t q = 0; q < num_queries; ++q) {
            this->randmults.push_back(vector<GF2E_Element>());
            vector<GF2E_Element>& query_randmults = this->randmults.back();
            for (nservers_t j = 0; j < this->num_servers; ++j) {
                GF2E_Element r = 0;
                while (r == 0) {
                    r = RandomLen_long(8*sizeof(GF2E_Element));
                }
                query_randmults.push_back(r);
            }
        }
    }
    
    // Construct the shares of the e_{index} vector
    dbsize_t num_blocks = this->params.num_blocks();
    GF2E_Element * shares = new GF2E_Element[num_queries * num_blocks * this->num_servers];

    for (q=0; q<num_queries; ++q) {
        for (dbsize_t i = 0; i < num_blocks; ++i) {
            genshares_GF2E<GF2E_Element>(this->t, this->num_servers, this->indices,
                    shares + (q * num_blocks + i) * this->num_servers, i == block_numbers[q]);
        }
    }
    
    swap_symbols<GF2E_Element>(this->t, this->num_servers, shares, num_blocks, this->sync_error_locs);

    // Multiply shares by random multiples
    if (this->randomize) {
        for (nservers_t p = 0; p < this->num_servers; ++p) {
            for (nqueries_t i = 0; i < num_queries; ++i) {
                for (dbsize_t j = 0; j < num_blocks; ++j) {
                    dbsize_t index = (i * num_blocks + j) * this->num_servers + p;
                    shares[index] = 
                        multiply_GF2E<GF2E_Element>(shares[index],
                                this->randmults[i+previous_queries][p]);
                }
            }
        }
    }
    
    // Send the params and query to each server
    for (nservers_t j = 0; j < this->num_servers; ++j) {
        unsigned char nq[2];  //this stores the number of queries, and is read by the server
        nq[0] = (q >> 8) & 0xff;
        nq[1] = (q) & 0xff;
        osvec[j]->write((char *)nq, 2);    
        for (q=0; q<num_queries; ++q) {
            for (dbsize_t i = 0; i < num_blocks; ++i) {
                dbsize_t index = (q * num_blocks + i) * this->num_servers + j;
                char *shareptr = (char *)&(shares[index]);
                osvec[j]->write(shareptr, sizeof(GF2E_Element));
            }
        }
        osvec[j]->flush();
    }
    
    // Add block numbers to requests_blocks
    for (q = 0; q < num_queries; ++q) {
        this->requested_blocks.push_back(block_numbers[q]);
    }
    
    delete[] shares;

    return 0;
}

// Send a request for the synchronization information to the
// servers connected with the ostreams in the given vector.
template <typename GF2E_Element>
int PercyClient_RS_Sync<GF2E_Element>::send_sync_request(std::vector<ostream*> &osvec)
{
    if (this->num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request_GF2E.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << this->num_servers
            << ".\n";
        return -1;
    }
    
    // Build the query polynomial ( q(X) = X^t + X^(t-1) + ... + X + 1 )
    GF2E_Element * coeffs = new GF2E_Element[this->t+1];
    coeffs[0] = 0;
    for (nservers_t i=1;i<=this->t;++i) {
        coeffs[i] = 1;
    }
    
    // Send the query to each server
    for (nservers_t j = 0; j < this->num_servers; ++j) {
        
        GF2E_Element q_x = evalpoly_GF2E<GF2E_Element>(coeffs, this->t, j+1);
        unsigned char query[2];  //this stores the value of the query polynomial
        query[0] = (q_x >> 8) & 0xff;
        query[1] = (q_x & 0xff);
        osvec[j]->write((char *)query, 2); 
        osvec[j]->flush();
    }

    delete[] coeffs;
    return 0;
}


// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
template <typename GF2E_Element>
nservers_t PercyClient_RS_Sync<GF2E_Element>::receive_replies (
	std::vector<istream*> &isvec)
{
    dbsize_t words_per_block = this->params.words_per_block();
    nqueries_t num_queries = this->requested_blocks.size();
    nqueries_t q;

    // The vector of servers that have responded properly
    this->goodservers.clear();

    // The responses from the servers
    nqueries_t previous_queries = this->received_blocks.size();
    for (q = 0; q < num_queries; ++q) {
        this->answers.push_back(new GF2E_Element[words_per_block * this->num_servers]);
    }
    nqueries_t total_queries = this->answers.size();
    
    // Read the replies
    for (nservers_t j = 0; j < this->num_servers; ++j) {
        bool isgood = true;
        for (q=previous_queries; q<total_queries; ++q) {
            for (dbsize_t i = 0; isgood && i < words_per_block; ++i) {
                isvec[j]->read((char *)(this->answers[q] + i * this->num_servers + j), sizeof(GF2E_Element));
                if (isvec[j]->eof()) {
                    std::cerr << "Server " << j+1 << " did not send complete reply.\n";
                    std::cerr << "Marking server " << j+1 << " as bad (first).\n";
                    isgood = false;
                    break;
                }
                if ((dbsize_t)(isvec[j]->gcount()) < 1) {
                    // Mark this server as bad
                    std::cerr << "Marking server " << j+1 << " as bad (second).\n";
                    isgood = false;
                    break;
                }
            }
        }
        if (isgood) {
            this->goodservers.push_back(j);
        }
    }
    
    // Add to unfinished_results and decoded
    for (q=0; q<num_queries; ++q) {
	this->unfinished_results.push_back(
		vector<DecoderResult<GF2E> >(1, DecoderResult<GF2E>(this->goodservers, map<dbsize_t, GF2E>())));
	this->decoded.push_back(std::set<dbsize_t>());
    }

    // Moved block numbers from requests_blocks to received_blocks
    for (q = 0; q < num_queries; ++q) {
	this->received_blocks.push_back(this->requested_blocks[q]);
    }
    this->requested_blocks.clear();

    // Remove random multiple
    if (this->randomize) {
        for (nservers_t j = 0; j < this->num_servers; ++j) {
            for (nqueries_t q = 0; q < num_queries; ++q) {
                GF2E_Element randmult_inv =
                        inverse_GF2E<GF2E_Element>(this->randmults[q][j]);
                for (dbsize_t i = 0; i < words_per_block; ++i) {
                    this->answers[q+previous_queries][i*this->num_servers + j] = multiply_GF2E<GF2E_Element>(
			    this->answers[q+previous_queries][i*this->num_servers + j], randmult_inv);
                }
            }
        }
	this->randmults.clear();
    }

    std::cerr << "Got replies...\n";

    return this->goodservers.size();
}

// Receive the server's replies to the sync request, return # of servers that
// gave complete (but not necessarily correct) replies.
template <typename GF2E_Element>
nservers_t PercyClient_RS_Sync<GF2E_Element>::receive_sync_replies (
	std::vector<istream*> &isvec)
{
    // The vector of servers that have responded properly
    nservers_t res = 0;
    dbsize_t words_per_block = this->params.words_per_block();
    dbsize_t max_unsynchronized = this->params.max_unsynchronized();
    float expansion_factor = this->params.expansion_factor();
    
    dbsize_t num_rows = (dbsize_t) max_unsynchronized * expansion_factor * NUM_RATIOS;
    
    
    // For each query, read the input vector, which is a sequence of
    // num_blocks entries, each of length sizeof(GF2E_Element) bytes
    // GF2E_Element *input = new GF2E_Element[words_per_block];
    // std::vector<std::array<GF2E_Element,WORDS_PER_BLOCK> > input(num_rows);
    std::vector<GF2E_Element*> replies;
    for (dbsize_t i=0; i<this->num_servers; i++) {
        replies.push_back(new GF2E_Element[num_rows*WORDS_PER_BLOCK]);
    }
    
    // Read in the output from the various servers!
    for (nservers_t j = 0; j < this->num_servers; ++j) {
        isvec[j]->read((char *) replies[j], num_rows * WORDS_PER_BLOCK * sizeof(GF2E_Element));
    }
    
    // Interpolate the results
    std::vector<GF2E_Element*> compressed_results;
    for (dbsize_t i=0; i<num_rows; i++) {
        compressed_results.push_back(new GF2E_Element[WORDS_PER_BLOCK]);
    }
    interpolate_results(replies, this->num_servers, num_rows, compressed_results);
    
    
    
    find_unsynchronized_files(this->num_servers, num_rows, compressed_results);
    
    // Free the input memory
    while (replies.size() > 0) {
        if (replies.back() != NULL) {
            delete[] replies.back();
        }
        replies.pop_back();
    }
    
    // Free the input memory
    while (compressed_results.size() > 0) {
        if (compressed_results.back() != NULL) {
            delete[] compressed_results.back();
        }
        compressed_results.pop_back();
    }
    
    return res;
}

// Receive the server's replies to the sync request, return # of servers that
// gave complete (but not necessarily correct) replies.
template <typename GF2E_Element>
void PercyClient_RS_Sync<GF2E_Element>::interpolate_results (
	std::vector<GF2E_Element*> &replies, nservers_t num_servers, 
    dbsize_t num_rows, std::vector<GF2E_Element*> &compressed_results)
{
    dbsize_t result_len = num_rows * WORDS_PER_BLOCK;
    GF2E_Element *outputvec = new GF2E_Element[result_len];
    memset(outputvec, '\0', result_len*sizeof(GF2E_Element));
    const GF216_Element *block;
    
    for (dbsize_t j = 0; j < num_servers; ++j) {
        GF216_Element inpv_j; // fill this in! should be [0 0 ... 1]*V^-1
        switch(num_servers) {
            case 3:
                inpv_j = GF216_V_inv_3servers[j];
                break;
            case 4:
                inpv_j = GF216_V_inv_4servers[j];
                break;
            case 5:
                inpv_j = GF216_V_inv_5servers[j];
                break;
            default:
                inpv_j = 1;
                std::cout << "You need to work out the bottom row of V inverse for " << num_servers << " servers.\n";
        } 
        if (inpv_j != 0) {
            block = replies[j];
            const GF216_Element *blockc = block;
            GF216_Element log_j = GF216_log_table[inpv_j];
            const GF216_Element *start = GF216_exp_table + log_j;
            GF216_Element *oc = outputvec;
            GF216_Element *oc_end = oc + (result_len & ~3);
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
            for (dbsize_t c = 0; c < (result_len & 3); ++c, ++oc) {
                block_c = *(blockc++);
                if (block_c != 0) {
                    GF216_Element log_c = GF216_log_table[block_c];
                    *oc ^= start[log_c];
                }
            }
        }
    }
    
    
    for (dbsize_t i=0; i<num_rows; i++) {
        std::copy(&outputvec[WORDS_PER_BLOCK*i],&outputvec[WORDS_PER_BLOCK*(i+1)-1],compressed_results[i]);
    }
    delete[] outputvec;
}

template <typename GF2E_Element>
int PercyClient_RS_Sync<GF2E_Element>::checkSingletonRatio(
    const std::vector<GF2E_Element*> &compressed_results, dbsize_t i) {
    
    int singleton_idx = -1;
    bool singleton = true;
    GF2E_Element ratio1, ratio2;
    GF2E_Element previous_ratio = 1;
    for (dbsize_t j=0; j<4; j++) {
        ratio1 = multiply_GF2E<GF216_Element>(compressed_results[i+1][j],inverse_GF2E<GF216_Element>(compressed_results[i][j]));
        ratio2 = multiply_GF2E<GF216_Element>(compressed_results[i+2][j],inverse_GF2E<GF216_Element>(compressed_results[i+1][j]));
        if (i==2) {
            // std::cerr << "The ratios are " << ratio1 << " and " << ratio2 <<std::endl;
            // std::cerr << "The ratio1 stuff is " << compressed_results[i+1][j] << " and " << inverse_GF2E<GF216_Element>(compressed_results[i][j]) <<std::endl;
        }
        if (ratio1 != ratio2) {
            singleton = false;
            break;
        }
        if ((j>0) && (previous_ratio != ratio1)) {
            singleton = false;
            break;
        }
        previous_ratio = ratio1;
    }
    
    
    if (singleton && (previous_ratio != 0)) {
        singleton_idx = GF216_log_table[previous_ratio];
    }
    return singleton_idx;
    
}

template <typename GF2E_Element>
void PercyClient_RS_Sync<GF2E_Element>::find_unsynchronized_files(
    nservers_t num_servers, dbsize_t num_rows, std::vector<GF2E_Element*> &compressed_results){
    
    // do some PULSE decoding
    bool done = false;
    const GF216_Element *block;
    // this->sync_error_locs.push_back(1);
    
    while (!done) {
        done = true;
        for (dbsize_t i=0; i<num_rows; i+=NUM_RATIOS) {
            // find the location of the singleton, if there is one
            int singleton_idx = checkSingletonRatio(compressed_results,i);
            if (singleton_idx > -1) {
                done = false;
                
                // update the list of nonzero db entries
                std::cerr << "The singleton is located at index " << singleton_idx << std::endl;
                this->sync_error_locs.push_back(singleton_idx);
                
                // construct the list of indices to remove from the compressed vector
                deque<GF216_Element> bins;
                for (int z=0; z<DEGREE; z++) {
                    GF216_Element bin;
                    switch (num_rows/NUM_RATIOS) {
                        case 4:
                            bin = GF216_pulse_mtx_4bins[singleton_idx][z];
                            break;
                        case 6: 
                            bin = GF216_pulse_mtx_6bins[singleton_idx][z];
                            break;
                        case 8: 
                            bin = GF216_pulse_mtx_8bins[singleton_idx][z];
                            break;
                        case 10: 
                            bin = GF216_pulse_mtx_10bins[singleton_idx][z];
                            break;
                        case 12: 
                            bin = GF216_pulse_mtx_12bins[singleton_idx][z];
                            break;
                        default:
                            std::cerr << "Cannot compute: This number of bins has not been considered yet!";
                            break;
                    }
                    if (bin == i/3) {
                        bins.push_back(bin);
                    }
                    else {
                        bins.push_front(bin);
                    }
                }
                
                // remove the singleton from the results vector
                for (GF216_Element j=0; j<DEGREE; j++){
                    GF216_Element bin = bins.front();
                    bins.pop_front();
                    for (GF216_Element k=0; k<NUM_RATIOS; k++) {
                        // subtract the singleton from compressed results
                        block = compressed_results[i + k];
                        const GF216_Element *blockc = block;
                        // GF216_Element log_j = GF216_log_table[inpv_j];
                        const GF216_Element *start = GF216_exp_table;
                        GF216_Element *oc = compressed_results[bin*NUM_RATIOS + k];
                        GF216_Element *oc_end = oc + (WORDS_PER_BLOCK & ~3);
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
                        for (dbsize_t c = 0; c < (WORDS_PER_BLOCK & 3); ++c, ++oc) {
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
    }
}

template <typename GF2E_Element>
nqueries_t PercyClient_RS_Sync<GF2E_Element>::process_replies (
	nservers_t h, vector<PercyBlockResults> &results)
{
    dbsize_t words_per_block = this->params.words_per_block();
    nservers_t tau = this->params.tau();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return false;
    }

    std::cerr << this->goodservers.size() << " of " << this->num_servers << " servers responded.\n";

    if (this->try_fast_recover(h, results)) {
	// Clean up
	this->received_blocks.clear();
	while (this->answers.size() > 0) {
	    if (this->answers.back() != NULL) {
		delete[] this->answers.back();
	    }
	    this->answers.pop_back();
	}
	this->answers_ntl.clear();
	this->unfinished_results.clear();
	this->decoded.clear();
	return 0;
    }

    // Create NTL answers not yet created
    for (nqueries_t q = this->answers_ntl.size(); q < this->received_blocks.size(); ++q) {
	this->answers_ntl.push_back(vector<vec_GF2E>(words_per_block));
	for (dbsize_t c = 0; c < words_per_block; ++c) {
	    this->answers_ntl[q][c].SetLength(this->num_servers);
	    for (nservers_t j = 0; j < this->num_servers; ++j) {
		conv(this->answers_ntl[q][c][j],
			GF2XFromBytes((unsigned char *)(this->answers[q] +
			    c * this->num_servers + j), sizeof(GF2E_Element)));
	    }
	}
    }

    RSDecoder_GF2E decoder;

    // Call the decoder's Recover method
    bool res = decoder.Recover(sizeof(GF2E_Element), this->t+tau, h, this->goodservers,
	    this->answers_ntl, this->indices_ntl, this->unfinished_results, this->decoded);

    // Convert the results to PercyBlockResults
    // Remove queries that are completely decoded
    for (nqueries_t q = 0; q < this->received_blocks.size(); ++q) {
	results.push_back(PercyBlockResults());
	results.back().block_number = this->received_blocks[q];
	vector<PercyResult> &block_results = results.back().results;

	if (this->decoded[q].size() == words_per_block) {
	    GF2E_Element *sigma = new GF2E_Element[words_per_block];
	    for (dbsize_t i = 0; i < this->unfinished_results[q].size(); ++i) {
		for (dbsize_t j = 0; j < words_per_block; ++j) {
		    BytesFromGF2X((unsigned char *)(sigma + j),
			    rep(this->unfinished_results[q][i].recovered[j]), 
			    sizeof(GF2E_Element));
		}
		block_results.push_back(PercyResult(this->unfinished_results[q][i].G, 
			string((char *)sigma, words_per_block * sizeof(GF2E_Element))));
	    }
	    delete[] sigma;

	    // Remove query
	    this->received_blocks.erase(this->received_blocks.begin() + q);
	    if (this->answers[q] != NULL) {
		delete[] this->answers[q];
	    }
	    this->answers.erase(this->answers.begin() + q);
	    this->answers_ntl.erase(this->answers_ntl.begin() + q);
	    this->unfinished_results.erase(this->unfinished_results.begin() + q);
	    this->decoded.erase(this->decoded.begin() + q);
	    --q;
	}
    }

    (void)res;
    return this->received_blocks.size();
}

#endif
