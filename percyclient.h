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

#ifndef __PERCYCLIENT_H__
#define __PERCYCLIENT_H__

#include <vector>
#include <set>
#include <iostream>
#include <string.h>
#include <vec_vec_ZZ_p.h>
#include <vec_GF2E.h>
#include "percyresult.h"
#include "percyparams.h"
#include "gf2e.h"
#include "rsdecoder.h"

NTL_CLIENT


class PercyClient {
public:
    virtual ~PercyClient ();

    // A factory method used to get a PercyClient for the given mode
    static PercyClient * make_client (PercyClientParams &params, 
	    nservers_t num_servers, nservers_t t, sid_t * sids = NULL);

    // Send a request for the given block number (0-based) to the
    // servers connected with the ostreams in the given vector.
    virtual int send_request(vector<dbsize_t> block_numbers, 
	    std::vector<ostream*> &osvec) = 0;

    // Receive the servers' replies for all requests.  Return k, the number of
    // servers that returned anything.
    virtual nservers_t receive_replies(std::vector<istream*> &isvec) = 0;

    // Process the server's replies and pass the results to results.  The
    // number of queries that did not successfully decode is returned.
    virtual nqueries_t process_replies(nservers_t h, 
	    vector<PercyBlockResults> &results) = 0;

    // Do all of the above in one shot.
    // current_results will be set to be the results of decoding the blocks
    //     requested with the current call.
    // previous_results will be set to be the results of decoding all blocks
    //     that were not decoded previously.
    nqueries_t fetch_blocks(vector<dbsize_t> block_numbers,
	    vector<ostream*> &osvec, vector<istream*> &isvec,
	    vector<PercyBlockResults> &current_results,
	    vector<PercyBlockResults> &previous_results);
    nqueries_t fetch_blocks(vector<dbsize_t> block_numbers,
	    vector<ostream*> &osvec, vector<istream*> &isvec,
	    vector<PercyBlockResults> &current_results);

protected:
    // The constructor is protected so that only make_client can be used to
    // create one of the derived client types
    PercyClient (PercyClientParams &params, nservers_t num_servers,
	    nservers_t t);

    PercyClient (const PercyClient &other);

    // Members needed for all modes
    PercyClientParams params;
    nservers_t num_servers, t;
    vector<nservers_t> goodservers;
    static const bool randomize = true;
    vector<dbsize_t> requested_blocks;
    vector<dbsize_t> received_blocks;
};


class PercyClient_ZZ_p : public PercyClient {
public:
    PercyClient_ZZ_p (PercyClientParams &params, nservers_t num_servers,
	    nservers_t t, sid_t * sids);

    PercyClient_ZZ_p (const PercyClient_ZZ_p &other);
    PercyClient_ZZ_p& operator= (PercyClient_ZZ_p other);

    virtual ~PercyClient_ZZ_p ();

    // Virtual members as described in PercyClient class
    virtual int send_request(vector<dbsize_t> block_numbers, 
	    std::vector<ostream*> &osvec);
    virtual nservers_t receive_replies(std::vector<istream*> &isvec);
    virtual nqueries_t process_replies(nservers_t h,
	    vector<PercyBlockResults> &results);

private:
    // Based on the sids of the servers, choose the index that will
    // correspond with each server.  This function is called by the
    // constructor.
    virtual void choose_indices(sid_t *sids);

    // Private ZZ_p members
    vec_ZZ_p indices;
    vector<vector<ZZ_p> > randmults;
    vector<vector<vec_ZZ_p> > answers;
    vector<vector<DecoderResult<ZZ_p> > > unfinished_results;
    vector<std::set<dbsize_t> > decoded;
};

class PercyClient_Chor : public PercyClient {
public:
    PercyClient_Chor (PercyClientParams &params, nservers_t num_servers,
	    nservers_t t);

    PercyClient_Chor (const PercyClient_Chor &other);
    PercyClient_Chor& operator= (PercyClient_Chor other);

    virtual ~PercyClient_Chor ();

    // Virtual members as described in PercyClient class
    virtual int send_request(vector<dbsize_t> block_numbers, 
	    std::vector<ostream*> &osvec);
    virtual nservers_t receive_replies(std::vector<istream*> &isvec);
    virtual nqueries_t process_replies(nservers_t h,
	    vector<PercyBlockResults> &results);

private:
    // Private Chor members
    vector<unsigned char *> answers;
};

template<typename GF2E_Element>
class PercyClient_GF2E : public PercyClient {
public:
    PercyClient_GF2E (PercyClientParams &params, nservers_t num_servers,
	    nservers_t t, sid_t * sids);

    PercyClient_GF2E (const PercyClient_GF2E &other);
    PercyClient_GF2E& operator= (PercyClient_GF2E other);

    virtual ~PercyClient_GF2E ();

    // Virtual members as described in PercyClient class
    virtual int send_request(vector<dbsize_t> block_numbers, 
	    std::vector<ostream*> &osvec);
    virtual nservers_t receive_replies(std::vector<istream*> &isvec);
    virtual nqueries_t process_replies(nservers_t h,
	    vector<PercyBlockResults> &results);

private:
    virtual void choose_indices(sid_t *sids);

    // A NTL-less method to attempt a fast recovery
    bool try_fast_recover (nservers_t h, vector<PercyBlockResults> &results);
    // Some helpers for it
    void construct_lagrange_coeffs(GF2E_Element *coeffs, GF2E_Element alpha,
	    nservers_t firstpoint, nservers_t numpoints);
    inline GF2E_Element interpolate(const GF2E_Element *word_answers, 
	    const GF2E_Element *coeffs, nservers_t firstpoint,
	    nservers_t numpoints);

    // Private GF2E members
    GF2E_Element * indices;
    vec_GF2E indices_ntl;
    vector<vector<GF2E_Element> > randmults;
    vector<GF2E_Element *> answers;
    vector<vector<vec_GF2E> > answers_ntl;
    vector<vector<DecoderResult<GF2E> > > unfinished_results;
    vector<std::set<dbsize_t> > decoded;
};

#include "percyclient_impl.h"

#endif
