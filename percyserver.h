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

#ifndef __PERCYSERVER_H__
#define __PERCYSERVER_H__

#include <iostream>
#include <vec_ZZ_p.h>
#include "datastore.h"
#include "percyparams.h"
#include "gf2e.h"

NTL_CLIENT

// The class for Percy++ servers
class PercyServer {
public:
    // Initialize a server with the given DataStore.
    PercyServer(DataStore * datastore) : byzantine(false),
                     synchronized(true),
					 datastore(datastore) {}
    ~PercyServer() {}

    // Tell the server to be Byzantine
    void be_byzantine() { byzantine = true; }

    // Tell the server to tolerate unsynchronized databases (and hence compute the content hashes)
    template <typename GF2E_Element>
    void tolerate_unsynchronized(PercyServerParams &params);
    
    // Handle a request.
    virtual bool handle_request(PercyServerParams &params, std::istream &is,
	    std::ostream &os);

private:
    bool handle_request_ZZ_p(PercyServerParams &params, std::istream &is, 
	    std::ostream &os);
    bool handle_request_GF28(PercyServerParams &params, std::istream &is,
	    std::ostream &os);
    template <typename GF2E_Element>
    bool handle_request_GF2E(PercyServerParams &params, 
    	    std::istream &is, std::ostream &os);
    bool handle_request_Chor(PercyServerParams &params, std::istream &is,
	   std::ostream &os);
    // template <typename GF2E_Element>
    template <typename GF2E_Element>
    bool handle_request_RS_Sync(PercyServerParams &params, std::istream &is,
	   std::ostream &os);
   template <typename GF2E_Element>
    bool handle_hash_request_RS_Sync(PercyServerParams &params, std::istream &is, std::ostream &os, nqueries_t num_queries);
       
    // hash computation methods
    void compute_hashes(PercyServerParams &params);
    
    bool byzantine;
    bool synchronized;
    void compute_one(ZZ_p *value, bool hybrid_protection,
	dbsize_t num_blocks, nqueries_t num_queries,
	const vec_ZZ_p *inputvector, dbsize_t c);

protected:
    DataStore *datastore;
};

#endif
