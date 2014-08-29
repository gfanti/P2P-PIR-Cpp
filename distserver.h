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

#ifndef __DISTSERVER_H__
#define __DISTSERVER_H__

#include <vector>
#include "datastore.h"
#include "percyparams.h"
#include "percyserver.h"
#include "percytypes.h"

class PercyDistServerParams : public PercyServerParams {
public:
    PercyDistServerParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    dbsize_t expansion_factor, nservers_t tau, ZZ modulus, PercyMode mode, bool be_byzantine,
	    char *pcparams_filename, bool do_spir, nservers_t sid,
	    dbsize_t vsplit = 1, dbsize_t hsplit = 1) :
	PercyServerParams(words_per_block, num_blocks, max_unsynchronized, expansion_factor, tau, modulus, mode,
		be_byzantine, pcparams_filename, do_spir, sid),
	vsplit(vsplit),
	hsplit(hsplit),
	workerio(std::vector<std::iostream*>(vsplit * hsplit))
    {}

    PercyDistServerParams(PercyServerParams& sparams,
	    dbsize_t vsplit = 1, dbsize_t hsplit = 1) :
	PercyServerParams(sparams),
	vsplit(vsplit),
	hsplit(hsplit),
	workerio(std::vector<std::iostream*>(vsplit * hsplit))
    {}

    dbsize_t get_vsplit () const {
	return vsplit;
    }

    dbsize_t get_hsplit () const {
	return hsplit;
    }

    bool set_workerio (dbsize_t workerid, std::iostream * wio) {
	if (workerid >= vsplit * hsplit) {
	    return false;
	}
	workerio[workerid] = wio;
	return true;
    }

    bool unset_workerio (dbsize_t workerid) {
	if (workerid >= vsplit * hsplit) {
	    return false;
	}
	workerio[workerid] = NULL;
	return true;
    }

    bool have_all_workerio () {
	for (dbsize_t i = 0; i < vsplit * hsplit; ++i) {
	    if (workerio[i] == NULL) {
		return false;
	    }
	}
	return true;
    }

    std::vector<std::iostream*>& get_workerio () {
	return workerio;
    }

    PercyServerParams sub_params () const {
	dbsize_t new_wpb = _words_per_block / hsplit;
	dbsize_t new_nb = _num_blocks / vsplit;
	// NOTE: We assume all workers will not do SPIR
	return PercyServerParams(new_wpb, new_nb, _tau, _max_unsynchronized, _expansion_factor, modulus, mode, be_byzantine, 
		pcparams_filename, false, sid);
    }

protected:
    dbsize_t vsplit, hsplit;
    std::vector<std::iostream*> workerio;
};

class PercyMasterServer : public PercyServer {
public:
    // Initialize a server with the given datastore
    PercyMasterServer() :
	PercyServer(0)
    {}
    virtual ~PercyMasterServer() {}

    virtual bool handle_request(PercyServerParams &params, std::istream &is, std::ostream &os);
};

#endif
