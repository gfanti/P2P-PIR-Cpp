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

#ifndef __MPISERVER_H__
#define __MPISERVER_H__

#include "percyparams.h"
#include "cmdtools.h"
#include "distserver.h"

enum MPIRole {MPI_MASTER, MPI_WORKER, MPI_TEST_CLIENT, MPI_NONE};

class PercyMPIServerParams : public PercyDistServerParams {
public:
    PercyMPIServerParams(dbsize_t words_per_block, dbsize_t num_blocks,
	    nservers_t tau, ZZ modulus, PercyMode mode, bool be_byzantine,
	    char *pcparams_filename, bool do_spir, nservers_t sid,
	    dbsize_t vsplit = 1, dbsize_t hsplit = 1, int rank = 0,
	    int num_nodes = 1, bool test_client = false, 
	    nqueries_t num_test_queries = 0);

    PercyMPIServerParams(PercyServerParams &sparams, dbsize_t vsplit = 1,
	    dbsize_t hsplit = 1, int rank = 0, int num_nodes = 1,
	    bool test_client = false, nqueries_t num_test_queries = 0);

    PercyMPIServerParams(PercyDistServerParams &dsparams, int rank = 0,
	    int num_nodes = 1, bool test_client = false, 
	    nqueries_t num_test_queries = 0);

    ~PercyMPIServerParams ();

    int rank () const;
    int num_nodes () const;
    bool has_test_client () const;
    nqueries_t num_test_queries () const;
    MPIRole role () const;

protected:
    int _rank;
    int _num_nodes;
    bool test_client;
    nqueries_t _num_test_queries;
    MPIRole _role;

private:
    void init ();
};

void print_usage_mpi (const char * bin);

PercyMPIServerParams * init_params_mpi (ParsedArgs &pargs);

DataStore * init_datastore_mpi (PercyMPIServerParams &params, ParsedArgs &pargs);

int mpi_test_client (PercyMPIServerParams& params);

#endif
