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

#include <sys/time.h>
#include "mpicomm.h"
#include "mpiserver.h"

PercyMPIServerParams::PercyMPIServerParams(dbsize_t words_per_block, dbsize_t num_blocks,
	nservers_t tau, ZZ modulus, PercyMode mode, bool be_byzantine,
	char *pcparams_filename, bool do_spir, nservers_t sid,
	dbsize_t vsplit, dbsize_t hsplit, int rank, int num_nodes,
	bool test_client, nqueries_t num_test_queries) :
    PercyDistServerParams(words_per_block, num_blocks, tau, modulus, mode,
	    be_byzantine, pcparams_filename, do_spir, sid, vsplit, hsplit),
    _rank(rank),
    _num_nodes(num_nodes),
    test_client(test_client),
    _num_test_queries(num_test_queries)
{
    init();
}

PercyMPIServerParams::PercyMPIServerParams(PercyServerParams &sparams, dbsize_t vsplit,
	dbsize_t hsplit, int rank, int num_nodes, bool test_client,
	nqueries_t num_test_queries) :
    PercyDistServerParams(sparams, vsplit, hsplit),
    _rank(rank),
    _num_nodes(num_nodes),
    test_client(test_client),
    _num_test_queries(num_test_queries)
{
    init();
}

PercyMPIServerParams::PercyMPIServerParams(PercyDistServerParams &dsparams, int rank,
	int num_nodes, bool test_client, nqueries_t num_test_queries) :
    PercyDistServerParams(dsparams),
    _rank(rank),
    _num_nodes(num_nodes),
    test_client(test_client),
    _num_test_queries(num_test_queries)
{
    init();
}

PercyMPIServerParams::~PercyMPIServerParams ()
{
    for (dbsize_t i = 0; i < vsplit * hsplit; ++i) {
	delete workerio[i]->rdbuf();
	delete workerio[i];
    }
}

void PercyMPIServerParams::init ()
{
    if (_rank == 0) {
	_role = MPI_MASTER;
    } else if ((dbsize_t)_rank <= vsplit * hsplit) {
	_role = MPI_WORKER;
    } else if ((dbsize_t)_rank == vsplit * hsplit + 1) {
	_role = MPI_TEST_CLIENT;
    } else {
	_role = MPI_NONE;
    }
}

int PercyMPIServerParams::rank () const {
    return _rank;
}

int PercyMPIServerParams::num_nodes () const {
    return _num_nodes;
}

bool PercyMPIServerParams::has_test_client () const {
    return test_client;
}

nqueries_t PercyMPIServerParams::num_test_queries () const {
    return _num_test_queries;
}

MPIRole PercyMPIServerParams::role () const {
    return _role;
}


void print_usage_options ();

void print_usage_mpi (const char * bin) {
    // Only print if MASTER
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
	return;
    }

    std::cerr << "Usage: " << bin << " [OPTIONS...] VSPLIT HSPLIT SUBDB_1 SUBDB_2 ... SUBDB_M" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Start a PIR server as a master for the specified workers." << std::endl;
    std::cerr << std::endl;

    // Explanation of arguments for Master and MPI usage
    std::cerr << "Arguments:" << std::endl;
    std::cerr << "   VPSLIT       The database is split into VSPLIT sets of rows." << std::endl;
    std::cerr << "   HSPLIT       The database is split into HSPLIT sets of columns." << std::endl;
    std::cerr << "   SUBDB_i      The subdatabase file at index i where M=VSPLIT*HSPLIT." << std::endl;
    std::cerr << "                The order is important as index i represents the subdatabase at position" << std::endl;
    std::cerr << "                (i % VSPLIT, i / HSPLIT) in the matrix of subdatabases." << std::endl;
    std::cerr << std::endl;

    print_usage_options();
}

// Forward declaration
PercyServerParams * init_params(ParsedArgs &pargs, bool checkdb);

PercyMPIServerParams * init_params_mpi (ParsedArgs &pargs)
{
    // Get basic server params
    PercyServerParams * sparams = init_params(pargs, false);
    if (sparams == NULL) {
	return NULL;
    }
    
    // Get MPI information
    int rank, num_nodes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Check if there is a test client
    bool test_client = false;
    nqueries_t num_test_queries = 0;
    optmap::iterator oiter = pargs.opts.find('M');
    if (oiter != pargs.opts.end()) {
	test_client = true;
	num_test_queries = strtoll(oiter->second, NULL, 10);
    }

    // Get splitsA
    dbsize_t vsplit = strtoull(pargs.nonoptv[0], NULL, 10);
    dbsize_t hsplit = strtoull(pargs.nonoptv[1], NULL, 10);
    if (vsplit == 1 && hsplit == 1) {
	std::cerr << "Warning: Splitting the database 1x1.\n";
    }

    // Check that there are enough nodes
    if ((dbsize_t)num_nodes < vsplit * hsplit + ( test_client ? 2 : 1 )) {
	std::cerr << "There are not enough MPI nodes\n";
	print_usage_mpi (pargs.exec);
	return NULL;
    }

    // Check that there are the right number of non-option arguments
    int args_needed = vsplit * hsplit + 2;
    if (pargs.nonoptc < args_needed) {
	std::cerr << "Not enough arguments\n";
	print_usage_mpi (pargs.exec);
	return NULL;
    }

    // Check if the database is evenly split
    dbsize_t num_blocks = sparams->num_blocks();
    dbsize_t words_per_block = sparams->num_blocks();
    if (num_blocks % vsplit != 0) {
	fprintf(stderr, "VSPLIT does not evenly split the database.\n");
	return NULL;
    }
    if (words_per_block % hsplit != 0) {
	fprintf(stderr, "HSPLIT does not evenly split the database.\n");
	return NULL;
    }

    // Create MPI params
    PercyMPIServerParams * params = new PercyMPIServerParams(*sparams,
	    vsplit, hsplit, rank, num_nodes, test_client, num_test_queries);
    delete sparams;

    // Add worker iostreams if MASTER
    if (params->role() == MPI_MASTER) {
	for (dbsize_t i = 0; i < vsplit * hsplit; ++i) {
	    MPIStreamBuf * sb = new MPIStreamBuf(i + 1);
	    std::iostream * sio = new std::iostream(sb);
	    params->set_workerio(i, sio);
	}
    }

    return params;
}


// Forward declaration
dbsize_t database_bytes (const char * database);

DataStore * init_datastore_mpi (PercyMPIServerParams &params, ParsedArgs &pargs)
{
    const char * database = pargs.nonoptv[params.rank() + 1];
    fprintf(stderr, "[WORKER %d] Using database '%s'\n", params.rank(), database);

    dbsize_t dbsize = database_bytes(database);
    if (dbsize == 0) {
	fprintf(stderr, "Error: The database must exist and be non-empty\n");
	return NULL;
    }

    DataStore * ds = new FileDataStore(database, params, params.tau());
    return ds;
}


/*
int check_params (const PercyServerParams& params, bool silent = false) {
    if (params.get_disttype() != DIST_TYPE_NORMAL && 
	    params.split_horizontal() == 1 && params.split_vertical() == 1) {
	if (!silent) {
	    std::cerr << "No database split!\n";
	}
	return 5;
    }
    if (params.words_per_block() % params.split_horizontal() != 0) {
	if (!silent) {
	    std::cerr << "Horizontal split does not evenly partition the database\n";
	}
	return 6;
    }
    if (params.num_blocks() % params.split_vertical() != 0) {
	if (!silent) {
	    std::cerr << "Vertical split does not evenly partition the database\n";
	}
	return 7;
    }
    return 0;
}

int master (int argc, char ** argv, int size) {
    //attach_gdb();
    PercyServerParams * params = 0;
    DataStore * datastore = 0;
    init(argc, argv, &params, &datastore);
    if (datastore == 0) {
	std::cerr << "DataStore not initialized\n";
	return 1;
    }
    if (params == 0) {
	std::cerr << "PercyServerParams not initialized\n";
	return 2;
    }
    int ret = check_params(*params);
    if (ret != 0) {
	return ret;
    }
    return run (*params, datastore, MASTER_RANK, size, true);
}
*/

int mpi_test_client (PercyMPIServerParams& params) {
    nqueries_t num_tests = params.num_test_queries();

    MPIStreamBuf siobuf(0);
    std::iostream sio(&siobuf);

    std::cerr << "CLIENT: Sending params...\n";
    sio << params;

    unsigned char failure;
    sio.read((char*)&failure, 1);
    fprintf(stderr, "CLIENT: Got failure: 0x%x\n", failure);
    if (failure & 1 || failure & 2) {
	std::cerr << "CLIENT: Master did not accept parameters\n";
	return 1;
    }

    std::cerr << "Doing " << num_tests << " queries.\n";
    struct timeval ts, te;
    unsigned char nq[2];
    for (unsigned int ti = 0; ti < num_tests; ++ti) {
	gettimeofday(&ts, NULL);

	std::cerr << "CLIENT: Sending number of queries...\n";
	nqueries_t num_queries = 1;
	nq[0] = (num_queries >> 8) & 0xff;
	nq[1] = num_queries & 0xff;
	sio.write((char *)nq, 2);

	dbsize_t words_per_block = params.words_per_block();
	dbsize_t bytes_per_word = params.bytes_per_word();

	std::cerr << "CLIENT: Sending queries...\n";
	GF28_Element one = 0x01;
	GF28_Element zero = 0x00;
	(void)one;
	(void)zero;
	dbsize_t num_blocks = params.num_blocks();
	GF28_Element query[num_blocks];
	for (dbsize_t i = 0; i < num_blocks; ++i) {
	    ZZ randzz;
	    RandomBits(randzz, 8 * bytes_per_word);
	    BytesFromZZ(query + i, randzz, 1);
	}
	sio.write((char*)query, num_blocks);

	std::cerr << "CLIENT: Waiting for reply...\n";
	GF28_Element result[words_per_block];
	sio.read((char*)result, words_per_block);

	gettimeofday(&te, NULL);
	int td_sec = te.tv_sec - ts.tv_sec;
	int td_usec = te.tv_usec - ts.tv_usec;
	if (td_usec < 0) {
	    td_sec -= 1;
	    td_usec += 1000000;
	}

	fprintf(stderr, "TIME: %d.%06d seconds\n", td_sec, td_usec);

	std::cerr << "CLIENT: Got reply.\n";
	/*
	for (dbsize_t i = 0; i < words_per_block*bytes_per_word; ++i) {
	    fprintf(stdout, "%c", result[i]);
	}
	fflush(stdout);
	*/
    }

    nq[0] = 0;
    nq[1] = 0;
    sio.write((char *)nq, 2);

    return 0;
}

/*
int worker (int argc, char ** argv, int rank, int size) {
    PercyServerParams * params = 0;
    DataStore * datastore = 0;
    init(argc, argv, &params, &datastore);
    if (datastore == 0) {
	std::cerr << "DataStore not initialized\n";
	return 1;
    }
    if (params == 0) {
	std::cerr << "PercyServerParams not initialized\n";
	return 2;
    }
    // Stop worker if --single
    if (params->get_disttype() == DIST_TYPE_NORMAL) {
	return 0;
    }
    int ret = check_params(*params, true);
    if (ret != 0) {
	return ret;
    }

    // Fix params
    dbsize_t new_wpb = params->words_per_block() / params->split_horizontal();
    dbsize_t new_nb = params->num_blocks() / params->split_vertical();
    nservers_t new_tau = params->tau();
    ZZ new_modulus = params->get_modulus();
    PercyMode new_mode = params->get_mode();
    bool new_be_byzantine = params->is_byzantine();
    char * new_pcparams_filename = params->get_pcparams_filename();
    bool new_do_spir = params->spir();
    nservers_t new_sid = params->get_sid();
    uint16_t new_port = params->get_port();
    PercyServerParams newparams(new_wpb, new_nb, new_tau, new_modulus, new_mode,
	    0, 0, DIST_TYPE_WORKER, new_be_byzantine, new_pcparams_filename, 
	    new_do_spir, new_sid, new_port);

    return run (newparams, datastore, rank, size, false);
}
*/

/*
int main (int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < MASTER_RANK + 1) {
	std::cerr << "Not enough processes.  Quitting process " << rank << "...\n";
	return 1;
    }

    // Initialize ZZ_p
    ZZ modinit = to_ZZ("257");
    ZZ_p::init(modinit);
    // Initialize random stream
    unsigned char randbuf[128];
    //std::ifstream urand("/dev/urandom");
    std::ifstream urand("/home/cjdevet/getdata.py");
    urand.read((char *)randbuf, sizeof(randbuf));
    urand.close();
    ZZ randzz = ZZFromBytes(randbuf, sizeof(randbuf));
    SetSeed(randzz);

    // Get number of tests from args
    char * newargv[argc-1];
    newargv[0] = argv[0];
    for (unsigned int i = 1; i < (unsigned int)(argc-1); ++i) {
	newargv[i] = argv[i+1];
    }
    unsigned int num_tests = atoi(argv[1]);

    int ret = 0;
    if (rank == MASTER_RANK) { // Master Server
	ret = master(argc-1, newargv, size);
    }
#ifdef WITH_CLIENT
    else if (rank == CLIENT_RANK) { // Client
	ret = client(argc-1, newargv, size, num_tests);
    }
#endif
    else { // Worker
	ret = worker(argc-1, newargv, rank, size);
    }

    MPI_Finalize();
    return ret;
}
*/
