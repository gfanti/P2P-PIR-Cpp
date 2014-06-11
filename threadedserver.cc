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
#include <sstream>
#include <fstream>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <NTL/vec_ZZ_p.h>
#include "percyio.h"
#include "threadedserver.h"
#include "xor.h"

#ifdef SPIR_SUPPORT
#include "spirserver.h"
#endif

// Close fds 3 and up, except for the one given (pass -1 to close them
// all)
static void close_highfds_except(int exceptfd)
{
    // Find the max fd number
    struct rlimit limit;
    getrlimit(RLIMIT_NOFILE, &limit);
    for (int fd = 3; fd < (int)(limit.rlim_cur); ++fd) {
        if (fd != exceptfd) {
            // There's no ill effect from closing a non-open fd, so just
            // do it
            close(fd);
        }
    }
}

// Returns the number of bytes communicated for each word
dbsize_t comm_bytes_per_word (const PercyServerParams &params)
{
    switch (params.get_mode()) {
    case MODE_ZZ_P:
	// For ZZ_p we need to use the size of the modulus
	if (params.hybrid()) {
	    params.mod_modulussq();
	    return params.modulussq_bytes();
	} else {
	    params.mod_modulus();
	    return params.modulus_bytes();
	}
	break;
    case MODE_CHOR:
	// Always use 1 byte words for Chor
	return 1;
	break;
    case MODE_GF28:
    case MODE_GF216:
    default:
	// Otherwise, just use the given value
	return params.bytes_per_word();
	break;
    }
    return 0;
}

// Returns the number of blocks that a worker will process.
dbsize_t worker_num_blocks (const PercyServerParams &params, dbsize_t workerid)
{
    dbsize_t num_blocks = params.num_blocks();
    PercyMode mode = params.get_mode();
    if (mode == MODE_CHOR) {
	// For Chor, we want to divide the blocks into sets of 8
	num_blocks /= 8;
    }
    dbsize_t worker_nb = 0;
    switch (params.get_ttype()) {
    case THREADING_ROWS: {
	dbsize_t num_blocks_div = num_blocks / params.get_num_threads();
	dbsize_t num_blocks_mod = num_blocks % params.get_num_threads();
	worker_nb = num_blocks_div + ( workerid < num_blocks_mod ? 1 : 0 );
	} break;
    case THREADING_QUERIES: {
	worker_nb = num_blocks;
	} break;
    }
    if (mode == MODE_CHOR) {
	worker_nb *= 8;
    }
    return worker_nb;
}

// Returns the size (in bytes) of each query for the given workerid
dbsize_t worker_query_bytes (const PercyServerParams &params, dbsize_t workerid)
{
    dbsize_t worker_nb = worker_num_blocks(params, workerid);
    if (params.get_mode() == MODE_CHOR) {
	worker_nb /= 8;
    }
    return worker_nb * comm_bytes_per_word(params);
}

// Returns the offset (in number of blocks) from the start of a query for the
// given workerid
dbsize_t worker_query_offset (const PercyServerParams &params, dbsize_t workerid)
{
    dbsize_t num_blocks = params.num_blocks();
    PercyMode mode = params.get_mode();
    if (mode == MODE_CHOR) {
	num_blocks /= 8;
    }
    dbsize_t worker_off = 0;
    switch (params.get_ttype()) {
    case THREADING_ROWS: {
	dbsize_t num_blocks_div = num_blocks / params.get_num_threads();
	dbsize_t num_blocks_mod = num_blocks % params.get_num_threads();
	if (workerid < num_blocks_mod) {
	    worker_off = (num_blocks_div + 1) * workerid;
	} else {
	    worker_off = num_blocks_div * workerid + num_blocks_mod;
	}
	} break;
    case THREADING_QUERIES: {
	worker_off = 0;
	} break;
    }
    if (mode == MODE_CHOR) {
	worker_off *= 8;
    }
    return worker_off;
}

// Returns the offset (in bytes) from the start of a query for each query for
// the given workerid
dbsize_t worker_bytes_offset (const PercyServerParams &params, dbsize_t workerid)
{
    dbsize_t worker_off = worker_query_offset(params, workerid);
    if (params.get_mode() == MODE_CHOR) {
	worker_off /= 8;
    }
    return worker_off * comm_bytes_per_word(params);
}

// Returns the number of words per block for a worker (Here for if we add a
// THREADING_COLUMNS type)
dbsize_t worker_words_per_block (const PercyServerParams &params, dbsize_t workerid)
{
    switch (params.get_ttype()) {
    case THREADING_ROWS:
    case THREADING_QUERIES: {
	return params.words_per_block();
	} break;
    }
    return 0;
}

// Returns the number of queries that a worker will process.
nqueries_t worker_num_queries (const PercyServerParams &params, dbsize_t workerid, 
	nqueries_t num_queries)
{
    switch (params.get_ttype()) {
    case THREADING_ROWS: {
	return num_queries;
	} break;
    case THREADING_QUERIES: {
	nqueries_t num_queries_div = num_queries / params.get_num_threads();
	nqueries_t num_queries_mod = num_queries % params.get_num_threads();
	return num_queries_div + ( workerid < num_queries_mod ? 1 : 0 );
	} break;
    }
    return 0;
}

// Returns the index of the first query that a worker will index
nqueries_t worker_first_query (const PercyServerParams &params, dbsize_t workerid, 
	nqueries_t num_queries)
{
    switch (params.get_ttype()) {
    case THREADING_ROWS: {
	return 0;
	} break;
    case THREADING_QUERIES: {
	nqueries_t num_queries_div = num_queries / params.get_num_threads();
	nqueries_t num_queries_mod = num_queries % params.get_num_threads();
	if ((nqueries_t)workerid < num_queries_mod) {
	    return (num_queries_div + 1) * workerid;
	}
	return num_queries_div * workerid + num_queries_mod;
	} break;
    }
    return 0;
}

// Returns the params for the given workerid
PercyServerParams worker_params (const PercyServerParams &params, dbsize_t workerid)
{
    return PercyServerParams(worker_words_per_block(params, workerid),
	    worker_num_blocks(params, workerid), params.tau(), params.get_modulus(), 
	    params.get_mode(), params.is_byzantine(), NULL, false, params.get_sid());
}


ThreadedDataStore::ThreadedDataStore (const char * filename, const PercyServerParams &params,
	bool tau_independent) :
    FileDataStore(filename, params, tau_independent),
    num_threads(num_threads),
    ttype(ttype)
{
    num_threads = params.get_num_threads();
    ttype = params.get_ttype();
    switch (ttype) {
    case THREADING_ROWS:
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    unsigned char * subdbptr = database + 
		    worker_query_offset(params, i) * words_per_block
		    * bytes_per_word;
	    PercyServerParams subparams = worker_params(params, i);
	    subdatastores.push_back(new MemoryDataStore(subdbptr, subparams));
	}
	break;
    case THREADING_QUERIES:
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    PercyServerParams subparams = params;
	    subdatastores.push_back(new MemoryDataStore(database, subparams));
	}
	break;
    }
}

ThreadedDataStore::~ThreadedDataStore ()
{
    for (dbsize_t i = 0; i < num_threads; ++i) {
	delete subdatastores[i];
    }
}

MemoryDataStore * ThreadedDataStore::get_subds (dbsize_t index)
{
    if (index >= num_threads) {
	std::cerr << "Invalid thread index: " << index << "\n";
	return NULL;
    }
    return subdatastores[index];
}

dbsize_t ThreadedDataStore::get_num_threads ()
{
    return num_threads;
}


ThreadStreamBuf::ThreadStreamBuf () :
    inbuffers(),
    outbuffers(),
    inbufferindex(0),
    outbufferindex(0),
    incharindex(0),
    outcharindex(0)
{}

ThreadStreamBuf::~ThreadStreamBuf () {}

void ThreadStreamBuf::add_inbuffer (char * addr, dbsize_t size)
{
    BufferInfo bi;
    bi.addr = addr;
    bi.size = size;
    inbuffers.push_back(bi);
}

void ThreadStreamBuf::add_outbuffer (char * addr, dbsize_t size)
{
    BufferInfo bi;
    bi.addr = addr;
    bi.size = size;
    outbuffers.push_back(bi);
}

std::streamsize ThreadStreamBuf::showmanyc ()
{
    std::streamsize total = 0;
    std::vector<BufferInfo>::iterator iter;
    for (iter = inbuffers.begin(); iter != inbuffers.end(); ++iter) {
	total += ( iter->size < 0 ? 0 : iter->size );
    }
    return total;
}

int ThreadStreamBuf::underflow ()
{
    if (in_eof()) {
	return EOF;
    }
    char * currbuf = inbuffers[inbufferindex].addr;
    return traits_type::to_int_type(currbuf[incharindex]);
}

int ThreadStreamBuf::uflow ()
{
    int retval = underflow();
    if (retval == EOF) {
	return EOF;
    }
    ++incharindex;
    if (incharindex >= inbuffers[inbufferindex].size) {
	++inbufferindex;
	incharindex = 0;
    }
    return retval;
}

int ThreadStreamBuf::overflow (int c)
{
    if (out_eof()) {
	return EOF;
    }
    char * currbuf = outbuffers[outbufferindex].addr;
    currbuf[outcharindex] = traits_type::to_char_type(c);
    ++outcharindex;
    if (outcharindex >= outbuffers[outbufferindex].size) {
	++outbufferindex;
	outcharindex = 0;
    }
    return c;
}

bool ThreadStreamBuf::in_eof ()
{
    return inbufferindex >= (std::streamsize)(inbuffers.size());
}

bool ThreadStreamBuf::out_eof ()
{
    return outbufferindex >= (std::streamsize)(outbuffers.size());
}


PercyThreadedServer::PercyThreadedServer(ThreadedDataStore * ds) :
    PercyServer(ds)
{
    num_threads = ds->get_num_threads();
    for (dbsize_t i = 0; i < num_threads; ++i) {
	subservers.push_back(new PercyServer(ds->get_subds(i)));
    }
}

PercyThreadedServer::~PercyThreadedServer ()
{
    for (dbsize_t i = 0; i < num_threads; ++i) {
	delete subservers[i];
    }
}


struct ThreadParams {
    ThreadParams (PercyThreadedServer &thisserver, dbsize_t thread_index,
	    PercyServerParams &params, nqueries_t num_queries, 
	    unsigned char * const queries, unsigned char * const worker_responses) :
	thisserver(thisserver),
	thread_index(thread_index),
	params(params),
	num_queries(num_queries),
	queries(queries),
	worker_responses(worker_responses),
	exitstatus(true)
    {}

    PercyThreadedServer &thisserver;
    dbsize_t thread_index;
    PercyServerParams &params;
    nqueries_t num_queries;
    unsigned char * const queries;
    unsigned char * const worker_responses;
    bool exitstatus;
};

void * PThreadWork (void * arg) 
{
    ThreadParams * tparams = (ThreadParams*)arg;
    PercyThreadedServer &thisserver = tparams->thisserver;
    dbsize_t thread_index = tparams->thread_index;
    PercyServerParams &params = tparams->params;
    nqueries_t num_queries = tparams->num_queries;
    unsigned char * const queries = tparams->queries;
    unsigned char * const worker_responses = tparams->worker_responses;

    bool ret = thisserver.thread_work(thread_index, params, num_queries, queries,
	    worker_responses);
    tparams->exitstatus = ret;

    pthread_exit(NULL);
}

bool PercyThreadedServer::thread_work (dbsize_t thread_index, PercyServerParams &params, 
	nqueries_t num_queries, unsigned char * const queries, 
	unsigned char *	const worker_responses)
{
    // Get parameters
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t num_blocks = params.num_blocks();
    dbsize_t bytes_per_word = comm_bytes_per_word(params);
    PercyMode mode = params.get_mode();
    
    ThreadStreamBuf tsb;
    PercyServerParams subparams = worker_params(params, thread_index);

    // Add the number of queries for the worker to the buffer
    unsigned char nq[2];
    nqueries_t worker_nq = worker_num_queries(params, thread_index, num_queries);
    nq[0] = (worker_nq >> 8) & 0xff;
    nq[1] = (worker_nq) & 0xff;
    tsb.add_inbuffer((char*)nq, 2);

    // Add the worker's queries to the buffer.
    dbsize_t worker_qs = worker_query_bytes(params, thread_index);
    dbsize_t master_qs = bytes_per_word * num_blocks;
    if (mode == MODE_CHOR) {
	master_qs = num_blocks / 8;
    }
    dbsize_t worker_boffset = worker_bytes_offset(params, thread_index);
    nqueries_t worker_query_offset = worker_first_query(params, thread_index, num_queries);
    for (nqueries_t k = worker_query_offset; k < worker_query_offset + worker_nq; ++k) {
	tsb.add_inbuffer((char*)(queries + (master_qs * k) + worker_boffset), worker_qs);
    }

    // Add response array to stream
    tsb.add_outbuffer(
	    (char*)(worker_responses + worker_query_offset * words_per_block * bytes_per_word), 
	    worker_nq * words_per_block * bytes_per_word);

    // Handle request on subserver
    std::iostream sio(&tsb);
    bool ret = subservers[thread_index]->handle_request(subparams, sio, sio);

    // Check if buffer filled
    if (!(tsb.out_eof())) {
	fprintf(stderr, "[WORKER %lu] Thread did not finish properly\n", (long int)thread_index);
	return false;
    }

    return ret;
}

bool PercyThreadedServer::handle_request(PercyServerParams &params, std::istream &is, std::ostream &os)
{
    // Get parameters
    dbsize_t num_threads = params.get_num_threads();
    PercyThreadMethod tmethod = params.get_tmethod();
    dbsize_t words_per_block = params.words_per_block();
    dbsize_t num_blocks = params.num_blocks();
    dbsize_t bytes_per_word = comm_bytes_per_word(params);
    PercyMode mode = params.get_mode();

#ifdef VERBOSE_THREADED
    PercyThreadingType ttype = params.get_ttype();
    std::cerr << "Using threading type ";
    switch (ttype) {
    case THREADING_ROWS:
	std::cerr << "ROWS\n";
	break;
    case THREADING_QUERIES:
	std::cerr << "QUERIES\n";
	break;
    default:
	std::cerr << "***UNKNOWN***\n";
	break;
    }
    std::cerr << "Using threading method ";
    switch (tmethod) {
    case THREAD_METHOD_PTHREAD:
	std::cerr << "PTHREAD\n";
	break;
    case THREAD_METHOD_FORK:
	std::cerr << "FORK\n";
	break;
    case THREAD_METHOD_NONE:
	std::cerr << "NONE\n";
	break;
    default:
	std::cerr << "**UNKNOWN***\n";
	break;
    }
    std::cerr << "The number of threads is " << num_threads << "\n";
#endif

    if (mode == MODE_ZZ_P && tmethod == THREAD_METHOD_PTHREAD) {
	std::cerr << "Cannot use pthread method with ZZ_p\n";
	return false;
    }

    // Read the number of queries
    unsigned char nq[2];
    is.read((char*)nq, 2);
    if (is.eof()) {
	return false;
    }
    nqueries_t num_queries = (nq[0] << 8) | nq[1];

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
    dbsize_t master_qs = num_blocks * bytes_per_word;
    if (mode == MODE_CHOR) {
	master_qs = num_blocks / 8;
    }
    unsigned char * queries = new unsigned char[num_queries * master_qs];
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
	is.read((char*)(queries + (q * master_qs)), master_qs);
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
		ZZFromBytes(inputz, (queries + (master_qs * q) + (i * bytes_per_word)), 
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

    // Make response array
    std::vector<unsigned char *> worker_responses;
    for (dbsize_t i = 0; i < num_threads; ++i) {
	worker_responses.push_back(
		new unsigned char[num_queries *	words_per_block * bytes_per_word]);
	memset(worker_responses[i], 0, num_queries * words_per_block * bytes_per_word);
    }

    // Run the threads
    bool totalret = true;
    switch (tmethod) {
    case THREAD_METHOD_NONE: {
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    bool ret = thread_work(i, params, num_queries, queries,
		    worker_responses[i]);
	    if (!ret) {
		totalret = false;
	    }
	}
	delete[] queries;
    } break;

    case THREAD_METHOD_PTHREAD: {
	pthread_t * threads = new pthread_t[num_threads];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int rc;

	vector<ThreadParams*> tparams;
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    tparams.push_back(new ThreadParams(*this, i, params, num_queries,
		    queries, worker_responses[i]));
	    rc = pthread_create(&threads[i], &attr, PThreadWork, (void*)tparams[i]);
	    if (rc) {
		totalret = false;
	    }
	}

	pthread_attr_destroy(&attr);

	for (dbsize_t i = 0; i < num_threads; ++i) {
	    void * status;
	    rc = pthread_join(threads[i], &status);
	    if (rc) {
		std::cerr << "Error joining to thread\n";
		delete[] queries;
		for (dbsize_t i = 0; i < num_threads; ++i) {
		    delete[] worker_responses[i];
		}
#ifdef SPIR_SUPPORT
		delete [] spir_queries;
#endif
		return false;
	    }
	}

	delete[] queries;
	delete[] threads;
    } break;

    case THREAD_METHOD_FORK: {
	int * childfds = new int[num_threads];
	pid_t * pid = new pid_t[num_threads];
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    int fds[2];
	    if (pipe(fds)) {
		std::cerr << "Pipe failed\n";
		return false;
	    }
	    pid[i] = fork();
	    if (pid[i] < 0) {
		std::cerr << "Fork failed\n";
		delete [] childfds;
		delete [] pid;
		for (dbsize_t i = 0; i < num_threads; ++i) {
		    delete[] worker_responses[i];
		}
		delete [] queries;
#ifdef SPIR_SUPPORT
		delete [] spir_queries;
#endif
		return false;
	    } else if (pid[i] == 0) {
		// Child
		int parentfd = fds[1];
		close(fds[0]);
		// For compatibility with MPI
		close_highfds_except(parentfd);
		bool ret = thread_work(i, params, num_queries, queries,
			worker_responses[i]);
		// Send responses to parent
		ssize_t written = write(parentfd, (char*)worker_responses[i], 
			num_queries * words_per_block * bytes_per_word);
		if ((dbsize_t)written != num_queries * words_per_block * bytes_per_word) {
		    fprintf(stderr, "Did not write enough to pipe\n");
		    delete [] childfds;
		    delete [] pid;
		    for (dbsize_t i = 0; i < num_threads; ++i) {
			delete[] worker_responses[i];
		    }
		    delete [] queries;
#ifdef SPIR_SUPPORT
		    delete [] spir_queries;
#endif
		    ret = false;
		}
		close(parentfd);
		exit(ret);
	    } else {
		// Parent
		childfds[i] = fds[0];
		close(fds[1]);
	    }
	}
	delete [] queries;

	// Get responses from children
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    dbsize_t rr = read(childfds[i], (char*)worker_responses[i], num_queries * words_per_block * bytes_per_word);
	    if (rr < num_queries * words_per_block * bytes_per_word) {
		std::cerr << "Did not get entire response from thread " << i << "\n";
		delete [] childfds;
		delete [] pid;
		for (dbsize_t i = 0; i < num_threads; ++i) {
		    delete[] worker_responses[i];
		}
#ifdef SPIR_SUPPORT
		delete [] spir_queries;
#endif
		return false;
	    }
	}

	// Wait for children
	int status;
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    waitpid(pid[i], &status, 0);
	    if (WEXITSTATUS(status) != 0) {
		totalret = false;
	    }
	}

	delete[] childfds;
	delete[] pid;
    } break;
    }

    // Array to combine responses into
    unsigned char * responses = new unsigned char[num_queries * words_per_block * bytes_per_word];
    memset(responses, 0, num_queries * words_per_block * bytes_per_word);

    // Combine the worker responses
    switch (mode) {
    case MODE_ZZ_P: {
	vec_ZZ_p responses_ZZ_p;
	responses_ZZ_p.SetLength(num_queries * words_per_block);
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    for (dbsize_t w = 0; w < num_queries * words_per_block; ++w) {
		responses_ZZ_p[w] += to_ZZ_p(
			ZZFromBytes(worker_responses[i] + (w * bytes_per_word), bytes_per_word));
	    }
	}
	for (dbsize_t i = 0; i < num_queries * words_per_block; ++i) {
	    BytesFromZZ(responses + (i * bytes_per_word), rep(responses_ZZ_p[i]), bytes_per_word);
	}
    } break;

    case MODE_GF28:
    case MODE_GF216:
    case MODE_CHOR:
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    XOR_equal(responses, worker_responses[i], 
		    num_queries * words_per_block * bytes_per_word);
	}
	break;
    }

#ifdef SPIR_SUPPORT
    if (params.spir()) {
	vec_ZZ_p response;
	response.SetLength(words_per_block);
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    for (dbsize_t i = 0; i < words_per_block; ++i) {
		ZZ single;
		ZZFromBytes(single, (responses + ((q * words_per_block) + i) * bytes_per_word),
			bytes_per_word);
		response[i] = to_ZZ_p(single);
	    }
	    spir_queries[q].randomize_response(response);
	    for (dbsize_t i = 0; i < words_per_block; ++i) {
		BytesFromZZ(responses + ((q * words_per_block) + i) * bytes_per_word, 
			rep(response[i]), bytes_per_word);
	    }
	}
    }
#endif

    // Send the responses
    os.write((char*)responses, num_queries * words_per_block * bytes_per_word);
    os.flush();

    // Free the arrays
#ifdef SPIR_SUPPORT
    delete[] spir_queries;
#endif
    for (dbsize_t i = 0; i < num_threads; ++i) {
	delete[] worker_responses[i];
    }
    delete[] responses;

    return totalret;
}

