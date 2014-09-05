//  Percy++ Copyright 2007,2012,2013 Ian Goldberg <iang@cs.uwaterloo.ca>,
//  Casey Devet <cjdevet@cs.uwaterloo.ca>,
//  Paul Hendry <pshdenry@uwaterloo.ca>,
//  Ryan Henry <rhenry@cs.uwaterloo.ca>,
//  Femi Olumofin <fgolumof@cs.uwaterloo.ca>
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
#include <iostream>
#include <fstream>
#include <sstream>
#include "datastore.h"
#include "percyserver.h"
#include "percyparams.h"
#include "config.h"
#include <sys/types.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <socket++/sockinet.h>
#include <unistd.h>
#include <getopt.h>
#include "distserver.h"
#include "threadedserver.h"
#include "cmdtools.h"

#define PERCY_DEFAULT_PORT 31337
#define PERCY_MAX_BIND_ATTEMPTS 10

// For now, no MPI
#ifdef MPI_DIST_SERVER
#error MPI is not currently supported!
#endif

#ifdef MPI_DIST_SERVER
#include <mpi.h>
#include "mpicomm.h"
#include "mpiserver.h"
#endif

// Global variables
// Distributed server information

#define PERCY_MAX_CONNECT_ATTEMPTS 5

// A method to connect to a server.  Returns a pointer to a socket if successful
// and a NULL pointer if unsuccessful.
iosockinet * connect_to_server (const char * addr, const uint16_t port) {
    std::cerr << "    Attempting to connect to " << addr << ":" << port << "...";
    bool connected = false;
    unsigned short attempts = 0;
    iosockinet *socket = new iosockinet(sockbuf::sock_stream);
    while (!connected && (attempts++ < PERCY_MAX_CONNECT_ATTEMPTS))
    {   
        try 
        {   
            (*socket)->connect(addr, port);
            connected = true;
        }   
        catch (sockerr e)
        {   
            cerr << ".";
            sleep(1);
        }   
    }   
    if (connected) {
        std::cerr << "succeeded!" << std::endl;
        return socket;
    } else {
        std::cerr << "failed!" << std::endl;
        return NULL;
    }
}

void print_usage_options () {
    std::cerr << "Mandatory arguments to long options are mandatory for short options too." << std::endl;
    std::cerr << "   -n DBBYTES             use only the first DBBYTES bytes of database (default: entire file)." << std::endl;
    std::cerr << "   -w WORDSIZE            use a word size of WORDSIZE bytes (default: 8)." << std::endl;
    std::cerr << "   -b BLOCKSIZE           use a block size of BLOCKSIZE bytes (default: sqrt(DBBYTES*WORDSIZE)/8)." << std::endl;
    std::cerr << "   -u MAX_UNSYNCHRONIZED  specifies the maximum number of db files that can be unsynchronized." << std::endl;
    std::cerr << "   -e EXPANSION_FACTOR    specifies the expansion factor in the number of bins to produce for synchronization." << std::endl;
    std::cerr << "   -t, -tau               specify that database is tau independent." << std::endl;
    std::cerr << "   -S, --SID SERVERID     use the specified SID." << std::endl;
    std::cerr << "   -p, --port PORTNO      listen for connections on the specified port." << std::endl;
    std::cerr << "   -m, --mode MODE        use the specified mode of operation. Supported modes are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          GF28        g            use fast arithmetic in GF(2^8)" << std::endl;
    std::cerr << "                          GF216       s            use fast arithmetic in GF(2^16)" << std::endl;
    std::cerr << "                          ZZ_P        z            use arithmetic in Z mod p" << std::endl;
    std::cerr << "                          CHOR        c            use Chor et al.'s lightweight protocol" << std::endl;
    std::cerr << "                          RS_SYNC     r            use unsynchronized database scheme with Reed-Solomon decoding" << std::endl;
    // std::cerr << "                          PULSE_SYNC  p            use unsynchronized database scheme with PULSE decoding" << std::endl;    
    std::cerr << "                          (default: ZZ_P)" << std::endl;
#ifdef SPIR_SUPPORT
    std::cerr << "   -s, --spir PCPARAMS    do symmetric PIR with specified PolyCommit" << std::endl;
#endif
    std::cerr << "                          parameters (a file)." << std::endl;
    std::cerr << "   -h, --hybrid           support hybrid security." << std::endl;
    std::cerr << "   -z, --byzantine        be byzantine." << std::endl;
    std::cerr << "   -1, --oneconn          accept only a single conncetion; do not fork." << std::endl;
    std::cerr << "       --help             display this help and exit." << std::endl;
    std::cerr << "       --version          output version information and exit." << std::endl;
    std::cerr << std::endl;
#ifndef DIST_MASTER
    std::cerr << "Distributed Server Options:" << std::endl;
    std::cerr << "   -T, --num-threads T    Distribute computation over T threads." << std::endl;
    std::cerr << "   -P, --thread-type P    Specify how the queries are split up between threads.  Supported" << std::endl;
    std::cerr << "                          types are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          row         r            each thread is assigned a subset of database rows" << std::endl;
    std::cerr << "                          queries     q            each thread is assigned a subset of the queries" << std::endl;
    std::cerr << "                          (default: row)" << std::endl;
    std::cerr << "   -Q, --thread-method Q  Specify the threading method to be used.  Supported methods are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          pthread     p            thread using the pthread library (not compatible" << std::endl;
    std::cerr << "                                                   with ZZ_P mode)" << std::endl;
    std::cerr << "                          fork        f            each thread is actually a forked child" << std::endl;
    std::cerr << "                          none        n            run the workers in series in one thread" << std::endl;
    std::cerr << "                          (default: 'pthread', except 'fork' when using ZZ_P)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Query Input/Response Output Options:" << std::endl;
#ifdef MPI_DIST_SERVER
    std::cerr << "   -M, --mpi-tests TESTS          Simulates TESTS number of tests by a client over MPI." << std::endl;
#endif
/* TODO
    std::cerr << "   -F, --queries-from-file FILE   Use the file FILE as the data sent from the client." << std::endl;
    std::cerr << "   -G, --responses-to-file FILE   Redirect the responses to FILE." << std::endl;
*/
    std::cerr << std::endl;
#endif

    std::cerr << "Report bugs to iang+percy@cs.uwaterloo.ca." << std::endl;
}

void print_usage_normal (const char * bin) {
    std::cerr << "Usage: " << bin << " [OPTIONS...] DATABASE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Start the PIR server on the specified database file DATABASE." << std::endl;
    std::cerr << std::endl;

    print_usage_options();
}

void print_usage_dist (const char * bin) {
    std::cerr << "Usage: " << bin << " [OPTIONS...] VSPLIT HSPLIT WORKERINFO" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Start a PIR server as a master for the specified workers." << std::endl;
    std::cerr << std::endl;

    // Explanation of arguments for Master and MPI usage
    std::cerr << "Arguments:" << std::endl;
    std::cerr << "   VPSLIT       The database is split into VSPLIT sets of rows." << std::endl;
    std::cerr << "   HSPLIT       The database is split into HSPLIT sets of columns." << std::endl;
    std::cerr << "   WORKERINFO   The addresses and port numbers for each worker.  Must be of the form" << std::endl;
    std::cerr << "                \"sid0:addr0:port0 sid1:addr1:port1 ... sidM:addrM:portM\" where M=VSPLIT*HSPLIT." << std::endl;
    std::cerr << "                The order is important as index i represents the subdatabase at position" << std::endl;
    std::cerr << "                (i % VSPLIT, i / HSPLIT) in the matrix of subdatabases." << std::endl;
    std::cerr << std::endl;

    print_usage_options();
}

void print_usage (const char * bin) {
#if defined(DIST_MASTER)
    print_usage_dist(bin);
#elif defined(MPI_DIST_SERVER)
    print_usage_mpi(bin);
#else
    print_usage_normal(bin);
#endif
}

// Get the database size.  Returns 0 if the database does not exist.
dbsize_t database_bytes (const char * database) {
    struct stat filestatus;
    int not_exists = stat(database, &filestatus);
    if (not_exists)
    {
        std::cerr << "Error: cannot find database file " << database << std::endl;
        return 0;
    }
    return filestatus.st_size;
}

PercyServerParams * init_params(ParsedArgs& pargs, bool checkdb = true)
{
    // Allow hybrid queries?
    bool do_hybrid = false;
    
    // Mode of operation selected (ZZ_p, GF28, GF216 or Chor)
    PercyMode mode = MODE_ZZ_P;

    // Do symmetric PIR?
    bool do_spir = false;

    // PolyCommit Params file to read for SPIR
    char *pcparams_file = NULL;
    
    // Is the database tau independent?
    bool is_tau = false;

    // Should we be byzantine?
    bool be_byzantine = false;

    // server must know its SID if and only if it supports spir or the
    // database is tau independent. A SID of 0 means "unknown SID".
    nservers_t sid = 0;
    
    dbbits_t n_bytes = 0;
    dbsize_t w = 0;
    dbsize_t b = 0;
    dbsize_t max_unsynchronized = 0; //max number of unsynchronized records in the database
    dbsize_t expansion_factor = 1; //max number of unsynchronized records in the database

    // Threading parameters
    dbsize_t num_threads = 0;
    PercyThreadingType ttype = THREADING_ROWS;
    PercyThreadMethod tmethod = THREAD_METHOD_PTHREAD;

    optmap::iterator oiter;
    for (oiter = pargs.opts.begin(); oiter != pargs.opts.end(); ++oiter) {
	char * optarg = oiter->second;
        switch(oiter->first) {
            case 'n':
                n_bytes = strtoull(optarg, NULL, 10);
                break;
            case 'b':
                b = strtoull(optarg, NULL, 10) * 8;
                break;
            case 'w':
                w = strtoull(optarg, NULL, 10);
                break;
            case 'u':
                max_unsynchronized = strtoull(optarg, NULL, 10);
                break;
            case 'e':
                expansion_factor = strtoull(optarg, NULL, 10);
                break;
            case 't':
                is_tau = true;
                break;
            case 'm':
                if(!strcmp(optarg, "ZZ_P") || !strcmp(optarg, "z")) {
                    mode = MODE_ZZ_P;
                }
                else if(!strcmp(optarg, "GF28") || !strcmp(optarg, "g")) {
                    mode = MODE_GF28;
                }
                else if(!strcmp(optarg, "GF216") || !strcmp(optarg, "s")) {
                    mode = MODE_GF216;
                }
                else if(!strcmp(optarg, "CHOR") || !strcmp(optarg, "c")) {
                    mode = MODE_CHOR;
                }
                else if(!strcmp(optarg, "RS_SYNC") || !strcmp(optarg, "r")) {
                    mode = MODE_RS_SYNC;
                }
                else {
                    std::cerr << "Unknown mode selected. Valid modes are ZZ_P, GF28, GF216 and CHOR.\n\n";
		    print_usage(pargs.exec);
		    return NULL;
                }
                break;
            case 'h':
                do_hybrid = true;
                break;
            case 'z':
                be_byzantine = true;
                break;
#ifdef SPIR_SUPPORT
            case 's':
                do_spir = true;
                pcparams_file = optarg;
                break;
#endif
            case 'S':
                sid = strtoul(optarg, NULL, 10);
                //std::cerr << "SID is " << sid << "." << std::endl;
                break;
            case 'v':
                std::cerr << "Percy++ pirserver version " << VERSION << std::endl;
                std::cerr << AUTHOR << std::endl;
		return NULL;
                break;
	    case 'T':
		num_threads = strtoll(optarg, NULL, 10);
		if (num_threads < 1) {
		    fprintf(stderr, "Must specify at least 1 threads!\n");
		    return NULL;
		}
		break;
	    case 'P':
		if (!strcmp(optarg, "row") || !strcmp(optarg, "r")) {
		    ttype = THREADING_ROWS;
		} else if (!strcmp(optarg, "queries") || !strcmp(optarg, "q")) {
		    ttype = THREADING_QUERIES;
		} else {
		    std::cerr << "Invalid threading type selected.\n";
		    print_usage(pargs.exec);
		    return NULL;
		}
		break;
	    case 'Q':
		if (!strcmp(optarg, "pthread") || !strcmp(optarg, "p")) {
		    tmethod = THREAD_METHOD_PTHREAD;
		} else if (!strcmp(optarg, "fork") || !strcmp(optarg, "f")) {
		    tmethod = THREAD_METHOD_FORK;
		} else if (!strcmp(optarg, "none") || !strcmp(optarg, "n")) {
		    tmethod = THREAD_METHOD_NONE;
		} else {
		    std::cerr << "Invalid threading method selected.\n";
		    print_usage(pargs.exec);
		    return NULL;
		}
		break;
            default:
            // Invalid options are handled in argument parsing (no error or
            // exit here)
            break;
        }
    }

    // Change threading method from pthread to fork if in ZZ_P
    if (num_threads > 0 && mode == MODE_ZZ_P && tmethod == THREAD_METHOD_PTHREAD) {
        fprintf(stderr, "The pthread library is not compatible with ZZ_p.  Using the fork method instead.\n");
        tmethod = THREAD_METHOD_FORK;
    }
    
    // TODO: Check compatible with DIST_SERVER
    if (do_hybrid && (mode != MODE_ZZ_P)) {
        fprintf(stderr, "Error: hybrid security can only be used with the integers mod p mode of operation.\n");
        return NULL;
    }
#ifdef SPIR_SUPPORT
    if (do_hybrid && do_spir) {
        fprintf(stderr, "Error: cannot use hybrid security with symmetric PIR.\n");
        return NULL;
    }
    if (do_spir && mode != MODE_ZZ_P) {
        fprintf(stderr, "Error: symmetric PIR can only be used with the integers mod p mode of operation.\n");
        return NULL;
    }
#endif
    if (is_tau && mode == MODE_CHOR) {
        fprintf(stderr, "Error: Chor et al.'s PIR scheme does not support tau independence.\n");
        return NULL;
    }

    // Make sure enough mandatory arguments are present.
    if (pargs.nonoptc < 1) {
        fprintf(stderr, "Not enough arguments\n");
        print_usage(pargs.exec);
        return NULL;
    }

    if (checkdb) {
        // Make sure the specified database file exists.
        const char *database = pargs.nonoptv[0];
        dbsize_t dbsize = database_bytes(database);
        if (dbsize == 0) {
            fprintf(stderr, "Error: the database must exist and be non-empty.\n");
            return NULL;
        }

        // If no value for "n" is specified, then use a default database
        // size of dbsize. Otherwise, just check that 0<n<=dbsize.
        if (!n_bytes) {
            n_bytes = dbsize;
        } else if (n_bytes > dbsize) {
            fprintf(stderr, "Error: n cannot be larger than database file.\n");
            return NULL;
        }

    } else {
        // If not checking database, database size must be specified
        if (!n_bytes) {
            fprintf(stderr, "Error: Database size (n) must be specified.\n");
            return NULL;
        }
    }

    dbbits_t n = n_bytes * 8;
    if (n_bytes > n)
    {
        fprintf(stderr, "Error: database file is too large for the current architecture!\n");
        return NULL;
    }

    // If no value for "w" is specified, then use a default word size
    // of 8 bits.
    if (!w)
    {
        w = 8;
    }
    
    // If no value for "b" is specified, then use a default block size
    // of \sqrt(n * w) bits.
    if (!b)
    {
        b = sqrt(n * w);
        if (n != b*b/w)
        {
            fprintf(stderr, "Error: optimal parameter choice is invalid for this database. Please specify a value for both of b and w.\n");
	    return NULL;
        }
    }
    
    //std::cerr << "Debug: (n,b,w) = (" << n << "," << b << "," << w << ") " << std::endl;
    
    // Sanity checks for (n,b,w).
    if (n % b != 0 || b % w != 0)
    {
        fprintf(stderr, "Error: b must divide n and w must divide b.\n");
	return NULL;
    }
    if (mode == MODE_CHOR) {
        if (w != 1) {
            fprintf(stderr, "Error: w must be 1 in Chor et al.'s PIR scheme.\n");
	    return NULL;
        }
    }
    else {
        if (w % 8 != 0) {
            fprintf(stderr, "Error: 8 must divide w.\n");
	    return NULL;
        }
    }
    if (mode == MODE_GF28 && w != 8) {
        fprintf(stderr, "Error: w must be 8 for gf28.\n");
	return NULL;
    }
    
    // Compute the number of blocks, and number of words per block.
    //std::cerr << "b = " << b << "\n";
    dbsize_t num_blocks = n / b;
    dbsize_t words_per_block = b / w;
    
    //std::cerr << "Number of blocks: " << num_blocks << std::endl;
    //std::cerr << "Words per block:  " << words_per_block << std::endl;
    //std::cerr << "Bits per block:  " << b << std::endl;
    //std::cerr << "Bits per word:  " << w << std::endl;
    if (num_blocks != words_per_block)
    {
        std::cerr << "Warning: non-optimal choice of blocksize detected." << std::endl;
    }

    // Choose an appropriate modulus.
    ZZ modulus;
    if (w == 2048)
    {
        modulus = to_ZZ("51162405833378812589599605953260132300166393994651819099454781579567509212081792013783783759303440508155949594262147212874957344953142209597742684263402581129339826752613431877280173074502314648334418584122460414512816448592261381117519846844295394134225624418756277265452922709245846828145574822031541004633366879073894273715489429502290966133193310966178373909137394353164436844312924586836474134940807305776164928781025210917912257206480517698118422827367766257579221703667784216949825206167241852365543481875593117676222875888924950402025039269210778276794873837063438751454865130720887819939394489366347567251243");
    }
    else if (w == 1536)
    {
        modulus = to_ZZ("4065256781338999183533854850423382625119065920051798531476300569026463202897155088318466013703570859212040475097762405522038651420119366364979939687154236065682459920101982590074846996306687236388206057475890613264408059472973401701686869808348910896596468985609043697525749128687318350246421674945679872669881805678484464202726328189280359385791023305618545788872763420795247846720674554774196715770302797683129209164871258189464484019233379849839076263862630987");
    }
    else if (w == 1024)
    {
        modulus = to_ZZ("343308946066366926839932845260501528909643718159825813630709694160026342456154871924497152436552679706642965502704642456637620829912957820221098686748075257358288200837461739492534713539606088624083011849535450485951774635526473457667739540374042376629835941950802202870595346459371144019363420985729553740241");
    }
    else if (mode == MODE_GF28)
    {
        modulus = to_ZZ("256");
    }
    else if (mode == MODE_GF216)
    {
        modulus = to_ZZ("65536");
    }
    else if (mode == MODE_CHOR)
    {
        //Important: for Chor we pretend as though a word is 1 byte. This is because many 
        //parts of the code rely on a word being a byte multiple (for example, where bytes_per_word 
        //is used). We set this here since the calculations for the optimal database shape need a 
        //word size of 1 bit.
        words_per_block /= 8;
        modulus = to_ZZ("256");
    }
    else if (w == 8 && !do_hybrid)
    {
        modulus = to_ZZ("257");
    }
    else if (w == 16 && !do_hybrid)
    {
        modulus = to_ZZ("65537");
    }
    else if (w == 32 && !do_hybrid)
    {
        modulus = to_ZZ("4294967311");
    }
    else if (w == 96 && !do_hybrid)
    {
        modulus = to_ZZ("79228162514264337593543950397");
    }
    else if (w == 128 && !do_hybrid)
    {
        modulus = to_ZZ("340282366920938463463374607431768211507");
    }
    else if (w == 160 && !do_hybrid)
    {
        // NOTE: p2s is the prime from the PolyCommit params; spir
        // will break if this value gets changed!
        //
        // TODO: read the prime from the PolyCommit params and check
        //          that it is consistent with w.
        modulus = to_ZZ("2425980306017163398341728799446792216592523285797");
    }
    else if (w == 192 && !do_hybrid)
    {
        modulus = to_ZZ("6277101735386680763835789423207666416102355444464034513029");
    }
    else if (w == 256 && !do_hybrid)
    {
        modulus = to_ZZ("115792089237316195423570985008687907853269984665640564039457584007913129640233");
    }
    else if (do_hybrid)
    {
        std::cerr << "Error: No hybrid-compatible modulus available for w = " << w << "." << std::endl;
	return NULL;
    }
    else
    {
        std::cerr << "Error: No modulus available for w = " << w << "." << std::endl;
	return NULL;
    }
#ifdef SPIR_SUPPORT
    if (do_spir && w!=160)
    {
        fprintf(stderr, "Error: symmetric PIR currently supports only w=160.\n");
	return NULL;
    }
#endif

    // Create the PercyServerParams object.
    PercyServerParams * params = new PercyServerParams(
	    words_per_block, num_blocks, max_unsynchronized, expansion_factor, is_tau, modulus, mode, 
	    be_byzantine, pcparams_file, do_spir, sid,
	    num_threads, ttype, tmethod);

    return params;
}


PercyDistServerParams * init_dist_params (ParsedArgs& pargs)
{
    // Get basic server params
    PercyServerParams * sparams = init_params(pargs, false);
    if (sparams == NULL) {
	return NULL;
    }

    // Check that there are at least 3 non-option arguments
    if (pargs.nonoptc < 3) {
	fprintf(stderr, "Not enough arguments\n");
	print_usage(pargs.exec);
	return NULL;
    }

    // Get splits
    dbsize_t vsplit = strtoull(pargs.nonoptv[0], NULL, 10);
    dbsize_t hsplit = strtoull(pargs.nonoptv[1], NULL, 10);
    if (vsplit == 1 && hsplit == 1) {
	std::cerr << "Warning: Splitting the database 1x1.\n";
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

    // Create dist params
    PercyDistServerParams * params = new PercyDistServerParams(*sparams,
	    vsplit, hsplit);
    delete sparams;

    // Get worker addresses and ports
    std::vector<char*> workerstr;
    char * token = strtok(pargs.nonoptv[2], " ");
    while (token != NULL) {
	workerstr.push_back(token);
	token = strtok(NULL, " ");
    }
    if (workerstr.size() != vsplit * hsplit) {
	std::cerr << "Information was given for " << workerstr.size()
		<< " workers.  Need " << vsplit + hsplit << " workers\n";
	return NULL;
    }

    // Create worker iostreams
    for (dbsize_t i = 0; i < vsplit * hsplit; ++i) {
	char * sidchar = strtok(workerstr[i], ":");
	char * addrchar = strtok(NULL, ":");
	char * portchar = strtok(NULL, ":");
	if (addrchar == NULL || portchar == NULL || portchar == NULL) {
	    fprintf(stderr, "Error: worker information was incorrectly formatted.\n");
	    return NULL;
	}
	unsigned long sidnum = strtoul(sidchar, NULL, 10);
	if (!sidnum || sidnum > params->get_modulus()) {
	    std::cerr << "Error: SID must be an integer greater than 0 and less than " << params->get_modulus() << ".\n";
	    return NULL;
	}
	unsigned long portnum = strtoul(portchar, NULL, 10);
	if (portnum < 1024 || portnum > 65535) {
	    fprintf(stderr, "Error: port number must be an integer greater than 1024 and less than 65535.\n");
	    return NULL;
	}
	// connect_to_server function is located in distserver.{h,cc}
	iosockinet * socket = connect_to_server(addrchar, portnum);
    std::cerr << "Hi there\n";
	if (socket == NULL) {
	    std::cerr << "Error: cannot connect to worker " << i << ".\n";
	    return NULL;
	}
	std::iostream * stream = socket;
	params->set_workerio(i, socket);

	// Send params to worker
	unsigned char failure;
	*stream << params->sub_params();
	// NOTE: We assume that all workers are not doing SPIR
	//       If they are, the test below will fail since sub_params returns
	//       params with do_spir=false.
	// Send worker sid if needed.
	if (params->tau()) {
	    unsigned char sidc = sidnum & 0xff;
	    stream->write((char*)(&sidc), 1);
	}
	stream->flush();
	// Check for failure by worker
	stream->read((char*)(&failure), 1);
	if (failure & 1) {
	    std::cerr << "Error: Worker " << i << " did not accept parameters.\n";
	    return NULL;
	} else if (failure & 2) {
	    std::cerr << "Error: Worker " << i << " did not accept sid.\n";
	    return NULL;
	}
    }

    return params;
}

// Initializes and allocates the datastore
// If num_threads > 1, a ThreadedDataStore is created.  Otherwise, a
// FileDataStore is created.
DataStore * init_datastore (PercyServerParams& serverparams, const char * database)
{
    dbsize_t dbsize = database_bytes(database);
    if (dbsize == 0) {
	fprintf(stderr, "Error: The database must exist and be non-empty\n");
	return NULL;
    }

    DataStore * ds;
    if (serverparams.get_num_threads() > 0) {
	ds = new ThreadedDataStore(database, serverparams, serverparams.tau());
    } else {
	ds = new FileDataStore(database, serverparams, serverparams.tau());
    }
    return ds;
}


// Handle all of the requests on a single connection
static void handle_requests(PercyServer * server, std::istream &is, std::ostream &os,
    PercyServerParams &serverparams)
{
    std::cerr << "Received an incoming connection." << std::endl;
    // Receive the parameters from the client (and the SID, if necessary).
    unsigned char failure = 0;
    nservers_t sid = serverparams.get_sid();

    // First, read the client's query parameters.
    PercyParams clientparams;
    //                std::cerr << "Receiving query parameters from client...";
    is >> clientparams;
    //                std::cerr << "done" << std::endl;
    // check if the two sets of params are compatible...
    failure = !serverparams.is_compatible(clientparams);
    if (failure) {
        std::cerr << "Client is not compatible with server\n";
    }
    
    // Now do the SID, but only if it's needed
    if (clientparams.spir() || clientparams.tau())
    {
        nservers_t client_sid;
        //                    std::cerr << "Receiving SID from client...";
        unsigned char sidc;
        is.read((char*)&sidc, 1);
        client_sid = (nservers_t) sidc;
        //                    std::cerr << "done" << std::endl;
        if (client_sid != sid)
        {
            std::cerr << "Received incorrect SID from client. (Expected " << sid << " but saw " << client_sid << ".)" << std::endl;
            failure &= 2;
        }
    }
    os.write((char*)&failure, 1);
    os.flush();
    
    // Finally, do the PIR query!
    // With probability $PIRS_FAIL/100, fail completely
    // With probability $PIRS_BYZ/100, be Byzantine
    unsigned long rndval = RandomBnd(100);
    unsigned long failat = 0, byzat = 0, byznum = 0;
    const char *failenv = getenv("PIRS_FAIL");
    if (failenv) failat = atoi(failenv);
    const char *byzenv = getenv("PIRS_BYZ");
    if (byzenv) byzat = atoi(byzenv);
    const char *byznumenv = getenv("PIRS_BYZN");
    if (byznumenv) byznum = atoi(byznumenv);

    if (rndval < failat)
    {
        std::cerr << "["<<sid<<"] Failing.\n";
        return;
    }
    if (serverparams.is_byzantine() || rndval < failat + byzat || sid <= byznum)
    {
        std::cerr << "["<<sid<<"] Going Byzantine.\n";
        server->be_byzantine();
    }

    // Handle the request(s)
    // This loop will run until all queries are read and eof is
    // read.  We then gracefully exit the child process.
    struct timeval ts, te;
    gettimeofday(&ts, NULL);
    while (server->handle_request(serverparams, is, os)) {
        gettimeofday(&te, NULL);
        int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
        fprintf(stderr, "%d.%03d msec computation + communication\n", td/1000, td%1000);
        gettimeofday(&ts, NULL);
    }
    std::cerr << "Query completed." << std::endl;
}

bool bind_to_port (sockinetbuf& sin, uint16_t& port)
{
    sin.reuseaddr(true);
    if (!port)
    {
        port = PERCY_DEFAULT_PORT;
        retry:
        try
        {
            sin.bind((unsigned long) INADDR_ANY, port);
        }
        catch (sockerr)
        {
            //std::cerr << "Debug: failed to bind to port " << port << "."<< std::endl;
            port++;
            if (port < PERCY_DEFAULT_PORT + PERCY_MAX_BIND_ATTEMPTS)
            {
                goto retry;
            }
            else
            {
                std::cerr << "Error: unable to bind socket to a port. (I tried using ports " << PERCY_DEFAULT_PORT << " through " << (PERCY_DEFAULT_PORT + PERCY_MAX_BIND_ATTEMPTS) << ".)" << std::endl;
		return false;
            }
        }
    }
    else
    {
        try
        {
            sin.bind((unsigned long) INADDR_ANY, port);
        }
        catch (sockerr)
        {
            std::cerr << "Error: unable to bind socket on port " << port << "." << std::endl;
	    return false;
        }
    }
    sin.listen();
    std::cerr << "Listening on port " << port << "." << std::endl;
    return true;
}

void init_NTL_and_rand ()
{
    // Initialize NTL and the random number stream
    ZZ modinit;
    modinit = to_ZZ(257);
    ZZ_p::init(modinit);
    unsigned char randbuf[128];
    ifstream urand("/dev/urandom");
    urand.read((char *)randbuf, sizeof(randbuf));
    urand.close();
    ZZ randzz = ZZFromBytes(randbuf, sizeof(randbuf));
    SetSeed(randzz);
}

// List of long (i.e. --) options:
// {"longoptname", no_argument|required_argument|optional_argument, 0, 'shortoptname'}
struct option longopts[] = {
    {"version",		    no_argument,	NULL, 'v'},
    {"help",		    no_argument,	NULL, 'a'},
    {"n",		    required_argument,	NULL, 'n'},
    {"b",		    required_argument,	NULL, 'b'},
    {"w",		    required_argument,	NULL, 'w'},
    {"u",		    required_argument,	NULL, 'u'},
    {"e",		    required_argument,	NULL, 'e'},
    {"tau",		    no_argument,	NULL, 't'},
    {"mode",		    required_argument,  NULL, 'm'},
    {"hybrid",		    no_argument,        NULL, 'h'},
    {"byzantine",	    no_argument,        NULL, 'z'},
    {"oneconn",		    no_argument,        NULL, '1'},
#ifdef SPIR_SUPPORT
    // spir: must specify PolyCommit parameters.
    {"spir",		    required_argument,  NULL, 's'},
#endif
    // SID: must specify a server ID.
    {"SID",		    required_argument,  NULL, 'S'},
    // port: must specify a port number.
    {"port",		    required_argument,  NULL, 'p'},
#ifdef MPI_DIST_SERVER
    {"mpi-tests",	    required_argument,	NULL, 'M'},
#endif
    {"num-threads",	    required_argument,	NULL, 'T'},
    {"thread-type",	    required_argument,  NULL, 'P'},
    {"thread-method",	    required_argument,  NULL, 'Q'},
    {"queries-from-file",   required_argument,	NULL, 'F'},
    {"responses-to-file",   required_argument,	NULL, 'G'},
    {NULL,		    0,                  NULL, 0},
};



#ifndef MPI_DIST_SERVER
// Normal main
int main (int argc, char ** argv)
{       
    // Ignore SIGPIPE
    signal(SIGPIPE, SIG_IGN);

    init_NTL_and_rand();

    // Parse arguments

    #if defined(DIST_MASTER) & defined(SPIR_SUPPORT)
        const char * shortopts = "d:n:b:u:e:w:tm:hz1S:p:F:G:s:";
    #elif defined(DIST_MASTER)
        const char * shortopts = "d:n:b:u:e:w:tm:hz1S:p:F:G:";
    #elif defined(SPIR_SUPPORT)
        const char * shortopts = "d:n:b:u:e:w:tm:hz1S:p:F:G:s:T:P:Q:";
    #else
        const char * shortopts = "d:n:b:u:e:w:tm:hz1S:p:F:G:T:P:Q:";
    #endif
    ParsedArgs pargs;
    if (!parse_long_opts(argc, argv, shortopts, longopts, pargs)) {
        print_usage(argv[0]);
        return -1;
    }

    // Check for help or version flags
    optmap::iterator oiter = pargs.opts.find('a');
    if (oiter != pargs.opts.end()) {
        print_usage(pargs.exec);
        return 0;
    }
    oiter = pargs.opts.find('v');
    if (oiter != pargs.opts.end()) {
        std::cerr << "Percy++ pirserver version " << VERSION << std::endl;
        std::cerr << AUTHOR << std::endl;
        return 0;
    }

    // Initialize the parameters
    #ifdef DIST_MASTER
        PercyDistServerParams * params = init_dist_params(pargs);
    #else
        PercyServerParams * params = init_params(pargs);
    #endif
    if (params == NULL) {
        return -1;
    }
    
    // Create datastore
    DataStore * datastore = NULL;
#ifndef DIST_MASTER
    datastore = init_datastore(*params, pargs.nonoptv[0]);
    if (datastore == NULL) {
        delete params;
        fprintf(stderr, "DataStore was not initialized.\n");
        return -1;
    }
#endif

    // Get port
    uint16_t port = 0;
    oiter = pargs.opts.find('p');
    if (oiter != pargs.opts.end()) {
        port = strtoul(oiter->second, NULL, 10);
    }

    // Create a socket for clients to connect to.
    sockinetbuf sin(sockbuf::sock_stream);
    if (!bind_to_port(sin, port)) {
        delete params;
        if (datastore != NULL) {
            delete datastore;
        }
        fprintf(stderr, "Did not successfully bind to port %d.\n", port);
        return -1;
    }

    // Create the server
    PercyServer * server = NULL;
    #ifdef DIST_MASTER
        server = new PercyMasterServer();
    #else
        if (params->get_num_threads() > 0) {
            server = new PercyThreadedServer(static_cast<ThreadedDataStore*>(datastore));
        } else {
            // Create the PIR server
            server = new PercyServer(datastore);
            if (port == 31338) {
                server->set_server_unsynchronized(*params);
            }
        }
    #endif
    if (server == NULL) {
        delete params;
        if (datastore != NULL) {
            delete datastore;
        }
        fprintf(stderr, "Server not created successfully\n");
        return -1;
    }

    // Get daemon_mode
    bool daemon_mode = true;
    oiter = pargs.opts.find('1');
    if (oiter != pargs.opts.end()) {
        daemon_mode = false;
    }

    // Get redirection options
    // TODO:

    if (daemon_mode) {
        // Daemon mode
        while(true) {
            iosockinet sio(sin.accept());

            pid_t childpid = fork();
            if (childpid) {
                waitpid(childpid, NULL, 0);
            } else {
                // spawn a grandchild and commit suicide so that the
                // parent doesn't have to wait()
                pid_t grandchildpid = fork();
                if (grandchildpid) {
                    break; // Will exit loop, clean up and exit
                } else {
                    // Handle request
                    handle_requests(server, sio, sio, *params);
                    break; // Will exit loop, clean up and exit
                }
            }
        }
    } else {
        // One connection
        /*
        std::ifstream ifs("test.in");
        std::ofstream ofs("test.out");
        std::istream is(ifs.rdbuf());
        std::ostream os(ofs.rdbuf());
        handle_requests(is, os, serverparams, datastore, be_byzantine);
        */
        // Get incoming socket connection.
        iosockinet sio(sin.accept());

        // Handle request
        handle_requests(server, sio, sio, *params);
    }

    // Clean up
    delete params;
    if (datastore != NULL) {
        delete datastore;
    }
    delete server;

    return 0;
}

#else
// MPI main
int main(int argc, char **argv)
{

    std::cerr << "in MPI main\n";
    // Ignore SIGPIPE
    signal(SIGPIPE, SIG_IGN);

    MPI_Init(&argc, &argv);

    init_NTL_and_rand();

    // Parse arguments
    const char * shortopts = "d:n:b:w:tm:hz1S:p:F:G:T:P:Q:M:";
    ParsedArgs pargs;
    if (!parse_long_opts(argc, argv, shortopts, longopts, pargs)) {
	print_usage(argv[0]);
	MPI_Finalize();
	return -1;
    }

    // Check for help or version flags
    optmap::iterator oiter = pargs.opts.find('a');
    if (oiter != pargs.opts.end()) {
	print_usage(pargs.exec);
	MPI_Finalize();
	return 0;
    }
    oiter = pargs.opts.find('v');
    if (oiter != pargs.opts.end()) {
	std::cerr << "Percy++ pirserver version " << VERSION << std::endl;
	std::cerr << AUTHOR << std::endl;
	MPI_Finalize();
	return 0;
    }

    // Initialize paramaters and datastore
    PercyMPIServerParams * params = init_params_mpi(pargs);
    if (params == NULL) {
	MPI_Finalize();
	return -1;
    }

    // Split based on role
    int ret = 0;
    switch (params->role()) {
    case MPI_MASTER: {
	PercyServer * server = new PercyMasterServer();
	if (server == NULL) {
	    fprintf(stderr, "[MASTER] The server was not initialized\n");
	    delete params;
	    MPI_Finalize();
	    return -1;
	}

	if (params->has_test_client()) {
	    MPIStreamBuf sb(params->num_nodes() + 1);
	    std::iostream sio(&sb);
	    handle_requests(server, sio, sio, *params);
	} else {
	    char * infile = NULL;
	    char * outfile = NULL;
	    // Check for redirection
	    oiter = pargs.opts.find('F');
	    if (oiter != pargs.opts.end()) {
		infile = oiter->second;
	    }
	    oiter = pargs.opts.find('G');
	    if (oiter != pargs.opts.end()) {
		outfile = oiter->second;
	    }

	    if (infile != NULL && outfile != NULL) {
		fprintf(stderr, "[MASTER] Reading input from '%s'\n", infile);
		fprintf(stderr, "[MASTER] Sending responses to '%s'\n", outfile);
		std::ifstream ifs((const char *)infile);
		std::ofstream ofs((const char *)outfile);
		std::istream is(ifs.rdbuf());
		std::ostream os(ofs.rdbuf());
		handle_requests(server, is, os, *params);
	    } else if (infile != NULL) {
		fprintf(stderr, "[MASTER] Reading input from '%s'\n", infile);
		fprintf(stderr, "[MASTER] Sending responses to stdout\n", outfile);
		std::ifstream ifs(infile);
		std::istream is(ifs.rdbuf());
		handle_requests(server, is, std::cout, *params);
	    } else {
		// Get port
		uint16_t port = 0;
		oiter = pargs.opts.find('p');
		if (oiter != pargs.opts.end()) {
		    port = strtoul(oiter->second, NULL, 10);
		}

		// Create a socket for clients to connect to.
		sockinetbuf sin(sockbuf::sock_stream);
		if (!bind_to_port(sin, port)) {
		    delete params;
		    fprintf(stderr, "Did not successfully bind to a port\n");
		    return -1;
		}

		iosockinet sio(sin.accept());
		handle_requests(server, sio, sio, *params);
	    }
	}

	delete server;
	} break;

    case MPI_WORKER: {
	DataStore * datastore = init_datastore_mpi(*params, pargs);
	if (datastore == NULL) {
	    fprintf(stderr, "[WORKER %d] The datastore was not initialized\n", params->rank());
	    delete params;
	    MPI_Finalize();
	    return -1;
	}

	PercyServer * server = new PercyServer(datastore);
	if (server == NULL) {
	    fprintf(stderr, "[WORKER %d] The server was not initialized\n", params->rank());
	    delete params;
	    delete datastore;
	    MPI_Finalize();
	    return -1;
	}

	MPIStreamBuf sb(0);
	std::iostream sio(&sb);

	PercyServerParams workerparams = params->sub_params();

	handle_requests(server, sio, sio, workerparams);

	delete datastore;
	delete server;
	} break;

    case MPI_TEST_CLIENT: {
	ret = mpi_test_client(*params);
	} break;

    case MPI_NONE: {
	} break;
    }

    // Clean up and return
    delete params;
    MPI_Finalize();

    return ret;
}

#endif

