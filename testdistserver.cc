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


#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/prctl.h>
#include <sstream>
#include <vector>
#include <string.h>
#include "percytypes.h"
#include "cmdtools.h"

void usage (const char * exec)
{
    fprintf(stderr,
	"\n"
        "Usage: %s [OPTIONS] DBSIZE BLOCKSIZE VSPLIT HSPLIT DB_1 ... DB_M\n"
        "Arguments:\n"
        "   DBSIZE       The size of the database in bytes.\n"
        "   BLOCKSIZE    The size of each block of the database in bytes.\n"
        "                BLOCKSIZE must evenly divide DBSIZE.  We set NUMBLOCKS\n"
        "                to DBSIZE/BLOCKSIZE.\n"
        "   VSPLIT       The database is split evenly into VSPLIT sets of rows.\n"
        "                VSPLIT must evenly divide NUMBLOCKS.\n"
        "   HSPLIT       The database is split evenly into HSPLIT sets of columns.\n"
        "                HSPLIT must evenly divide BLOCKSIZE.\n"
        "   DB_i         The subdatabase at indices (i %% VSPLIT, i / HSPLIT).\n"
        "                Each subdatabase must have size at least DBSIZE/M where\n"
        "                M = VSPLIT * HSPLIT.\n"
        "Options:\n"
        "   -m METHOD    The method to use. (Default: gf28)\n"
        "   -b           Be byzantine.\n"
        "   -t           Use tau-independence.\n"
        "   -S SID       Specify the SID of the master.\n"
        "   -p PORT      Specify the port of the master.\n"
        "   -T THREADS   Specify the number of threads for each worker. (Default: not set)\n"
        "   --verbose    Print the output of the workers.\n"
        "   --gdb        Run the master server under gdb.\n"
        "   --no-kill    Do not kill the workers after running.\n"
	"   --test-cmds  Do not actually run commands.\n"
        "   --help       Display this message.\n"
	"\n"
        "Available methods are {w160, w32, w16, w8, gf28, gf216, chor, spir}.\n"
	"\n", exec);
}

#define NUM_METHODS 8

enum Method {
    METHOD_GF28,
    METHOD_GF216,
    METHOD_W160,
    METHOD_W32,
    METHOD_W16,
    METHOD_W8,
    METHOD_CHOR,
    METHOD_SPIR
};

const char * method_str[NUM_METHODS] = {
    "gf28",
    "gf216",
    "w160",
    "w32",
    "w16",
    "w8",
    "chor",
    "spir"
};

int main (int argc, char ** argv)
{
    // Handle arguments and options
    struct option longopts[] = {
	{"method",	required_argument,  NULL,   'm'},
	{"byzantine",	no_argument,	    NULL,   'b'},
	{"tau",		no_argument,	    NULL,   't'},
	{"sid",		required_argument,  NULL,   'S'},
	{"port",	required_argument,  NULL,   'p'},
	{"num-threads",	required_argument,  NULL,   'T'},
	{"verbose",	no_argument,	    NULL,   'v'},
	{"gdb",		no_argument,	    NULL,   'g'},
	{"no-kill",	no_argument,	    NULL,   'n'},
	{"help",	no_argument,	    NULL,   'h'},
	{"test-cmds",	no_argument,	    NULL,   'c'},
	{NULL,		0,		    NULL,   0},
    };

    const char * shortopts = "m:btS:p:T:";

    ParsedArgs pargs;
    if (!parse_long_opts(argc, argv, shortopts, longopts, pargs)) {
	usage(argv[0]);
	return 1;
    }

    // Print usage
    if (argc == 1 || pargs.opts.find('h') != pargs.opts.end()) {
	usage(pargs.exec);
	return 1;
    }

    // Sanity checks
    struct stat st_pscc, st_ps, st_dscc, st_psm;
    bool exists_pscc = ! stat("pirserver.cc", &st_pscc);
    bool exists_ps = ! stat("pirserver", &st_ps);
    bool exists_dscc = ! stat("distserver.cc", &st_dscc);
    bool exists_psm = ! stat("pirserver_master", &st_psm);

    // Make sure we've run make
    if ((exists_pscc && !exists_ps) || (exists_dscc && !exists_psm)) {
	fprintf(stderr,
	    "Error:\n"
	    "You must build the Percy++ tools before running $EXE.\n"
	    "Please see the README for instructions.\n");
	return 1;
    }

    // Make sure we can see the percy tools
    if (!exists_ps || !exists_psm) {
	fprintf(stderr,
	    "Error:\n"
	    "Cannot find local instances of Percy++ tools.\n"
	    "Percy++ %s must be invoked in the same directory\n"
	    "as the Percy++ tools, usually percy++-(version number).\n",
	    argv[0]);
	return 1;
    }

    // Make sure we can call the percy tools
    if (!(st_ps.st_mode & S_IXUSR) || !(st_psm.st_mode & S_IXUSR)) {
	fprintf(stderr,
	    "Error:\n"
	    "Percy++ tools are not executable; please check permissions.\n");
	return 1;
    }



    // Mandatory arguments
    if (pargs.nonoptc < 4) {
	usage(argv[0]);
	return 1;
    }

    dbsize_t dbsize = strtoull(pargs.nonoptv[0], NULL, 10);
    dbsize_t blocksize = strtoull(pargs.nonoptv[1], NULL, 10);
    if (dbsize % blocksize != 0) {
	fprintf(stderr, "Error: BLOCKSIZE must divide DBSIZE\n");
	return 1;
    }
    dbsize_t numblocks = dbsize / blocksize;

    dbsize_t vsplit = strtoull(pargs.nonoptv[2], NULL, 10);
    if (numblocks % vsplit != 0) {
	fprintf(stderr, "Error: VSPLIT must divide NUMBLOCKS\n");
	return 1;
    }
    dbsize_t hsplit = strtoull(pargs.nonoptv[3], NULL, 10);
    if (blocksize % hsplit != 0) {
	fprintf(stderr, "Error: HSPLIT must divide BLOCKSIZE\n");
	return 1;
    }
    dbsize_t numworkers = vsplit * hsplit;

    // Get database files
    if (pargs.nonoptc < numworkers + 4) {
	fprintf(stderr, "Error: Not enough subdatabases specified\n");
	return 1;
    }
    std::vector<char*> databases;
    for (dbsize_t i = 0; i < numworkers; ++i) {
	databases.push_back(pargs.nonoptv[i + 4]);
    }

    Method method = METHOD_GF28;
    bool byz = false;
    bool tauind = false;
    nservers_t sid = 1;
    uint16_t port = 31337;
    bool threaded = false;
    dbsize_t numthreads = 0;
    bool verbose = false;
    bool gdb = false;
    bool nokill = false;
    bool testcmd = false;

    // Get options
    optmap::iterator oiter;
    bool flag = false;
    for (oiter = pargs.opts.begin(); oiter != pargs.opts.end(); ++oiter) {
	const char opt = oiter->first;
	const char * optarg = oiter->second;
	switch (opt) {
	case 'm':
	    flag = true;
	    for (int i = 0; i < NUM_METHODS; ++i) {
		if (!strcmp(optarg, method_str[i])) {
		    method = (Method)i;
		    flag = false;
		    break;
		}
	    }
	    if (flag) {
		fprintf(stderr, "Error: invalid method '%s'\n", optarg);
		return 1;
	    }
	    break;
	case 'b':
	    byz = true;
	    break;
	case 't':
	    tauind = true;
	    break;
	case 'S':
	    sid = strtoul(optarg, NULL, 10);
	    break;
	case 'p':
	    port = strtoul(optarg, NULL, 10);
	    break;
	case 'T':
	    threaded = true;
	    numthreads = strtoull(optarg, NULL, 10);
	    break;
	case 'v':
	    verbose = true;
	    break;
	case 'g':
	    gdb = true;
	    break;
	case 'n':
	    nokill = true;
	    break;
	case 'c':
	    testcmd = true;
	    break;
	case 'h':
	    usage(pargs.exec);
	    return 0;
	    break;
	default:
	    usage(pargs.exec);
	    return 1;
	    break;
	}
    }

    // Prepare worker commands
    dbsize_t sub_dbsize = dbsize / numworkers;
    dbsize_t sub_blocksize = blocksize / hsplit;
    std::stringstream worker_cmd;
    worker_cmd << "./pirserver"
	    << " -n " << sub_dbsize
	    << " -b " << sub_blocksize;

    // Prepare master commands
    std::stringstream master_cmd;
    master_cmd << "./pirserver_master"
	    << " -n " << dbsize
	    << " -b " << blocksize;

    // Manage methods
    dbbits_t wordsize = 8;
    switch (method) {
    case METHOD_GF28:
	worker_cmd << " -w 8 -m g";
	master_cmd << " -w 8 -m g";
	break;
    case METHOD_GF216:
	worker_cmd << " -w 16 -m s";
	master_cmd << " -w 16 -m s";
	break;
    case METHOD_W160:
	worker_cmd << " -w 160 -m z";
	master_cmd << " -w 160 -m z";
	break;
    case METHOD_W32:
	worker_cmd << " -w 32 -m z";
	master_cmd << " -w 32 -m z";
	break;
    case METHOD_W16:
	worker_cmd << " -w 16 -m z";
	master_cmd << " -w 16 -m z";
	break;
    case METHOD_W8:
	worker_cmd << " -w 8 -m z";
	master_cmd << " -w 8 -m z";
	break;
    case METHOD_CHOR:
	worker_cmd << " -w 1 -m c";
	master_cmd << " -w 1 -m c";
	break;
    case METHOD_SPIR:
	worker_cmd << " -w 160 -m z";
	master_cmd << " -w 160 -m z -s polycommit.params";
	break;
    }

    // Other options
    if (byz) {
	worker_cmd << " -z";
	master_cmd << " -z";
    }
    if (tauind) {
	worker_cmd << " -t";
	master_cmd << " -t";
    }
    master_cmd << " -S " << sid << " -p " << port;
    if (threaded) {
	worker_cmd << " -T " << numthreads;
    }

    master_cmd << " " << vsplit << " " << hsplit << " '";

    // Start the workers
    nservers_t sidno = 1;
    uint16_t portno = 31337;
    for (dbsize_t i = 0; i < numworkers; ++i) {
	if (sid == sidno) ++sidno;
	if (port == portno) ++portno;
	if (i != 0) master_cmd << " ";
	master_cmd << sidno << ":localhost:" << portno;
	pid_t pid = fork();
	if (pid < 0) {
	    fprintf(stderr, "Fork %lu failed\n", i);
	    return 1;
	} else if (pid == 0) {
	    // Child i
	    if (!nokill && !testcmd) {
		prctl(PR_SET_PDEATHSIG, SIGKILL, 0, 0, 0);
	    }
	    std::stringstream curr_worker_cmd;
	    curr_worker_cmd << worker_cmd.str()
		    << " -S " << sidno
		    << " -p " << portno
		    << " " << databases[i];
	    printf("Starting worker %lu ...\n", i);
	    if (verbose || testcmd) {
		fprintf(stderr, "%s\n", curr_worker_cmd.str().c_str());
	    } else {
		curr_worker_cmd << " > /dev/null 2>&1";
	    }
	    if (!testcmd) {
		system(curr_worker_cmd.str().c_str());
	    }
	    return 0;
	}
	++sidno;
	++portno;
    }
    master_cmd << "'";

    // Parent
    // exec the master
    printf("Starting master ...\n");
    if (verbose || testcmd) {
	fprintf(stderr, "%s\n", master_cmd.str().c_str());
    }
    if (!testcmd) {
	system(master_cmd.str().c_str());
    }

    return 0;
}

