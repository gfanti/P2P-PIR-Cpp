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
#include "cmdtools.h"

bool parse_long_opts (int argc, char ** argv, const char * shortopts,
	const struct option * longopts, ParsedArgs& retargs)
{
    // Initialize retargs
    retargs.exec = argv[0];
    retargs.opts = optmap();
    retargs.nonoptc = 0;
    retargs.nonoptv = argvec();

    // Go throught the options and record in the map
#ifdef VERBOSE_PARSE_OPTS
    fprintf(stderr, "Options:\n");
#endif
    int option_index = -1;
    opterr = 0;
    int opt = getopt_long(argc, argv, shortopts, longopts, &option_index);
    while (opt != -1) {
	// Check if we have already seen this option
	if (retargs.opts.find(opt) != retargs.opts.end()) {
	    fprintf(stderr, "Multiple instances of option '%c'\n", opt);
	    return false;
	}

	// Check for invalid option
	if (opt == '?') {
	    fprintf(stderr, "Invalid option '%s'\n", argv[optind-1]);
	    return false;
	}

#ifdef VERBOSE_PARSE_OPTS
	if (option_index >= 0) {
	    fprintf(stderr, "    --%s", longopts[option_index].name);
	} else {
	    fprintf(stderr, "    -%c", opt);
	}
	if (optarg != NULL) {
	    fprintf(stderr, " = %s", optarg);
	}
	fprintf(stderr, "\n");
#endif

	// Add to resulting map
	retargs.opts[opt] = optarg;

	// Next opt
	option_index = -1;
	opt = getopt_long(argc, argv, shortopts, longopts, &option_index);
    }

#ifdef VERBOSE_PARSE_OPTS
    fprintf(stderr, "Non-option Arguments:\n");
#endif
    while (optind < argc) {
#ifdef VERBOSE_PARSE_OPTS
	fprintf(stderr, "    %s\n", argv[optind]);
#endif
	retargs.nonoptv.push_back(argv[optind++]);
	++retargs.nonoptc;
    }
    std::cerr << "Args: " << retargs.nonoptv[0] << " size: " << retargs.nonoptv.size() << std::endl;

    return true;
}

