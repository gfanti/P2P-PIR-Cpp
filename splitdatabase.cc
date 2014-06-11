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
#include <fstream>
#include <sstream>
#include <string.h>
#include <NTL/ZZ_pX.h>
#include "percyio.h"
#include "config.h"

NTL_CLIENT

// Split a given database into l pieces with tau-independence: no
// coalition of up to tau of the pieces has any information about the
// contents of the original database.  You must supply a modulus to use
// for the splitting, and that same modulus must be used during the PIR
// protocol to query the pieces.
//
// If the original database name is "dbname", the pieces will be called
// "dbname.1", "dbname.2", etc.
int main(int argc, char **argv)
{
    int c;
    for (c = 1; c < argc; c++) {
	if (strcmp("--version", argv[c]) == 0) {
	    std::cerr << "Percy++ splitdatabase version " << VERSION
		<< std::endl;
	    std::cerr << AUTHOR << std::endl;
	    exit(0);
	}
    }

    if (argc < 4) {
	std::cerr << "Usage: " << argv[0] << " database tau l [modulus]\n";
	exit(1);
    }

    // Initialize the global ZZ_p context with the given modulus
    stringstream ss (stringstream::in | stringstream::out);
    ZZ modulus;
    if (argc > 4) {
	ss << argv[4];
    } else {
	// Use the standard 1024-bit modulus from pirclient
	ss << "343308946066366926839932845260501528909643718159825813630709694160026342456154871924497152436552679706642965502704642456637620829912957820221098686748075257358288200837461739492534713539606088624083011849535450485951774635526473457667739540374042376629835941950802202870595346459371144019363420985729553740241";
    }
    ss >> modulus;
    ZZ_p::init(modulus);

    streamsize outbytes = NumBytes(modulus);
    streamsize inbytes = outbytes-1;
    if (outbytes < 2) {
	std::cerr << "Error: modulus is too small; must be at least 256.\n";
	exit(1);
    }

    int tau = strtoul(argv[2], NULL, 10);
    int l = strtoul(argv[3], NULL, 10);

    std::cout << "Converting " << inbytes << "-byte inputs to " << outbytes
	<< "-byte outputs.\n";

    // Open the input database
    ifstream infile(argv[1]);
    int flen = strlen(argv[1]);

    // Create l output databases
    ofstream * outfile = new ofstream[l];
    int i;
    char * ofile = new char[flen+13];
    for(i=1;i<=l;++i) {
	sprintf(ofile, "%s.%d", argv[1], i);
	outfile[i-1].open(ofile);
	// Write the tau-independence header for each output database
	outfile[i-1].write("PIRD\x01\x00", 6);
	percy_write_ZZ(outfile[i-1], modulus);
    }
    delete[] ofile;

    unsigned char * inbuf = new unsigned char[inbytes];
    unsigned char * outbuf = new unsigned char[outbytes];

    while(1) {
	// Read the input in chunks of inbytes bytes
	infile.read((char *)inbuf, inbytes);

	if (infile.gcount() == 0) break;

	// Convert the chunk to a ZZ
	ZZ Wz = ZZFromBytes(inbuf, infile.gcount());
	// cout << "reading " << Wz << "\n";

	// Pick a random polynomial of degree tau, and set the constant
	// coefficient to the ZZ read from the input database
	ZZ_pX randpoly = random_ZZ_pX(tau+1);
	SetCoeff(randpoly, 0, to_ZZ_p(Wz));
	// cout << "poly(" << tau << ") = " << randpoly << "\n";

	// For each output file i, output the value of the polynomial
	// evaluated at i.
	for(i=1;i<=l;++i) {
	    ZZ_p value;
	    eval(value, randpoly, to_ZZ_p(i));
	    BytesFromZZ(outbuf, rep(value), outbytes);
	    outfile[i-1].write((char *)outbuf, outbytes);
	    // cout << "writing " << i << "/" << tau << ": " << rep(value) << "\n";
	}
    }

    for(i=0;i<l;++i) {
	outfile[i].close();
    }
    infile.close();

    delete[] outbuf;
    delete[] inbuf;
    delete[] outfile;

    return 0;
}
