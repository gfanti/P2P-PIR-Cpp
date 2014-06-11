//  Percy++ Copyright 2007,2012 Ian Goldberg <iang@cs.uwaterloo.ca>,
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

#include "subset_iter.h"

// Produce all subsets of a given vector<nservers_t> of a given size
//
// Usage:
//    vector<nservers_t> G;
//    ...
//    subset_iterator iter(G, k);
//    while(!iter.atend()) {
//	do_something_with(*iter);
//	++iter;
//    }

subset_iterator::subset_iterator(const std::vector<nservers_t>& source,
	nservers_t subset_size)
{
    this->source = source;
    int slen = source.size();
    if (subset_size > slen) subset_size = slen;
    this->subset_size = subset_size;
    for(int i=0; i<subset_size; ++i) {
	this->indices.push_back(i);
	this->output.push_back(source[i]);
    }
    this->finished = false;
}

subset_iterator& subset_iterator::operator++()
{
    if (!finished) {
	// Find the rightmost index not at its maximal position
	int slen = source.size();
	int pos;
	for (pos=subset_size-1; pos>=0; --pos) {
	    if (indices[pos] != slen-(subset_size-pos)) break;
	}
	if (pos >= 0) {
	    // Increment that index, and set all subsequent elements
	    // to subsequent values.  For example, if slen == 10,
	    // subset_size == 5, and indices == [ 2 4 5 8 9 ], "pos"
	    // here will be 2, pointing to the 5, which we will
	    // increment to 6, and change the subsequent elements to 7
	    // and 8, resulting in [ 2 4 6 7 8 ].
	    nservers_t newval = indices[pos] + 1;
	    for (; pos<subset_size; ++pos) {
		indices[pos] = newval;
		output[pos] = source[newval];
		++newval;
	    }
	} else {
	    finished = true;
	}
    }
    return *this;
}

const std::vector<nservers_t>& subset_iterator::operator*() const
{
    return output;
}

bool subset_iterator::atend() const
{
    return finished;
}

#ifdef TEST_SUBSET_ITER
#include <iostream>

static void show(const std::vector<nservers_t>& v)
{
    std::cout << "[ ";
    for (nservers_t i=0; i<v.size(); ++i) {
	std::cout << v[i] << " ";
    }
    std::cout << "]\n";
}

int main(int argc, char **argv)
{
    std::vector<nservers_t> G;
    G.push_back(1);
    G.push_back(4);
    G.push_back(9);
    G.push_back(16);
    G.push_back(25);
    G.push_back(36);
    G.push_back(49);
    G.push_back(64);
    G.push_back(81);
    G.push_back(100);

    int k = argc > 1 ? atoi(argv[1]) : 5;

    subset_iterator iter(G, k);
    while(!iter.atend()) {
	show(*iter);
	++iter;
    }
}
#endif
