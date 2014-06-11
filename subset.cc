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

#include <NTL/ZZ.h>  // For RandomBnd()
#include "subset.h"

NTL_CLIENT

// Operations on subsets of the set of servers.  These subsets are
// stored as sorted vectors of nservers_t.

// Append k random elements of G[offset..end] to I
void random_subset(const vector<nservers_t> &G,
	vector<nservers_t> &I, nservers_t k, nservers_t offset)
{
    if (k == 0) return;

    // How many elements are there to choose from?
    nservers_t n = G.size() - offset;

    if (n == 0) return;

    // With probability k/n, choose G[offset]
    if (RandomBnd(n) < k) {
	I.push_back(G[offset]);
	random_subset(G, I, k-1, offset+1);
    } else {
	random_subset(G, I, k, offset+1);
    }
}

// Compute the intersection of the two *strictly sorted* vectors v1 and v2
vector<nservers_t> intersect(const vector<nservers_t> &v1,
	const vector<nservers_t> &v2)
{
    vector<nservers_t> intersection;
    vector<nservers_t>::const_iterator v1p, v2p;
    v1p = v1.begin();
    v2p = v2.begin();
    while (v1p != v1.end() && v2p != v2.end()) {
	if (*v1p == *v2p) {
	    // This is an element of the intersection
	    intersection.push_back(*v1p);
	    ++v1p;
	    ++v2p;
	} else if (*v1p > *v2p) {
	    ++v2p;
	} else {  // *v1p < *v2p
	    ++v1p;
	}
    }
    return intersection;
}
