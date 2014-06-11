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

#ifndef __SUBSET_ITER_H__
#define __SUBSET_ITER_H__

#include <vector>
#include "percytypes.h"

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

class subset_iterator {
private:
    std::vector<nservers_t> source;
    std::vector<nservers_t> indices;
    std::vector<nservers_t> output;
    nservers_t subset_size;
    bool finished;
public:
    subset_iterator(const std::vector<nservers_t>& source,
	    nservers_t subset_size);
    subset_iterator& operator++();
    const std::vector<nservers_t>& operator*() const;
    bool atend() const;
};

#endif
