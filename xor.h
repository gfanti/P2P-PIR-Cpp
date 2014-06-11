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

typedef unsigned int v4ui __attribute__ ((vector_size (16)));

inline void XOR_equal(unsigned char *dst, const unsigned char *src, size_t len) {
    v4ui *dl = (v4ui *)dst;
    const v4ui *sl = (const v4ui *)src;
    size_t ll = len / (8*sizeof(v4ui));
    size_t lr = len % (8*sizeof(v4ui));
    for(size_t i=0; i<ll; i++) {
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
    }
    if (lr) {
        dst = (unsigned char *)dl;
        src = (const unsigned char *)sl;
        for(size_t i=0; i<lr; i++) {
            *(dst++) ^= *(src++);
        }
    }
}

