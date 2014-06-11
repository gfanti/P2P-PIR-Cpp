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

#include <vector>
#include "subset.h"
#include "subset_iter.h"
#include "recover.h"

NTL_CLIENT

bool EasyRecover(dbsize_t bytes_per_word, nservers_t t,
	nservers_t h, vector<DecoderResult<GF2E> > &results,
	dbsize_t word_number, const vec_GF2E &values_vec, 
	const vec_GF2E &indices_vec)
{
    //vector<DecoderResult<GF2E> > newresults;
    vector<DecoderResult<GF2E> >::const_iterator Riter;
    vector<GF2E> newF;
    vector<vector<nservers_t> > newG;

#ifdef VERBOSE_RECOVER
    std::cerr << "EasyRecover:\n";
#endif
    for (Riter = results.begin(); Riter != results.end(); ++Riter) {
#ifdef VERBOSE_RECOVER
    std::cerr << "  iteration:\n";
#endif
        // Pick a random subset I of G, of size t+1
        vector<nservers_t> I;
        vector<nservers_t>::const_iterator Iiter;
        random_subset(Riter->G, I, t+1);

        // Use Lagrange interpolation to find the unique polynomial phi
        // of degree t which matches the points indexed by I
        vec_GF2E I_indices, I_values;
        I_indices.SetLength(t+1);
        I_values.SetLength(t+1);
        nservers_t i = 0;
        for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
            I_indices[i] = indices_vec[*Iiter];
            I_values[i] = values_vec[*Iiter];
        }
        GF2EX phi;
        phi = interpolate(I_indices, I_values);

        // Count the number of points in G that agree, and that
        // disagree, with phi
        nservers_t numagree = 0;
        nservers_t numdisagree = 0;
        newG.push_back(vector<nservers_t>());
        vector<nservers_t>::const_iterator Giter;
        for (Giter = Riter->G.begin(); Giter != Riter->G.end(); ++Giter) {
            if (eval(phi, indices_vec[*Giter]) == values_vec[*Giter]) {
                ++numagree;
                newG.back().push_back(*Giter);
            } else {
                ++numdisagree;
            }
        }

        // If at least h agreed, and less than h-t disagreed, then phi
        // can be the *only* polynomial of degree t matching at least
        // h points.
        if (numagree >= h && numdisagree < h-t) {
            newF.push_back(eval(phi, GF2E::zero()));
#ifdef VERBOSE_RECOVER
	    std::cerr << "        " << phi << "\n";
#endif
        } else {
            // This either isn't the right polynomial, or there may be
            // more than one.  Abort.
            return false;
        }
    }

    for (unsigned int i = 0; i < newF.size(); ++i) {
        results[i].G = newG[i];
        results[i].recovered[word_number] = newF[i];
    }

    return true;
}


bool EasyRecover(dbsize_t bytes_per_word, nservers_t t,
	nservers_t h, vector<DecoderResult<ZZ_p> > &results, 
	dbsize_t word_number, const vec_ZZ_p &values, 
	const vec_ZZ_p &indices)
{
    vector<ZZ_p> newF;
    vector<vector<nservers_t> > newG;
    vector<DecoderResult<ZZ_p> >::const_iterator Riter;
    for (Riter = results.begin(); Riter != results.end();
            ++Riter) {
        // Pick a random subset I of G, of size t+1
        vector<nservers_t> I;
        random_subset(Riter->G, I, t+1);

        nservers_t numagree, numdisagree;
        newG.push_back(vector<nservers_t>());

        ZZ_pX phi;
        RSDecoder_ZZ_p::test_interpolate(t, values, indices, I, Riter->G,
                numagree, numdisagree, newG.back(), phi);

        // If at least h agreed, and less than h-t disagreed, then phi
        // can be the *only* polynomial of degree t matching at least
        // h points.
        if (numagree >= h && numdisagree < h-t) {
            // Find the secret determined by phi
            ZZ_p wz;
            eval(wz, phi, ZZ_p::zero());
            newF.push_back(wz);
        } else {
            // This either isn't the right polynomial, or there may be
            // more than one.  Abort, and we'll use HardRecover.
            return false;
        }
    }

    for (unsigned int i = 0; i < newF.size(); ++i) {
        results[i].G = newG[i];
        results[i].recovered[word_number] = newF[i];
    }

    return true;
}

