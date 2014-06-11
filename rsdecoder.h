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

#ifndef __GSDECODER_H__
#define __GSDECODER_H__

#include <vector>
#include <map>
#include <set>
#include <string>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_GF2E.h>
#include "percyresult.h"
#include "FXY.h"
#include "portfolio.h"

// An implementation of the Guruswami-Sudan decoder.  This decoder will
// produce all possible codewords of degree at most t that match the
// given codeword with k non-erasures in more than sqrt(k*t) places.
//
// It runs in polynomial time, but in the worst case, it's O(k^12).
// Luckily, the worst case is when t=k-2, which is easily solved by
// brute force.  (You will need to do this brute force outside this
// class.)
//
// This implementation is based on the K\"otter and Roth-Ruckenstein
// algorithms, as described in:  R. J. McEliece. The Guruswami-Sudan
// Decoding Algorithm for Reed-Solomon Codes. IPN Progress Report
// 42-153, May 15, 2003.
// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html

NTL_CLIENT

// Note: In most places in this file, F (which is usually a field ZZ_p
// or GF2E) may also be a ring ZZ_p modulo a product of two large
// primes.

// A struct containing a reconstructed polynomial over the field F
template <class FX>
struct RecoveryPoly {
    RecoveryPoly(vector<nservers_t> G, FX phi) : G(G), phi(phi) {}
    vector<nservers_t> G;
    FX phi;
};

// A struct containing a reconstructed polynomial over the field F
template <class FX>
struct RecoveryPolyMulti {
	RecoveryPolyMulti(vector<nservers_t> G, vector<FX> phis) : G(G), phis(phis) {}
    vector<nservers_t> G;
    vector<FX> phis;
};

// Contains the results (so far) of decoding words
template <class F>
struct DecoderResult {
    DecoderResult (vector<nservers_t> G, map<dbsize_t, F> recovered) : G(G), recovered(recovered) {}
    vector<nservers_t> G;
    map<dbsize_t, F> recovered;
};


// A class for caching the results of C() computations.  Note that you
// need to keep a separate Ccache around for each different field you
// work in.

typedef pair<unsigned int, unsigned int> pairint;
#define Ccache_t map<pairint, F>

template <class F>
static F C(Ccache_t &Ccache, unsigned int n, unsigned int k)
{
    // Do we already know the answer to this one?
    pairint nk;
    nk.first = n;
    nk.second = k;
    typename Ccache_t::const_iterator Ci = Ccache.find(nk);
    if (Ci != Ccache.end()) {
	// Yes; return it
	return Ci->second;
    }
    
    // Calculate in ZZ in case we're working in a low-characteristic
    // field.
    ZZ num, dem;
    num = 1;
    dem = 1;
    unsigned int i;

    for (i=0;i<k;++i) {
	num *= (n-i);
	dem *= (k-i);
    }

    num /= dem;

    // Now convert it into the right field, cache the answer, and return
    // it
    F res;
    conv(res, num);
    Ccache[nk] = res;
    return res;
}


// Comparison functions used to let polynomials be keys for maps
bool cmp (const GF2EX& a, const GF2EX& b);
bool cmp (const ZZ_pX& a, const ZZ_pX& b);


// The main RSDecoder class, genericized to work over different types of
// fields.  It's been tested on ZZ_p and GF2E.  I haven't tried it over
// ZZ_pE, but I don't know why it wouldn't work there too; it should
// only be necessary to specialize DT (see below) appropriately.
template <class F, class vec_F, class FX, class FXY, class mat_F>
class RSDecoder {
    public:
	// Call this constructor for GF2E.  Be sure to call GF2E::init
	// as well, with an appropriate polynomial modulus.
	RSDecoder() {}

	// Call this constructor for ZZ_p.  p2 should be prime, and p1
	// should be either prime or 1.  The modulus used is p1*p2.
	RSDecoder(const ZZ &p1, const ZZ &p2): p1(p1), p2(p2) {}

	// The main method to do the Reed-Soloman Decoding
	bool Recover (dbsize_t bytes_per_word, nservers_t t, nservers_t h,
		const vector<nservers_t> &goodservers, const vector<vector<vec_F> > &values,
		const vec_F &indices, vector<vector<DecoderResult<F> > > &results,
		vector<std::set<dbsize_t> > &decoded, const nqueries_t multi_only = 0);

	// A helper function to append an element of F to a string.
	// It's public and static so that EasyRecover can use it.
	static string append(const string &s, const F &wz,
		unsigned int bytes_per_word);

	// Set phi to the degree-t polynomial which interpolates the
	// (indices,values) pairs indexed by I.  Check which
	// (indices,values) pairs indexed by G agree and disagree with
	// phi, and set numagree, numdisagree, and vecagree accordingly.
	static void test_interpolate(unsigned short t,
		const vec_F &values, const vec_F &indices,
		const vector<unsigned short> &I,
		const vector<unsigned short> &G,
		unsigned short &numagree, unsigned short &numdisagree,
		vector<unsigned short> &vecagree, FX &phi);

#ifdef TEST_RR
    public:
#else
    private:
#endif
	// Return a list of roots for y of the bivariate polynomial
	// P(x,y).  If degreebound >= 0, only return those roots with
	// degree <= degreebound.  This routine handles the case where F
	// is the integers mod p1*p2 as well as fields.  This routine
	// may also return some spurious values, so you need to validate
	// the answers it gives.
	//
	// This is an implementation of the Roth-Ruckenstein algorithm
	// for finding roots of bivariate polynomials, as explained in
	// section IX of:  R. J. McEliece. The Guruswami-Sudan Decoding
	// Algorithm for Reed-Solomon Codes. IPN Progress Report 42-153,
	// May 15, 2003.
	// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
	vector<FX> findroots(const FXY &P, int degreebound);

#if defined(TEST_FINDPOLYS) || defined(TIME_FINDPOLYS)
    public:
#else
    private:
#endif

	// findpolys wrapper function
	vector< RecoveryPoly<FX> > findpolys(int k,
		unsigned int t, const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares, TestType &testType, DPType
        &dpType, int &gord);
	vector< RecoveryPoly<FX> > findpolys(int k,
		unsigned int t, const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares, TestType &testType)
    {
        DPType dptype; int gord;
        return findpolys(k, t, goodservers, indices, shares, testType, dptype,
                gord);
    }
	vector< RecoveryPoly<FX> > findpolys(int k,
		unsigned int t, const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares)
    {
        TestType testtype; DPType dptype; int gord;
        return findpolys(k, t, goodservers, indices, shares, testtype, dptype,
                gord);
    }

    // Portfolio algorithm for single-polynomial decoding
	vector< RecoveryPoly<FX> > findpolys_best(unsigned int n, int ell,
		unsigned int h, const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares, TestType &testType,
        DPType &dpType, int &gord);

    // Wrapper for multi-polynomial decoding
	vector< RecoveryPolyMulti<FX> > findpolys_multi(unsigned int k,
		unsigned int t, const vector<unsigned short> &goodservers,
		const vec_F &indices, const vector<vec_F> &shares, TestType testType = CH_TK1);

	// Use the K\"otter algorithm (below) to construct an
	// appropriate  bivariate polynomial and use the
	// Roth-Ruckenstein algorithm (above) to extract its roots.
	// Filter the results of the above to end up with just the
	// answers we're looking for.
	vector< RecoveryPoly<FX> > findpolys_gs(unsigned int k, unsigned int t,
		const vector<unsigned short> &goodservers,
		const vec_F& indices, const vec_F& shares, TestType testType = KOTTER);

	// Find all polynomials of degree at most k that agree with the
	// given (index,share) pairs at at least t of the n points.
	// This version simply uses brute force and interpolates all
	// subsets of size k+1 of the n points.  Note that in general,
	// this version does *not* run in polynomial time, but for some
	// cases (with n-k small, for example) it is faster than
	// Guruswami-Sudan.
	//
	vector<RecoveryPoly<FX> > findpolys_brute(
		int k_signed, unsigned int t,
		const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares);

	// Find a polynomial of degree ell that agrees with at least h shares.  Uses
	// the Berlekamp-Welch algorithm.  This algorithm is only guaranteed to work
	// when 2h > n + ell.
	vector<RecoveryPoly<FX> > findpolys_bw(
			unsigned int ell, unsigned int h,
			const vector<unsigned short> &goodservers,
			const vec_F &indices, const vec_F &shares);

    // A dynamic programming approach to decoding.  We can assume that some
    // servers are correct, some are wrong or a combination of the two
    vector<RecoveryPoly<FX> > findpolys_dp (unsigned int ell, unsigned int h,
            const vector<unsigned short> &goodservers,
            const vec_F &indices, const vec_F &shares,
            const DPType, unsigned int gord);

#ifdef INCLUDE_GENERAL_MULTI
	// (Needed for findpolys_ch_multi)
	vector<map<unsigned int, FX> > solve_partial_solution (
			unsigned int ell, unsigned int h, map<unsigned int, FX> partial,
			const vector<vector<FX> >& rows, vector<unsigned int> row_order,
			const vector<vector<unsigned int> >& order_of_expts);

    // The general Cohn-Heninger multi-polynomial decoding algorithm as
    // described in "Approximate common divisors via lattices." (2011)
	vector<RecoveryPolyMulti<FX> > findpolys_ch_multi(
			unsigned int ell, unsigned int h, unsigned int m, unsigned int& t,
			unsigned int& k,
			const vector<unsigned short> &goodservers,
			const vec_F &indices, const vector<vec_F> &shares);
#endif

    // The linear version of the Cohn-Heninger multi-polynomial decoding
    // algorithm as described in "Approximate common divisors via lattices."
    // (2011).  This is the general algorithm with t=k=1.
	vector<RecoveryPolyMulti<FX> > findpolys_ch_tk1(
			unsigned int ell, unsigned int h, unsigned int m,
			const vector<unsigned short> &goodservers,
			const vec_F &indices, const vector<vec_F> &shares);

	// Solve the linear system of equations Ax=b
	vector<FX> solve_linsys_FX (vector<vector<FX> >& A, vector<FX>& b,
            unsigned int modulus = 0);

    private:
	// A structure for holding derived types.  It needs to be
	// manually specialized for every kind of base field (ZZ_p,
	// GF2E, etc.) we may want to work in.
	struct DT {};

	// Construct a bivariate polynomial P(x,y) such that, for any
	// polynomial f(x) of degree at most v that agrees with at least
	// t of the given points, (y-f(x)) is a factor of P(x,y).  This
	// version is K"otter's Interpolation Algorithm, as described in
	// Section VII of R. J.  McEliece. The Guruswami-Sudan Decoding
	// Algorithm for Reed-Solomon Codes. IPN Progress Report 42-153,
	// May 15, 2003.
	// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
	FXY interpolate_kotter(unsigned int v, unsigned int t,
		const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares);

private:
	// Construct a bivariate polynomial Q(x,y) such that, for any
	// polynomial g(x) of degree at most ell that agrees with at least
	// h of the given points, (y-g(x)) is a factor of Q(x,y).  This
	// version is Cohn and Heninger's improvement on the Guruswami-Sudan
	// Decoding Algorithm as described in sections 1.2, 2.2 and 4 of
	// Cohn, H. and Heninger, N.  Ideal forms of Coppersmith's theorem and
	// Guruswami-Sudan list decoding.  August 6, 2010.
	FXY interpolate_cohn_heninger(unsigned int v, unsigned int h,
		const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares);

	// Perform lattice basis reduction on a lattice with m basis vectors in F[z].
	// The lattice is represented by a vector 'lattice' in (F[z][x])^m such that
	// the coefficient vectors of each element in 'lattice' is a basis for the
	// lattice.  The result is a polynomial in F[z][x] whose coefficient vector
	// is the smallest vector in the lattice.  This implementation uses the
	// algorithm found in Mulders, T. and Storjohann, A.  On Lattice Reduction for
	// Polynomial Matrices.  December 2000.
	vector<vector<FX> > reduce_lattice_MS(vector<vector<FX> > lattice,
			// The auxiliary matrix must have the same number of rows as lattice
			// If not, no operations will be performed on aux.
			vector<vector<FX> >& aux,
			const int reduceTo = -1, const unsigned int reduceNumber = 1);

	FXY reduce_lattice_opt(vector<FXY> lattice, const unsigned int n,
			const unsigned int ell, const unsigned int k);

	// Evaluate the (r,s)th Hasse mixed partial derivative of g at
	// the point (alpha, beta), which is:
	//   \sum_{i,j} C(i,r) C(j,s) a_{i,j} alpha^{i-r} beta^{j-s}
	// = \sum_j C(j,s) beta^{j-s} \sum_i C(i,r) a_{i,j} alpha^{i-r}
	// where g(x,y) = \sum_{i,j} a_{i,j} x^i y^j
	// See page 14 of R. J. McEliece. The Guruswami-Sudan Decoding
	// Algorithm for Reed-Solomon Codes. IPN Progress Report 42-153,
	// May 15, 2003.
	// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
	F evalhasse(const FXY &g, unsigned int r, unsigned int s,
		F alpha, F beta, Ccache_t &Ccache);

	// Depth-first search on the tree of coefficients; used in
	// Roth-Ruckenstein.
	void dfs(vector<FX> &res,
		int u, 
		vector<int> &pi, 
		vector<int> &Deg,
		vector<F> &Coeff, 
		vector<FXY> &Q, 
		int &t, 
		int degreebound);

	// Return a list of roots for y of the bivariate polynomial
	// P(x,y).  If degreebound >= 0, only return those roots with
	// degree <= degreebound.  This routine only works over fields F
	// (not rings).
	vector<FX> rr_findroots(const FXY &P, int degreebound);

	// In the case of ZZ_p, store the factors of the modulus p1*p2.
	// If it's a field, p1 should be 1 and p2 should be prime.
	ZZ p1, p2;
};

// The two types of RSDecoders I've tested.  These are the types you
// would use to create a RSDecoder; i.e.:
//
// RSDecoder_ZZ_p decoder(p1, p2);
//
// or
//
// RSDecoder_GF2E decoder();
//
// Remember to also call ZZ_p::init(p1*p2) or GF2E::init(modpoly) as
// appropriate before creating the RSDecoder.
typedef RSDecoder<ZZ_p, vec_ZZ_p, ZZ_pX, ZZ_pXY, mat_ZZ_p> RSDecoder_ZZ_p;
typedef RSDecoder<GF2E, vec_GF2E, GF2EX, GF2EXY, mat_GF2E> RSDecoder_GF2E;

// Comparison functions used to let polynomials be keys for maps
bool cmp (const GF2EX& a, const GF2EX& b);
bool cmp (const ZZ_pX& a, const ZZ_pX& b);

#include "recover.h"

#endif
