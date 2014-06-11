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

#ifndef __GSDECODER_IMPL_H__
#define __GSDECODER_IMPL_H__

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_GF2E.h>
#include "subset.h"
#include "subset_iter.h"
#include "percyresult.h"
#include "FXY.h"

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

// This file contains the (templated) implementations of the class
// methods.

// Specializations for the derived types.  Only vec_FX and
// vec_pair_FX_long are in here, because they're only used once.  The
// rest of the derived types (vec_F, FX, FXY) are used many times, so
// it's more consice to make them explicit instead of having to write
// DT::vec_F every time.
template<>
struct RSDecoder_ZZ_p::DT {
    typedef vec_ZZ_pX vec_FX;
    typedef vec_pair_ZZ_pX_long vec_pair_FX_long;
};

template<>
struct RSDecoder_GF2E::DT {
    typedef vec_GF2EX vec_FX;
    typedef vec_pair_GF2EX_long vec_pair_FX_long;
};

// Set phi to the degree-t polynomial which interpolates the
// (indices,values) pairs indexed by I.  Check which
// (indices,values) pairs indexed by G agree and disagree with
// phi, and set numagree, numdisagree, and vecagree accordingly.
template<class F, class vec_F, class FX, class FXY, class mat_F>
void RSDecoder<F,vec_F,FX,FXY,mat_F>::test_interpolate(unsigned short t,
	const vec_F &values, const vec_F &indices,
	const vector<unsigned short> &I, const vector<unsigned short> &G,
	unsigned short &numagree, unsigned short &numdisagree,
	vector<unsigned short> &vecagree, FX &phi)
{
    numagree = numdisagree = 0;

    // Use Lagrange interpolation to find the unique polynomial phi
    // of degree t which matches the points indexed by I
    vec_F I_indices, I_values;
    I_indices.SetLength(t+1);
    I_values.SetLength(t+1);
    vector<unsigned short>::const_iterator Iiter;
    unsigned short i = 0;
    for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
	I_indices[i] = indices[*Iiter];
	I_values[i] = values[*Iiter];
    }
    interpolate(phi, I_indices, I_values);

    // Count the number of points in G that agree, and that
    // disagree, with phi
    vector<unsigned short>::const_iterator Giter;
    for (Giter = G.begin(); Giter != G.end(); ++Giter) {
	F phival;
	eval(phival, phi, indices[*Giter]);
	if (phival == values[*Giter]) {
	    ++numagree;
	    vecagree.push_back(*Giter);
	} else {
	    ++numdisagree;
	}
    }
}

extern uint64_t hasseop;

// Evaluate the (r,s)th Hasse mixed partial derivative of g at the point
// (alpha, beta), which is:
//   \sum_{i,j} C(i,r) C(j,s) a_{i,j} alpha^{i-r} beta^{j-s}
// = \sum_j C(j,s) beta^{j-s} \sum_i C(i,r) a_{i,j} alpha^{i-r}
// where g(x,y) = \sum_{i,j} a_{i,j} x^i y^j
// See page 14 of R. J. McEliece. The Guruswami-Sudan Decoding
// Algorithm for Reed-Solomon Codes. IPN Progress Report 42-153, May 15,
// 2003.  http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
template<class F, class vec_F, class FX, class FXY, class mat_F>
F RSDecoder<F,vec_F,FX,FXY,mat_F>::evalhasse(const FXY &g,
	unsigned int r, unsigned int s,
	F alpha, F beta, Ccache_t &Ccache)
{
    F res;
    res = 0;
    int ydeg = deg(g);
    for (int j = ydeg; j >= (int)s; --j) {
	const FX &gj = coeff(g,j);
	int xdeg = deg(gj);
	// Use Horner's method to evaluate the inner sum (which is the
	// coefficient of beta^{j-s})
	F resx;
	resx = 0;
	for (int i = xdeg; i >= (int) r; --i) {
	    resx *= alpha;
	    resx += C(Ccache, i,r) * coeff(gj, i);
	    ++hasseop;
	}
	// Use Horner's method to accumulate the results into the outer
	// sum
	res *= beta;
	res += C(Ccache, j,s) * resx;
	++hasseop;
    }

    return res;
}

extern uint64_t kotter_usec;

// Construct a bivariate polynomial P(x,y) such that, for any polynomial
// f(x) of degree at most v that agrees with at least t of the given
// points, (y-f(x)) is a factor of P(x,y).  This version is K"otter's
// Interpolation Algorithm, as described in Section VII of R. J.
// McEliece. The Guruswami-Sudan Decoding Algorithm for Reed-Solomon
// Codes. IPN Progress Report 42-153, May 15, 2003.
// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
template<class F, class vec_F, class FX, class FXY, class mat_F>
FXY RSDecoder<F,vec_F,FX,FXY,mat_F>::interpolate_kotter(
	unsigned int v, unsigned int t,
	const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
    struct timeval st, et;
    gettimeofday(&st, NULL);

    unsigned int n = goodservers.size();

    // Compute the m and L parameters
    unsigned int m = 1 + (unsigned int)(floor( v*n / (t*t-v*n)));
    unsigned int L = (m*t - 1)/v;

    if (getenv("PIRC_L")) {
    	L = atoi(getenv("PIRC_L"));
    }
    if (getenv("PIRC_m")) {
    	m = atoi(getenv("PIRC_m"));
    }

    std::cerr << "Constructing (1," << v << ")-degree " << L*v << " polynomial...\n";
    std::cerr << "Estimated work: " << n * m * (m+1) / 2 * (L+1) << "\n";
    std::cerr << "L = " << L << "\n";
    std::cerr << "m = " << m << "\n";
    std::cerr << "Min matches: " << t << "\n";
    std::cerr << "Max degree: " << v << "\n";
#if 0
    double Km = v * n * (m+1);
    Km /= (double) m;
    Km = floor(sqrt(Km));
    std::cerr << "Km ~= " << Km << "\n";
    unsigned int C = n * m * (m+1) / 2;
    for (int K=0;;++K) {
	cerr << (K*(K+v)+(K%v)*(v-(K%v)))/(2*v) << " " << C << "\n";
	if ( ((K*(K+v)+(K%v)*(v-(K%v)))/(2*v)) > C ) {
	    std::cerr << "Km: " << (K-1)/m + 1 << "\n";
	    break;
	}
    }
#endif

    // Initialize the g vector
    typedef pair<FXY, unsigned int> polydeg;
    polydeg * g = new polydeg[L+1];
    for (unsigned int j = 0; j <= L; ++j) {
	SetCoeff(g[j].first, j);
	g[j].second = j * v;
    }

    Ccache_t Ccache;

    F * Delta = new F[L+1];
    for (unsigned int i = 0; i < n; ++i) {
	F alpha = indices[goodservers[i]];
	F beta = shares[goodservers[i]];
	for (unsigned int r = 0; r < m; ++r) {
	    for (unsigned int s = 0; s < m - r; ++s) {
		int seennonzero = 0;
		unsigned int seendeg = 0, jstar = 0;
		for (unsigned int j = 0; j <= L; ++j) {
		    Delta[j] = evalhasse(g[j].first, r, s, alpha, beta,
			    Ccache);
		    // cerr << i << " " << r << " " << s << " " << j;
		    if (Delta[j] != 0) {
			// cerr << " nonzero";
			seennonzero = 1;
			seendeg = g[j].second;
			jstar = j;
		    }
		    // cerr << "\n";
		}
		if (seennonzero) {
		    for (unsigned int j = 0; j <= L; ++j) {
			if (Delta[j] != 0 && g[j].second <= seendeg) {
			    seendeg = g[j].second;
			    jstar = j;
			}
		    }
		    F Deltajstar = Delta[jstar];
		    FXY f = g[jstar].first;
		    // cerr << "Deltajstar = " << Deltajstar << "\n";
		    // cerr << "f = " << f << "\n";
		    for (unsigned int j = 0; j <= L; ++j) {
			if (Delta[j] != 0) {
			    if (j != jstar) {
				// cerr << "g["<<j<<"] = " << Deltajstar << " * " << g[j].first << " - " << Delta[j] << " * " << f << " = ";

				g[j].first = Deltajstar * g[j].first -
				    Delta[j] * f;
				// cerr << g[j].first << "\n";
			    } else {
				FX xminusalpha;
				SetCoeff(xminusalpha, 1, 1);
				SetCoeff(xminusalpha, 0, -alpha);
				// cerr << "g["<<j<<"] = " << Deltajstar << " * " << xminusalpha << " * " << f << " = ";
				g[j].first = Deltajstar * xminusalpha * f;
				// cerr << g[j].first << "\n";
				g[j].second += 1;
			    }
			}
			// cerr << "Now -> " << evalhasse(g[j].first, r, s, alpha, beta, Ccache) << "\n";
		    }
		}
	    }
	}
    }
    delete[] Delta;

    // Return the poly of least weighted degree from g
    unsigned int minweight = g[0].second;
    unsigned int minindex = 0;
    for (unsigned int i=1; i<=L; ++i) {
	if (g[i].second <= minweight) {
	    minweight = g[i].second;
	    minindex = i;
	}
    }

    gettimeofday(&et, NULL);
    kotter_usec = ((uint64_t)(et.tv_sec - st.tv_sec))*1000000
	+ (et.tv_usec - st.tv_usec);

    delete[] g;

    return g[minindex].first;
}

// Construct a bivariate polynomial Q(x,y) such that, for any
// polynomial g(x) of degree at most ell that agrees with at least
// h of the given points, (y-g(x)) is a factor of Q(x,y).  This
// version is Cohn and Heninger's improvement on the Guruswami-Sudan
// Decoding Algorithm as described in sections 1.2, 2.2 and 4 of
// Cohn, H. and Heninger, N.  Ideal forms of Coppersmith's theorem and
// Guruswami-Sudan list decoding.  August 6, 2010.
template<class F, class vec_F, class FX, class FXY, class mat_F>
FXY RSDecoder<F,vec_F,FX,FXY,mat_F>::interpolate_cohn_heninger(
	unsigned int ell, unsigned int h,
	const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
#ifdef VERBOSE_CH
	std::cerr << "STARTING COHN-HENINGER...\n";
#endif

	unsigned int n = goodservers.size();
	unsigned int m = ((h - ell) * n) / (h * h - ell * n) + 1;
	unsigned int k = (h * m) / n;
	unsigned int t = m - k;

	unsigned int i, j;

	FX p;
	SetCoeff(p, 0, 1);
	for (i = 0; i < n; ++i) {
		FX newTerm;
		SetCoeff(newTerm, 0, -(indices[goodservers[i]]));
		SetCoeff(newTerm, 1);
		p *= newTerm;
	}

	FX constCoeff;
	for (i = 0; i < n; ++i) {
		FX zPoly;
		SetCoeff(zPoly, 0, 1);
		for (j = 0; j < n; ++j) {
			if (i == j) {
				continue;
			}
			FX newTerm;
			SetCoeff(newTerm, 0, -(indices[goodservers[j]]));
			SetCoeff(newTerm, 1);
			zPoly *= newTerm;
			zPoly /= (indices[goodservers[i]] - indices[goodservers[j]]);
		}
		constCoeff += (zPoly * shares[goodservers[i]]);
	}
	FXY f;
	SetCoeff(f, 1, 1);
	SetCoeff(f, 0, -constCoeff);

#ifdef VERBOSE_CH
	std::cerr << "ell = " << ell << std::endl;
	std::cerr << "h = " << h << std::endl;
	std::cerr << "n = " << n << std::endl;
	std::cerr << "m = " << m << std::endl;
	std::cerr << "k = " << k << std::endl;
	std::cerr << "t = " << t << std::endl;
	std::cerr << "goodservers =";
	for (i = 0; i < n; ++i) {
		std::cerr << " " << goodservers[i];
	}
	std::cerr << std::endl;
	std::cerr << "indices =";
	for (i = 0; i < n; ++i) {
		std::cerr << " " << indices[goodservers[i]];
	}
	std::cerr << std::endl;
	std::cerr << "shares =";
	for (i = 0; i < n; ++i) {
		std::cerr << " " << shares[goodservers[i]];
	}
	std::cerr << std::endl;
	std::cerr << "p(z) = " << p << std::endl;
	std::cerr << "f(x) = " << f << std::endl;
#endif

	vector<FXY> basis; // Polynomials whose coefficient vectors are a lattice basis

	FX z_toell;
	SetCoeff(z_toell, ell);
	FXY z_toell_x;
	SetCoeff(z_toell_x, 1, z_toell);
	FXY f_of_z_toell_x = f;
	SetCoeff(f_of_z_toell_x, 1, z_toell);

	vector<FXY> powers_of_f;
	FXY one2;
	SetCoeff(one2, 0, 1);
	powers_of_f.push_back(one2);
	powers_of_f.push_back(f_of_z_toell_x);
	for (i = 2; i <= k; ++i) {
		powers_of_f.push_back(powers_of_f.back() * f_of_z_toell_x);
	}

	vector<FX> powers_of_p;
	FX one1;
	SetCoeff(one1, 0, 1);
	powers_of_p.push_back(one1);
	powers_of_p.push_back(p);
	for (i = 2; i <= k; ++i) {
		powers_of_p.push_back(powers_of_p.back() * p);
	}

	for (i = 0; i <= k; ++i) {
		FXY basis_poly = powers_of_f[i] * powers_of_p[k - i];
		basis.push_back(basis_poly);
	}

	char *opt_ch_env = getenv("PIRC_OPT_CH");
	if (opt_ch_env && atoi(opt_ch_env)) {
		// k+1-th basis vector (for if ell=n-2)
		if (ell == n-2) {
			FXY last_poly = basis[k] * z_toell_x;
			basis.push_back(last_poly);
		}

	} else {
		FXY current_basis_poly = powers_of_f[k];
		for (j = 1; j < t; ++j) {
			current_basis_poly *= z_toell_x;
			basis.push_back(current_basis_poly);
		}
	}

	vector<vector<FX> > lattice;
	for (i = 0; i < m; ++i) {
		vector<FX> row;
		for (j = 0; j < m; ++j) {
			row.push_back(coeff(basis[i], j));
		}
		lattice.push_back(row);
	}


#ifdef VERBOSE_CH

	std::cerr << "ELEMENT DEGREES:\n";
	for (i = 0; i < m; ++i) {
		std::cerr << "\t";
		for (j = 0; j < m; ++j) {
			std::cerr << deg(lattice[i][j]) << "\t";
		}
		std::cerr << std::endl;
	}
#endif

	vector<vector<FX> > aux (m, vector<FX>());
	vector<vector<FX> > reducedLattice = reduce_lattice_MS(lattice, aux, k*h);
	int mindeg = 0;
	int minindex = -1;
	for (i = 0; i < m; ++i) {
		int maxdeg = -1;
		for (j = 0; j < m; ++j) {
			FX c = reducedLattice[i][j];
			if (deg(c) >= maxdeg) {
				maxdeg = deg(c);
			}
		}
		if (minindex < 0 || maxdeg < mindeg) {
			minindex = i;
			mindeg = maxdeg;
		}
	}
//	FXY Q = reducedBasis[minindex];
	FXY Q;
	for (i = 0; i < m; ++i) {
		SetCoeff(Q, i, reducedLattice[minindex][i]);
	}

#ifdef VERBOSE_CH
	std::cerr << "Q = [" << std::endl;
	for (i = 0; i < m; ++i) {
		FX c = coeff(Q, i);
		std::cerr << "          " << c << " (degree " << deg(c) + 1 << ")"
				<< std::endl;
	}
	std::cerr << "          ]" << std::endl;
#endif

	FXY new_Q;
	for (i = 0; i < m; ++i) {
		SetCoeff(new_Q, i, RightShift(coeff(Q, i), i*ell));
	}

#ifdef VERBOSE_CH
	std::cerr << "new_Q = [" << std::endl;
	for (i = 0; i < m; ++i) {
		FX c = coeff(new_Q, i);
		std::cerr << "              " << c << " (degree " << deg(c) + 1 << ")"
				<< std::endl;
	}
	std::cerr << "         ]" << std::endl;
#endif

	return new_Q;
}

// Perform lattice basis reduction on a lattice with m basis vectors in F[z].
// The lattice is represented by a vector 'lattice' in (F[z][x])^m such that
// the coefficient vectors of each element in 'lattice' is a basis for the
// lattice.   If set, the algorithm stops after finding reduceNumber rows
// of degree less than reduceTo.  This implementation uses the
// algorithm found in Mulders, T. and Storjohann, A.  On Lattice Reduction for
// Polynomial Matrices.  December 2000.
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<vector<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::reduce_lattice_MS(vector<vector<FX> > lattice,
		// The auxiliary matrix must have the same number of rows as lattice
		// If not, no operations will be performed on aux.
		vector<vector<FX> >& aux,
		const int reduceTo, const unsigned int reduceNumber)
{

	// Keep track of number of rows with small enough degree
	map<unsigned int, bool> small_rows;

	// Keeps track of pivot information (row, degree and leading coefficient)
	// for a particular column of a lattice matrix.
	struct PivotInfo {
		bool used;
		unsigned int row;
		unsigned int degree;
		F lc;

		PivotInfo () {
		    used = false;
		    row = 0;
		    degree = 0;
		    lc = 0;
		}
	};

	unsigned int m = lattice.size();

	bool use_aux = aux.size() == m && aux[0].size() > 0;

	// Create the pivots array
	PivotInfo * pivots = new PivotInfo[m];

	for (unsigned int row = 0; row < m; ++row) {
		unsigned int workrow = row;
		// Find pivot
		unsigned int pIndex = 0;
		unsigned int degree = deg(lattice[workrow][0]) + 1;
		for (unsigned int i = 0; i < m; ++i) {
			if ((unsigned int)(deg(lattice[workrow][i]) + 1) >= degree) {
				pIndex = i;
				degree = deg(lattice[workrow][i]) + 1;  // Use +1 to avoid -1 degree
			}
		}
		F lc = LeadCoeff(lattice[workrow][pIndex]);

		// Check if we can stop
		if (reduceTo >= 0 && degree <= (unsigned int)reduceTo) {
			small_rows[workrow] = true;
		}
		if (small_rows.size() >= reduceNumber) {
			delete[] pivots;
			return lattice;
		}

#ifdef VERBOSE_MS
		std::cerr << "Element Degrees:" << std::endl;
		for (unsigned int i = 0; i < m; ++i) {
			std::cerr << "\t";
			for (unsigned int j = 0; j < m; ++j) {
				std::cerr << deg(lattice[i][j]) + 1 << "\t";
			}
			std::cerr << std::endl;
		}
		std::cerr << "Pivot Info:" << std::endl;
		for (unsigned int i = 0; i < m; ++i) {
			std::cerr << "pivot[" << i << "]:\tused = " << pivots[i].used
					<< "\trow = " << pivots[i].row << "\tdegree = "
					<< pivots[i].degree << "\tlc = " << pivots[i].lc << std::endl;
		}
		std::cerr << "workrow = " << workrow << std::endl;
		std::cerr << "pIndex =  " << pIndex << std::endl;
		std::cerr << "degree =  " << degree << std::endl;
		std::cerr << "lc =      " << lc << std::endl;
#endif

		while (pivots[pIndex].used) {
			if (degree < pivots[pIndex].degree) {
				// Swap current pivot with pivots[pIndex]
				unsigned int tmp_workrow = workrow;
				unsigned int tmp_degree = degree;
				F tmp_lc = lc;
				workrow = pivots[pIndex].row;
				degree = pivots[pIndex].degree;
				lc = pivots[pIndex].lc;
				pivots[pIndex].row = tmp_workrow;
				pivots[pIndex].degree = tmp_degree;
				pivots[pIndex].lc = tmp_lc;
			}

			// Perform the tranformation
			FX transConst (degree - pivots[pIndex].degree, lc / pivots[pIndex].lc);
			for (unsigned int i = 0; i < m; ++i) {
				lattice[workrow][i] -= lattice[pivots[pIndex].row][i] * transConst;
			}
			if (use_aux) {
				for (unsigned int i = 0; i < aux[0].size(); ++i) {
					aux[workrow][i] -= aux[pivots[pIndex].row][i] * transConst;
				}
			}

			// Find new pivot
			pIndex = 0;
			degree = deg(lattice[workrow][0]) + 1;
			for (unsigned int i = 1; i < m; ++i) {
				if ((unsigned int)(deg(lattice[workrow][i]) + 1) >= degree) {
					pIndex = i;
					degree = deg(lattice[workrow][i]) + 1;  // Use +1 to avoid -1 degree
				}
			}
			lc = LeadCoeff(lattice[workrow][pIndex]);

			// Check if we can stop
			if (reduceTo >= 0 && degree <= (unsigned int)reduceTo) {
				small_rows[workrow] = true;
			}
			if (small_rows.size() >= reduceNumber) {
				delete[] pivots;
				return lattice;
			}

#ifdef VERBOSE_MS
			std::cerr << "Element Degrees:" << std::endl;
			for (unsigned int i = 0; i < m; ++i) {
				std::cerr << "\t";
				for (unsigned int j = 0; j < m; ++j) {
					std::cerr << deg(lattice[i][j]) + 1 << "\t";
				}
				std::cerr << std::endl;
			}
			std::cerr << "Pivot Info:" << std::endl;
			for (unsigned int i = 0; i < m; ++i) {
				std::cerr << "pivot[" << i << "]:\tused = " << pivots[i].used
						<< "\trow = " << pivots[i].row << "\tdegree = "
						<< pivots[i].degree << "\tlc = " << pivots[i].lc << std::endl;
			}
			std::cerr << "workrow = " << workrow << std::endl;
			std::cerr << "pIndex =  " << pIndex << std::endl;
			std::cerr << "degree =  " << degree << std::endl;
			std::cerr << "lc =      " << lc << std::endl;
#endif
		}
		// Add PivotInfo to pivots
		pivots[pIndex].used = true;
		pivots[pIndex].row = workrow;
		pivots[pIndex].degree = degree;
		pivots[pIndex].lc = lc;
	}

#ifdef VERBOSE_MS
	std::cerr << "Element Degrees:" << std::endl;
	for (unsigned int i = 0; i < m; ++i) {
		std::cerr << "\t";
		for (unsigned int j = 0; j < m; ++j) {
			std::cerr << deg(lattice[i][j]) + 1 << "\t";
		}
		std::cerr << std::endl;
	}
	std::cerr << "Pivot Info:" << std::endl;
	for (unsigned int i = 0; i < m; ++i) {
		std::cerr << "pivot[" << i << "]:\tused = " << pivots[i].used
				<< "\trow = " << pivots[i].row << "\tdegree = "
				<< pivots[i].degree << "\tlc = " << pivots[i].lc << std::endl;
	}
#endif

	delete[] pivots;
	return lattice;
}

// Find all polynomials of degree at most k that agree with the given
// (index,share) pairs at at least t of the n points.  The notation is
// from Venkatesan Guruswami and Madhu Sudan, "Improved Decoding of
// Reed-Solomon and Algebraic-Geometry Codes".
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector< RecoveryPoly<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_gs(unsigned int k,
	unsigned int t, const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares, TestType testType)
{
	FXY P;

	// Check that method is KOTTER or CH_MS
	switch (testType) {
	case KOTTER:
    	P = interpolate_kotter(k, t, goodservers, indices, shares);
    	break;
	case CH_MS:
    	P = interpolate_cohn_heninger(k, t, goodservers, indices, shares);
    	break;
	default:
		std::cerr << "ERROR: Invalid test method for findpolys().\n";
		return vector<RecoveryPoly<FX> >(); // Fail
	}

#if 0
    // cerr << "factor(poly(0";
    for(int j=0;j<=deg(P);++j) {
	FX x = coeff(P,j);
	for(int i=0; i<=deg(x); ++i) {
	    // cerr << " + " << coeff(x,i) << "*x^" << i << "*y^" << j;
	}
    }
#endif
    // cerr << ", [x,y], IntMod(" << p1*p2 << ")));\n\n";
    std::cerr << "Finding roots of resulting polynomial...\n";

    // It turns out that any polynomial phi(x) that we should be
    // returning (since phi is of degree at most k and agrees with the
    // input data on at least t points) is such that (y - phi(x)) is a
    // factor of P(x,y).  So we first find all roots for y of P(x,y)
    // which are polynomials of degree at most k.
    vector<FX> roots = findroots(P, k);
    // cerr << "roots = " << roots << "\n";

    // For each of these roots, check how many input points it agrees
    // with.  If it's at least t, add it to the list of polys to return.
    vector< RecoveryPoly<FX> > polys;
    unsigned int numroots = roots.size();
    for (unsigned int i=0; i<numroots; ++i) {
	if (deg(roots[i]) > (long)k) continue;
	vector<unsigned short>::const_iterator gooditer;
	vector<unsigned short> vecagree;
	unsigned short numagree = 0;
	for (gooditer = goodservers.begin(); gooditer != goodservers.end();
		++gooditer) {
	    F phival;
	    eval(phival, roots[i], indices[*gooditer]);
	    if (phival == shares[*gooditer]) {
		++numagree;
		vecagree.push_back(*gooditer);
	    }
	}
	if (numagree >= t) {
	    RecoveryPoly<FX> n(vecagree, roots[i]);
	    polys.push_back(n);
	}
    }

    return polys;
}

// Find all polynomials of degree at most k that agree with the given
// (index,share) pairs at at least t of the n points.  This version
// simply uses brute force and interpolates all subsets of size k+1 of
// the n points.  Note that in general, this version does *not* run in
// polynomial time, but for some cases (with n-k small, for example) it
// is faster than Guruswami-Sudan.
//
// Runtime is about C(k,t+1) *
// [2.55712398139188+1.17564033117597*k+4.20999565858028*t+1.21270558067138*t^2]
// microseconds, according to observations and regression analysis.
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<RecoveryPoly<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_brute(
	int k_signed, unsigned int t,
	const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
    //std::cerr << "RUNNING BRUTE\n";

	if (k_signed < -1) {
		return vector<RecoveryPoly<FX> >();  // Fail
	}

	// k == -1
	if (k_signed == -1) {
		// Check if at least h shares are zero
		vector<unsigned short>::const_iterator gooditer;
		vector<unsigned short> vecagree;
		unsigned short numagree = 0;
		for (gooditer = goodservers.begin(); gooditer != goodservers.end();
			++gooditer) {
			if (F::zero() == shares[*gooditer]) {
				++numagree;
				vecagree.push_back(*gooditer);
			}
		}
		if (numagree >= t) {
			return vector< RecoveryPoly<FX> >(1, RecoveryPoly<FX>(vecagree, FX::zero()));
		} else {
			return vector< RecoveryPoly<FX> >();  // Fail
		}
	}

    //std::cerr << "CHECK 1\n";

	// k >= 0
	unsigned int k = (unsigned int)k_signed;
    vector<RecoveryPoly<FX> > polys;

unsigned int numcombs = 0;
    // Iterate over all subsets of goodservers of size k+1
    subset_iterator iter(goodservers, k+1);

    //std::cerr << "CHECK 2\n";

    while(!iter.atend()) {
        ++numcombs;
    
        //std::cerr << "CHECK 2 a\n";

        unsigned short numagree, numdisagree;
        vector<unsigned short> vecagree;
        FX phi;

        // We should probably avoid calling this if polys[i].G is a
        // superset of *iter for some i, but "is a superset" is a
        // non-trivial test with the current data structure, so we'll
        // just run it anyway for now, and check for uniqueness of phi
        // on the way out.
        test_interpolate(k, shares, indices, *iter, goodservers,
            numagree, numdisagree, vecagree, phi);

        //std::cerr << "CHECK 2 b\n";

        if (numagree >= t) {
            // As above: check to see if we've seen this phi before
            typename vector<RecoveryPoly<FX> >::iterator polysiter;
            bool matched = false;
            for (polysiter = polys.begin(); polysiter != polys.end();
                ++polysiter) {
            if (polysiter->phi == phi) {
                matched = true;
                break;
            }
            }
            if (!matched) {
            RecoveryPoly<FX> n(vecagree, phi);
            polys.push_back(n);
            }
            // See if this is the only possible solution.  If so, stop
            if (t > k + numdisagree) {
            break;
            }
        }

        //std::cerr << "CHECK 2 c\n";

        ++iter;
    }

    //std::cerr << "CHECK 3\n";

    return polys;
}

// Find a polynomial of degree ell that agrees with at least h shares.  Uses
// the Berlekamp-Welch algorithm.  This algorithm is only guaranteed to work
// when 2h > n + ell.
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<RecoveryPoly<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_bw(
	unsigned int ell, unsigned int h,
	const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
	unsigned int n = goodservers.size();
	unsigned int numcols = n + n - h - h + ell + 2;
	mat_F M(INIT_SIZE, n, numcols);

	for (unsigned int i = 0; i < n; ++i) {
		F multiplier = indices[goodservers[i]];
		M[i][0] = 1;
		for (unsigned int j = 1; j <= n - h + ell; ++j) {
			M[i][j] = M[i][j-1] * multiplier;
		}
		M[i][n-h+ell+1] = -(shares[goodservers[i]]);
		for (unsigned int j = n - h + ell + 2; j < numcols; ++j) {
			M[i][j] = M[i][j-1] * multiplier;
		}
	}

#ifdef VERBOSE_BW
	std::cerr << "n = " << n << std::endl;
	std::cerr << "ell = " << ell << std::endl;
	std::cerr << "h = " << h << std::endl;
	std::cerr << "indices =";
	for (unsigned int i = 0; i < n; ++i) {
		std::cerr << " " << indices[goodservers[i]];
	}
	std::cerr << std::endl;
	std::cerr << "shares =";
	for (unsigned int i = 0; i < n; ++i) {
		std::cerr << " " << shares[goodservers[i]];
	}
	std::cerr << std::endl;
	std::cerr << "numcols = " << numcols << std::endl;
	std::cerr << "M = " << M << std::endl;
#endif

	unsigned int rank = gauss(M);

#ifdef VERBOSE_BW
	std::cerr << "rank = " << rank << std::endl;
	std::cerr << "reduced M = " << M << std::endl;
#endif

	if (rank == numcols) {
#ifdef VERBOSE_BW
		std::cerr << "FAIL: y_i*E(x_i)=Q)x_i) has no solution" << std::endl;
#endif
		return vector<RecoveryPoly<FX> >();
	}

	unsigned int findex;
	for (findex = 0; findex < n; ++findex) {
		if (IsZero(M[findex][findex])) {
			break;
		}
	}

#ifdef VERBOSE_BW
	std::cerr << "findex = " << findex << std::endl;
#endif

	vec_F coeffs(INIT_SIZE, numcols);
	coeffs[findex] = 1;
	for (unsigned int i = findex; i > 0; ) {
		--i;  // To avoid checking that an unsigned int is >=0
		F c;
		for (unsigned int j = i + 1; j <= findex; ++j) {
			c -= M[i][j] * coeffs[j];
		}
		c /= M[i][i];
		coeffs[i] = c;
	}

#ifdef VERBOSE_BW
	std::cerr << "coeffs:" << std::endl;
	for (unsigned int i = 0; i < numcols; ++i) {
		std::cerr << "\t" << i << "\t" << coeffs[i] << std::endl;
	}
#endif

	FX Q, E;
	for (unsigned int i = 0; i <= n - h + ell; ++i) {
		SetCoeff(Q, i, coeffs[i]);
	}
	for (unsigned int i = 0; i <= n - h; ++i) {
		SetCoeff(E, i, coeffs[n - h + ell + 1 + i]);
	}

#ifdef VERBOSE_BW
	std::cerr << "Q = " << Q << std::endl;
	std::cerr << "E = " << E << std::endl;

	// Check that y_i*E(x_i)=Q(x_i)
	for (unsigned int i = 0; i < n; ++i) {
		if (shares[goodservers[i]] * eval(E, indices[goodservers[i]]) != eval(Q, indices[goodservers[i]])) {
			std::cerr << "ERROR: check failed for i = " << i << std::endl;
		}
	}
#endif

	FX result;
	if (!divide(result, Q, E)) {
#ifdef VERBOSE_BW
		std::cerr << "FAIL: Q % E != 0" << std::endl;
#endif
		return vector<RecoveryPoly<FX> >();
	}

	if (deg(result) > 0 && ((unsigned int)deg(result)) > ell) {
#ifdef VERBOSE_BW
		std::cerr << "FAIL: deg(result) > " << ell << std::endl;
#endif
		return vector<RecoveryPoly<FX> >();
	}

#ifdef VERBOSE_BW
	std::cerr << "result = " << result << std::endl;
#endif

	// For result, check how many input points it agrees
	// with.  If it's at least t, add it to the list of polys to return.
	vector< RecoveryPoly<FX> > polys;
	vector<unsigned short>::const_iterator gooditer;
	vector<unsigned short> vecagree;
	unsigned short numagree = 0;
	for (gooditer = goodservers.begin(); gooditer != goodservers.end();
		++gooditer) {
		F phival;
		eval(phival, result, indices[*gooditer]);
		if (phival == shares[*gooditer]) {
			++numagree;
			vecagree.push_back(*gooditer);
		}
	}
	if (numagree >= h && 2*numagree > n+ell) {
		RecoveryPoly<FX> n(vecagree, result);
		polys.push_back(n);
	}

	return polys;
}

/*
// Gaussian elimination of a matrix over FX.  Result leaves matrix in row
// echelon form.  This method is in-place (changes the given matrix).  The rank
// of the matrix is returned.
template<class F, class vec_F, class FX, class FXY, class mat_F>
int gaussian_FX (vector<vector<FX> >& A, vector<FX>& b) {
	vector<vector<FX> > mat = matrix;
    std::cerr << "CHECKPOINT 2a\n";
	unsigned n = mat.size();
    std::cerr << "CHECKPOINT 2b\n";
	unsigned m = mat[0].size();
    std::cerr << "CHECKPOINT 2c\n";
	unsigned int h = 0;
    std::cerr << "CHECKPOINT 2d\n";
	for (unsigned int i = 0; i < n && h < m; ++i) {
        std::cerr << "CHECKPOINT 2d 1\n";
		FX elem_ih = mat[i][h];
		// Find pivot
        std::cerr << "CHECKPOINT 2d 2\n";
		while (IsZero(elem_ih)) {
			if (h == m) {
				return mat;
			}
			unsigned int j;
			for (j = i + 1; j < n; ++j) {
				if (!IsZero(mat[j][h])) {
					vector<FX> tmp = mat[j];
					mat[j] = mat[i];
					mat[i] = tmp;
					break;
				}
			}
			if (j == n) {
				++h;
			}
			elem_ih = mat[i][h];
		}
        std::cerr << "CHECKPOINT 2d 3\n";
		for (unsigned int j = i + 1; j < n; ++j) {
			FX elem_jh = mat[j][h];
			if (elem_jh != 0) {
				// (Row j) <- (Row j) * M_ih / gcd(M_ih, M_jh) - (Row i) * M_jh / gcd(M_ih, M_jh)
				FX gcd = GCD(elem_ih, elem_jh);
				for (unsigned int k = h; k < m; ++k) {
					mat[j][k] = mat[j][k] * (elem_ih / gcd) - mat[i][k] * (elem_jh / gcd);
				}
			}
			elem_ih = mat[i][h];
		}
        std::cerr << "CHECKPOINT 2d 4\n";
		++h;
	}
	return mat;
}
*/

// Solve the linear system of equations Ax=b.  
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<FX> RSDecoder<F,vec_F,FX,FXY,mat_F>::solve_linsys_FX (
        vector<vector<FX> >& A, vector<FX>& b, unsigned int modulus) {

	// Sanity checks
	unsigned int n = A.size();  // Number of rows
	if (n == 0) {
		std::cerr << "FAIL: invalid_linsys: A is empty\n";
		return vector<FX>();  // Fail
	}
	unsigned int m = A[0].size();  // Number of cols
	if (b.size() != n) {
		std::cerr << "FAIL: invalid_linsys: Size of b does not match dimensions of A\n";
		return vector<FX>();  // Fail
	}

    std::vector<unsigned int> freeCols;
    unsigned int h = 0;
    for (unsigned int i = 0; i < n && h < m; ++i) {
        while (IsZero(A[i][h])) {
            //std::cerr << "(i, h) = (" << i << ", " << h << ")\n";
            // Find pivot for column h and put it in row i
            unsigned int j;
            for (j = i+1; j < n; ++j) {
                if (!IsZero(A[j][h])) {
                    // Swap rows i and j
                    vector<FX> tmp_row = A[i];
                    A[i] = A[j];
                    A[j] = tmp_row;
                    FX tmp_elem = b[i];
                    b[i] = b[j];
                    b[j] = tmp_elem;
                    break;
                }
            }

            // Did we find a pivot?
            if (j == n) {
                freeCols.push_back(h);
                ++h;
                // What if h == m now?
                if (h == m) {
                    goto done_upper_triangular;
                }
            }
        }

        // Do row operations
        for (unsigned int j = i+1; j < n; ++j) {
            if (!IsZero(A[j][h])) {
                // (Row j) <-- (Row j) * A_ih / gcd(A_ih, Ajh) 
                //             - (Row i) * A_jh / gcd(A_ih, A_jh)
                FX gcd = GCD(A[i][h], A[j][h]);
                FX mult1 = A[i][h] / gcd;
                FX mult2 = A[j][h] / gcd;
                if (modulus > 0) {
                    for (unsigned int k = h; k < m; ++k) {
                        A[j][k] = MulTrunc(A[j][k], mult1, modulus)
                                    - MulTrunc(A[i][k], mult2, modulus);
                    }
                    b[j] = MulTrunc(b[j], mult1, modulus) 
                            - MulTrunc(b[i], mult2, modulus);
                } else {
                    for (unsigned int k = h; k < m; ++k) {
                        A[j][k] = (A[j][k] * mult1) - (A[i][k] * mult2);
                    }
                    b[j] = (b[j] * mult1) - (b[i] * mult2);
                }
            }
        }
        ++h;
    }

done_upper_triangular:
    //std::cerr << "A | b =\n";
    //for (unsigned int i = 0; i < n; ++i) {
    //    std::cerr << "\t";
    //    for (unsigned int j = 0; j < n; ++j) {
    //        std::cerr << deg(A[i][j]) << "\t";
    //    }
    //    std::cerr << "|\t" << b[i] << "\n";
    //}
    //std::cerr << "freeCols = { ";
    //for (unsigned int i = 0; i < freeCols.size(); ++i) {
    //    std::cerr << freeCols[i] << " ";
    //}
    //std::cerr << "}\n";

	// Check that there is a solution (consistency)
	unsigned int rank = n - freeCols.size();
	for (unsigned int i = rank; i < n; ++i) {
		if (!IsZero(b[i])) {
			std::cerr << "FAIL: no_solution: inconsistent matrix\n";
			return vector<FX>();  // Fail
		}
	}

	// Get solution (free variables have value zero)
	vector<FX> solution (m, FX::zero());
    unsigned int sub = freeCols.size();
    for (unsigned int i = m; i > 0;) {
        --i;
        if (sub > 0 && freeCols[sub-1] == i) {
            --sub;
        } else {
            solution[i] = b[i-sub];
            for (unsigned int j = i+1; j < m; ++j) {
                solution[i] -= A[i-sub][j] * solution[j];
            }
            solution[i] = solution[i] / A[i-sub][i];
        }
    }

    //std::cerr << "solution =\n";
    //for (unsigned int i = 0; i < solution.size(); ++i) {
    //    std::cerr << "\t" << solution[i] << "\n";
    //}

	return solution;
}

#ifdef INCLUDE_GENERAL_MULTI
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<map<unsigned int, FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::solve_partial_solution (
		unsigned int ell, unsigned int h, map<unsigned int, FX> partial,
		const vector<vector<FX> >& rows, vector<unsigned int> row_order,
		const vector<vector<unsigned int> >& order_of_expts)
{
//	std::cerr << "solving partial: ";
//	typename map<unsigned int, FX>::iterator it;
//	for (it = partial.begin(); it != partial.end(); ++it) {
//		std::cerr << "x" << it->first << " = " << it->second << ", ";
//	}
//	std::cerr << "\n";

	unsigned int dim = order_of_expts.size();
	if (dim == 0) {
#ifdef VERBOSE_CH_MULTI
		std::cerr << "FAIL: order_of_expts is empty\n";
#endif
		return vector<map<unsigned int, FX> >();  // Fail
	}
	unsigned int m = order_of_expts[0].size();

	// Have a solution.  Check that it is consistent.
	if (partial.size() == m) {
//		std::cerr << "** BASE **\n";
		for (unsigned int i = 0; i < row_order.size(); ++i) {
			vector<FX> row = rows[row_order[i]];
			FX cterm;
			for (unsigned int j = 0; j < dim; ++j) {
				if (IsZero(row[j])) {
					continue;
				}
				FX current = row[j];
				vector<unsigned int> expts = order_of_expts[j];
				for (unsigned int p = 0; p < m; ++p) {
					current *= power(partial[p], expts[p]);
				}
				cterm += current;
			}
//			std::cerr << "cterm = " << cterm << "\n";
			if (!IsZero(cterm)) {
				return vector<map<unsigned int, FX> >();  // Not a solution.
			}
		}
		return vector<map<unsigned int, FX> >(1, partial);
	}

	// Solve another variable.
	for (unsigned int i = 0; i < row_order.size(); ++i) {
//		std::cerr << "i = " << i << "\n";
		vector<FX> row = rows[row_order[i]];
		// See if row with partial solution is in F[z][xi] for some i
		FXY tosolve;
		unsigned int tosolve_var = 0;
		bool has_tosolve_var = false;
		bool too_many_vars = false;
		for (unsigned int j = 0; j < dim; ++j) {
//			std::cerr << "j = " << j << "\n";
			if (IsZero(row[j])) {
				continue;
			}
			vector<unsigned int> expts = order_of_expts[j];
			FX coeff = row[j];
			for (unsigned int p = 0; p < m; ++p) {
//				std::cerr << "p = " << p << "\n";
				if (expts[p] == 0) {
					continue;
				}
				typename map<unsigned int, FX>::iterator partialIter = partial.find(p);
				if (partialIter == partial.end()) {
					// Unknown variable
					if (has_tosolve_var && tosolve_var != p) {
						too_many_vars = true;
						break;
					}
					has_tosolve_var = true;
					tosolve_var = p;
				} else {
					coeff *= power(partialIter->second, expts[p]);
				}
			}
			if (too_many_vars) {
				break;
			}
			if (has_tosolve_var) {
				tosolve += FXY(expts[tosolve_var], coeff);
			} else {
				tosolve += FXY(0, coeff);
			}
			if (tosolve == 0) {
				// Row equals zero
				has_tosolve_var = false;
			}
		}
//#ifdef VERBOSE_CH_MULTI
//		std::cerr << "row # = " << row_order[i] << "\n";
//#endif
		if (too_many_vars) {
			continue;
		}
//#ifdef VERBOSE_CH_MULTI
//		std::cerr << "tosolve = " << tosolve << "\n";
//#endif
		row_order.erase(row_order.begin() + i);
		if (!has_tosolve_var) {
			// No unknowns
			if (!IsZero(tosolve)) {
				return vector<map<unsigned int, FX> >();  // Not a solution.
			}
			--i;
			continue;
		}
		// Can solve
		vector<FX> roots = findroots(tosolve, ell);
//#ifdef VERBOSE_CH_MULTI
//		std::cerr << "roots =";
//		for (unsigned int j = 0; j < roots.size(); ++j) {
//			std::cerr << " " << roots[j];
//		}
//		std::cerr << "\n";
//#endif
		vector<map<unsigned int, FX> > result;
		for (unsigned int j = 0; j < roots.size(); ++j) {
			partial[tosolve_var] = roots[j];
			vector<map<unsigned int, FX> > root_result =
					solve_partial_solution(ell, h, partial, rows, row_order, order_of_expts);
			result.insert(result.end(), root_result.begin(), root_result.end());
		}
		return result;
	}
	if (partial.size() > 0) {
		return vector<map<unsigned int, FX> >(1, partial);  // Return partial solution
	}
	return vector<map<unsigned int, FX> >();  // No vars could be solved
}

template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<RecoveryPolyMulti<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_ch_multi(
		unsigned int ell, unsigned int h, unsigned int m, unsigned int& t,
		unsigned int& k,
		const vector<unsigned short> &goodservers,
		const vec_F &indices, const vector<vec_F> &shares) {

#ifdef VERBOSE_CH_MULTI
	std::cerr << "Starting Cohn-Heninger Multi-Polynomial algorithm...\n";
#endif

	unsigned int n = goodservers.size();
	unsigned int dim = 1;
	for (unsigned int i = 0; i < (m<t?m:t); ++i) {
		dim *= (m + t - i);
		dim /= (i + 1);
	}

#ifdef VERBOSE_CH_MULTI
	std::cerr << "n = " << n << "\n";
	std::cerr << "ell = " << ell << "\n";
	std::cerr << "h = " << h << "\n";
	std::cerr << "m = " << m << "\n";
	std::cerr << "t = " << t << "\n";
	std::cerr << "k = " << k << "\n";
	std::cerr << "dim = " << dim << "\n";
#endif

#ifdef TIME_PARTS
	std::stringstream tsout;
	tsout << n << "," << ell << "," << h << ",";
#endif

	// Create N(z)
	FX N;
	SetCoeff(N, 0, 1);
	for (unsigned int i = 0; i < n; ++i) {
		FX newTerm;
		SetCoeff(newTerm, 0, -(indices[goodservers[i]]));
		SetCoeff(newTerm, 1);
		N *= newTerm;
	}

#ifdef VERBOSE_CH_MULTI
	std::cerr << "N(z) = " << N << " (" << deg(N) << ")\n";
#endif

	// Create L[i](z), the interpolation of points in shares[i]
	vector<FX> L;
	vec_F goodIndices(INIT_SIZE, n);
	for (unsigned int j = 0; j < n; ++j) {
		goodIndices[j] = indices[goodservers[j]];
	}
	for (unsigned int i = 0; i < m; ++i) {
		vec_F goodShares(INIT_SIZE, n);
		for (unsigned int j = 0; j < n; ++j) {
			goodShares[j] = shares[i][goodservers[j]];
		}
		FX L_i;
		interpolate(L_i, goodIndices, goodShares);
		L.push_back(L_i);
	}

#ifdef VERBOSE_CH_MULTI
	for (unsigned int i = 0; i < m; ++i) {
		std::cerr << "L[" << i << "](z) = " << L[i] << " (" << deg(L[i]) << ")\n";
	}
#endif

	vector<FX> powers_of_N;
	powers_of_N.push_back(FX(0, 1));
	for (unsigned int i = 1; i <= k; ++i) {
		powers_of_N.push_back(powers_of_N[i-1] * N);
	}

#ifdef VERBOSE_CH_MULTI
	for (unsigned int i = 0; i <= k; ++i) {
		std::cerr << "powers_of_N[" << i << "](z) = " << powers_of_N[i] << ")\n";
	}
#endif

	// Create lattice basis
	vector<vector<FX> > lattice(m+1, vector<FX>(m+1));
	lattice[0][0] = N;
	FX z_toell(ell, 1);
	for (unsigned int i = 0; i < m; ++i) {
		lattice[i+1][0] = -(L[i]);
		lattice[i+1][i+1] = z_toell;
	}

#ifdef VERBOSE_CH_MULTI
	std::cerr << "lattice =\n";
	for (unsigned int i = 0; i <= m; ++i) {
		for (unsigned int j = 0; j <= m; ++j) {
			std::cerr << "\t" << i << "," << j << ":\t" << lattice[i][j] << "\n";
		}
	}
#endif

	vector<vector<FX> > aux (m+1, vector<FX>());

#ifdef TIME_PARTS
	struct timeval st, et;
	gettimeofday(&st, NULL);
#endif

	vector<vector<FX> > reducedLattice = reduce_lattice_MS(lattice, aux);

#ifdef TIME_PARTS
	gettimeofday(&et, NULL);
	unsigned long long elapsedus = ((unsigned long long)(et.tv_sec -
			st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
	tsout << "lattice reduction," << elapsedus << ",";
#endif

#ifdef VERBOSE_CH_MULTI
	std::cerr << "reducedLattice =\n";
	for (unsigned int i = 0; i < m+1; ++i) {
		for (unsigned int j = 0; j < m+1; ++j) {
			std::cerr << "\t" << i << "," << j << ":\t" << reducedLattice[i][j] << "\n";
		}
	}
	std::cerr << "reducedLattice degrees:\n";
	for (unsigned int i = 0; i < m+1; ++i) {
		std::cerr << "\t";
		for (unsigned int j = 0; j < m+1; ++j) {
			std::cerr << deg(reducedLattice[i][j]) << "\t";
		}
		std::cerr << "\n";
	}
#endif

	vector<unsigned int> smallRows;
	for (unsigned int i = 0; i < m+1; ++i) {
		int maxdeg = -2;
		for (unsigned int j = 0; j < m+1; ++j) {
			if (deg(reducedLattice[i][j]) > maxdeg) {
				maxdeg = deg(reducedLattice[i][j]);
			}
 		}
//		std::cerr << "i = " << i << ", maxdeg = " << maxdeg << ", h = " << h << "\n";
		if (maxdeg < (int)(h)) {
			smallRows.push_back(i);
		}
	}

#ifdef VERBOSE_CH_MULTI
	std::cerr << "smallRows =";
	vector<unsigned int>::iterator sriter;
	for (sriter = smallRows.begin(); sriter != smallRows.end(); ++sriter) {
		std::cerr << " " << *sriter;
	}
	std::cerr << "\n";
#endif

	if (smallRows.size() < m) {

		// If t=k=1, then fail
		if (t == 1 && k == 1) {
#ifdef VERBOSE_CH_MULTI
			std::cerr << "FAIL: explicit t=k=1 failed\n";
#endif
			return vector<RecoveryPolyMulti<FX> >();  // Fail
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "Not enough small vectors with t=k=1.  Trying larger values...\n";
#endif

		vector<vector<vector<FX> > > coeff_by_expts;  // coeff_by_expts[i][j][p] = [x^p](x_i - L_i(z))^j
		for (unsigned int i = 0; i < m; ++i) {
			vector<vector<FX> > coeffs_for_var;
			coeffs_for_var.push_back(vector<FX>(1, FX(0, 1)));
			for (unsigned int j = 1; j <= t; ++j) {
				vector<FX> current_j;
				for (unsigned int p = 0; p <= j; ++p) {
					FX current;
					if (p == 0) {
						current = coeffs_for_var[j-1][p] * (-(L[i]));
					} else if (p == j) {
						current = 1;
					} else {
						current = coeffs_for_var[j-1][p] * (-(L[i]));
						current += coeffs_for_var[j-1][p-1];
					}
					current_j.push_back(current);
				}
				coeffs_for_var.push_back(current_j);
			}
			coeff_by_expts.push_back(coeffs_for_var);
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "coeff_by_expts =\n";
		for (unsigned int i = 0; i < m; ++i) {
			for (unsigned int j = 0; j <= t; ++j) {
				for (unsigned int p = 0; p <= j; ++p) {
					std::cerr << "\t" << i << "," << j << "," << p << "\t" << coeff_by_expts[i][j][p] << "\n";
				}
			}
		}
#endif

		vector<vector<unsigned int> > order_of_expts;
		order_of_expts.push_back(vector<unsigned int>(m, 0));
		vector<vector<vector<unsigned int> > > prev_t;
		for (unsigned int i = 0; i < m; ++i) {
			vector<unsigned int> current (m, 0);
			current[i] = 1;
			prev_t.push_back(vector<vector<unsigned int> >(1, current));
			order_of_expts.push_back(current);
		}
		for (unsigned i = 2; i <= t; ++i) {
			vector<vector<vector<unsigned int> > > curr_t;
			for (unsigned int j = 0; j < m; ++j) {
				vector<vector<unsigned int> > current_j;
				for (unsigned int p = j; p < m; ++p) {
					vector<vector<unsigned int> >::iterator iter;
					for (iter = prev_t[p].begin(); iter != prev_t[p].end(); ++iter) {
						vector<unsigned int> current = *iter;
						current[j] += 1;
						current_j.push_back(current);
						order_of_expts.push_back(current);
					}
				}
				curr_t.push_back(current_j);
			}
			prev_t = curr_t;
		}

		vector<unsigned int> exptSum;
		for (unsigned int i = 0; i < dim; ++i) {
			unsigned int current = 0;
			vector<unsigned int> expts = order_of_expts[i];
			for (unsigned int p = 0; p < m; ++p) {
				current += expts[p];
			}
			exptSum.push_back(current);
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "order_of_expts, exptSum =\n";
		for (unsigned int i = 0; i < dim; ++i) {
			std::cerr << "\t";
			for (unsigned int j = 0; j < m; ++j) {
				std::cerr << order_of_expts[i][j] << " ";
			}
			std::cerr << ", " << exptSum[i] << "\n";
		}
#endif

		vector<vector<FX> > lattice2;
		for (unsigned int i = 0; i < dim; ++i) {
			vector<unsigned int> rowExpts = order_of_expts[i];
			vector<FX> row;
			for (unsigned j = 0; j < dim; ++j) {
				vector<unsigned int> termExpts = order_of_expts[j];
				bool nonZero = true;
				for (unsigned int p = 0; p < m; ++p) {
					if (termExpts[p] > rowExpts[p]) {
						nonZero = false;
						break;
					}
				}
				if (!nonZero) {
					row.push_back(FX());
					continue;
				}
				FX termCoeff = FX(0, 1);
				for (unsigned int p = 0; p < m; ++p) {
					termCoeff *= coeff_by_expts[p][rowExpts[p]][termExpts[p]];
				}
				termCoeff *= power(z_toell, exptSum[j]);
				if (exptSum[i] < k) {
					termCoeff *= powers_of_N[k - exptSum[i]];
				}
				row.push_back(termCoeff);
			}
			lattice2.push_back(row);
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "lattice2 =\n";
		for (unsigned int i = 0; i < dim; ++i) {
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << "\t" << i << "," << j << "\t" << lattice2[i][j] << "\n";
			}
		}
		std::cerr << "lattice2 non-zero elements:\n";
		for (unsigned int i = 0; i < dim; ++i) {
			std::cerr << "\t";
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << (IsZero(lattice2[i][j]) ? " " : "X");
			}
			std::cerr << "\n";
		}
#endif

		// Reduce lattice basis using Mulders-Storjohann
		vector<vector<FX> > aux (dim, vector<FX>());
		//vector<vector<FX> > reducedLattice2WSpoon = reduce_lattice_MS(lattice2, aux, h*k, m);
		vector<vector<FX> > reducedLattice2WSpoon = reduce_lattice_MS(lattice2, aux);

#ifdef VERBOSE_CH_MULTI
		std::cerr << "reducedLattice2WSpoon =\n";
		for (unsigned int i = 0; i < dim; ++i) {
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << "\t" << i << "," << j << ":\t" << reducedLattice2WSpoon[i][j] << "\n";
			}
		}

		std::cerr << "reducedLattice2WSpoon degrees:\n";
		for (unsigned int i = 0; i < dim; ++i) {
			std::cerr << "\t";
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << deg(reducedLattice2WSpoon[i][j]) << "\t";
			}
			std::cerr << "\n";
		}
#endif

		// Take out the wooden spoon
		vector<vector<FX> > reducedLattice2;
		for (unsigned int i = 0; i < dim; ++i) {
			vector<FX> row;
			for (unsigned int j = 0; j < dim; ++j) {
//				std::cerr << "Removing spoon from " << i << "," << j << "\n";
//				std::cerr << "\torder_of_expts[" << j << "], exptSum[" << j << "] =";
//				for (unsigned int p = 0; p < m; ++p) {
//					std::cerr << " " << order_of_expts[j][p];
//				}
//				std::cerr << ", " << exptSum[j] << "\n";
//				std::cerr << "\treducedLattice2WSpoon[" << i << "][" << j << "] = " << reducedLattice2WSpoon[i][j] << "\n";
//				std::cerr << "\tz_toell^exptSum[" << j << "] = " << power(z_toell, exptSum[j]) << "\n";
//				std::cerr << "\treducedLattice2WSpoon[" << i << "][" << j << "] / power(z_toell, exptSum[" << j << "]) = " << reducedLattice2WSpoon[i][j] / power(z_toell, exptSum[j]) << "\n";
				row.push_back(reducedLattice2WSpoon[i][j] / power(z_toell, exptSum[j]));
			}
			reducedLattice2.push_back(row);
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "reducedLattice2 =\n";
		for (unsigned int i = 0; i < dim; ++i) {
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << "\t" << i << "," << j << ":\t" << reducedLattice2[i][j] << "\n";
			}
		}

		std::cerr << "reducedLattice2 degrees:\n";
		for (unsigned int i = 0; i < dim; ++i) {
			std::cerr << "\t";
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << deg(reducedLattice2[i][j]) << "\t";
			}
			std::cerr << "\n";
		}
#endif

		map<int, unsigned int> reducedDegrees;
		vector<unsigned int> smallRows;
		for (unsigned int i = 0; i < dim; ++i) {
			int maxdeg = -2;
			for (unsigned int j = 0; j < dim; ++j) {
				if (deg(reducedLattice2[i][j]) > maxdeg) {
					maxdeg = deg(reducedLattice2[i][j]);
				}
			}
			if (maxdeg < (int)(k*h)) {
				smallRows.push_back(i);
			}
			map<int, unsigned int>::iterator it = reducedDegrees.find(maxdeg);
			if (it != reducedDegrees.end()) {
				reducedDegrees[maxdeg] += 1;
			} else {
				reducedDegrees[maxdeg] = 1;
			}
		}

//		std::cerr << "reducedDegrees =\n";
//		for (map<int, unsigned int>::iterator it = reducedDegrees.begin();
//				it != reducedDegrees.end(); ++it) {
//			std::cerr << "\t" << it->first << "\t" << it->second << "\n";
//		}

		std::cerr << "reducedLattice2 max degrees:";
		for (map<int, unsigned int>::iterator it = reducedDegrees.begin();
				it != reducedDegrees.end(); ++it) {
			for (unsigned int j = 0; j < it->second; ++j) {
				std::cerr << " " << it->first;
			}
		}
		std::cerr << "\n";

#ifdef VERBOSE_CH_MULTI
		std::cerr << "smallRows =";
		vector<unsigned int>::iterator sriter;
		for (sriter = smallRows.begin(); sriter != smallRows.end(); ++sriter) {
			std::cerr << " " << *sriter;
		}
		std::cerr << " (" << smallRows.size() << ")\n";
#endif

		if (smallRows.size() < m) {
			std::cerr << "FAIL: few_small_vectors: Not enough small vectors\n";
			return vector<RecoveryPolyMulti<FX> >();  // Fail
		}

		std::stringstream tosage;
		for (vector<unsigned int>::iterator iter = smallRows.begin();
				iter != smallRows.end();
				++iter) {
			tosage << "{" << reducedLattice2[*iter][0];
			for (unsigned int j = 1; j < dim; ++j) {
				tosage << "," << reducedLattice2[*iter][j];
			}
			tosage << "}\n";
		}

		const std::string tosageFile = "/tmp/tosage.txt";
		const std::string fromsageFile = "/tmp/fromsage.txt";
		std::ofstream tosageStream(tosageFile.c_str());
		tosageStream << tosage.str();
		tosageStream.close();

		std::stringstream cmd;
		cmd << "sage solveGroebnerBasis.sage " << m << " " << t;
// This needs to be changed.  Possibly use typeid from typeinfo
#ifdef USE_GF28
		cmd << " gf28";
#elif defined USE_GF24
		cmd << " gf24";
#elif defined USE_W8
		cmd << " w8";
#elif defined USE_W32
		cmd << " w32";
#else
		cmd << " w128";
#endif
		cmd << " < " << tosageFile << " > " << fromsageFile;

		if (system(cmd.str().c_str()) != 0) {
#ifdef VERBOSE_CH_MULTI
			std::cerr << "FAIL: sage failed\n";
#endif
			return vector<RecoveryPolyMulti<FX> >();  // Fail
		}

		vector<vector<FX> > fromsage;
		unsigned int rows;
		std::ifstream fromsageStream(fromsageFile.c_str());
		fromsageStream >> rows;
		for (unsigned int i = 0; i < rows; ++i) {
			vector<FX> row;
			for (unsigned int j = 0; j < dim; ++j) {
				FX current;
				fromsageStream >> current;
				row.push_back(current);
			}
			fromsage.push_back(row);
		}
		fromsageStream.close();

#ifdef VERBOSE_CH_MULTI
		std::cerr << "rows = " << rows << "\n";
		std::cerr << "fromsage =\n";
		for (unsigned int i = 0; i < rows; ++i) {
			for (unsigned int j = 0; j < dim; ++j) {
				std::cerr << "\t" << i << "," << j << "\t" << fromsage[i][j] << "\n";
			}
		}
#endif

		// Vars in each row
		vector<map<unsigned int, bool> > hasVars(rows);
		for (unsigned int i = 0; i < dim; ++i) {
			vector<unsigned int> expts = order_of_expts[i];
			for (unsigned int j = 0; j < rows; ++j) {
				if (fromsage[j][i] != 0) {
					for (unsigned int p = 0; p < m; ++p) {
						if (expts[p] != 0) {
							hasVars[j][p] = true;
						}
					}
				}
			}
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "hasVars =\n";
		for (unsigned int i = 0; i < rows; ++i) {
			std::cerr << "\tRow " << i << "\t";
			typename map<unsigned int, bool>::iterator varIter;
			for (varIter = hasVars[i].begin(); varIter != hasVars[i].end(); ++varIter) {
				std::cerr << varIter->first << " ";
			}
			std::cerr << "\n";
		}
#endif

		// Sort by # of vars
		vector<unsigned int> row_order;
		for (unsigned int i = 0; i <= m; ++i) {
			for (unsigned int j = 0; j < rows; ++j) {
				if (hasVars[j].size() == i) {
					row_order.push_back(j);
				}
			}
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "row_order =";
		for (unsigned int i = 0; i < rows; ++i) {
			std::cerr << " " << row_order[i];
		}
		std::cerr << "\n";
#endif

		// Work to find solutions
		vector<map<unsigned int, FX> >	solutions =
				solve_partial_solution(ell, h, map<unsigned int, FX>(), fromsage, row_order, order_of_expts);

		if (solutions.size() == 0) {
#ifdef VERBOSE_CH_MULTI
			std::cerr << "FAIL: no solutions found\n";
#endif
			return vector<RecoveryPolyMulti<FX> >();  // Fail
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "solutions =\n";
		for (unsigned int i = 0; i < solutions.size(); ++i) {
			for (unsigned int j = 0; j < m; ++j) {
				std::cerr << "\tx" << j << "\t";
				typename map<unsigned int, FX>::iterator it = solutions[i].find(j);
				if (it != solutions[i].end()) {
					std::cerr << solutions[i][j] << " [ ";
					for (unsigned int p = 0; p < n; ++p) {
						F phival;
						eval(phival, solutions[i][j], indices[goodservers[p]]);
						if (phival == shares[j][goodservers[p]]) {
							std::cerr << p << " ";
						}
					}
					std::cerr << "]";

				}
				std::cerr << "\n";
			}
			std::cerr << "\n";
		}
#endif

		vector<RecoveryPolyMulti<FX> > polys;
		for (unsigned int j = 0; j < solutions.size(); ++j) {
			// For each polynomial in the solution, check how many input points it agrees
			// with.  If it's at least h, add it to the list of polys to return.
			map<unsigned int, FX> solution = solutions[j];

			// For now
			if (solution.size() != m) {
				continue;
			}

			vector<unsigned short>::iterator gooditer;
			vector<unsigned short> vecagree = goodservers;
			vector<FX> vecsolution;
			for (unsigned int i = 0; i < m; ++i) {
				for (gooditer = vecagree.begin(); gooditer != vecagree.end(); ) {
					F phival;
					eval(phival, solution[i], indices[*gooditer]);
					if (phival == shares[i][*gooditer]) {
						++gooditer;
					} else {
						gooditer = vecagree.erase(gooditer);
					}
				}
				vecsolution.push_back(solution[i]);
			}
			if (vecagree.size() >= h) {
				RecoveryPolyMulti<FX> n(vecagree, vecsolution);
				polys.push_back(n);
			}
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "polys = \n";
		for (unsigned int i = 0; i < polys.size(); ++i) {
			for (unsigned int j = 0; j < m; ++j) {
				std::cerr << "\t" << j << "\t" << polys[i].phis[j] << "\n";
			}
			std::cerr << "\n";
		}
#endif

		return polys;

	} else {
		// Continue with t=k=1 work
		vector<vector<FX> > A;
		vector<FX> b;
		vector<unsigned int>::iterator smallIter = smallRows.begin();
		for (unsigned int i = 0; smallIter != smallRows.end() && i < m; ++smallIter, ++i) {
			vector<FX> currRow;
			for (unsigned int j = 1; j <= m; ++j) {
				currRow.push_back(RightShift(reducedLattice[*smallIter][j], ell));
			}
			A.push_back(currRow);
			b.push_back(-(reducedLattice[*smallIter][0]));
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "A =\n";
		for (unsigned int i = 0; i < m; ++i) {
			for (unsigned int j = 0; j < m; ++j) {
				std::cerr << "\t" << i << "," << j << ":\t" << A[i][j] << "\n";
			}
		}
		std::cerr << "b =\n";
		for (unsigned int i = 0; i < m; ++i) {
			std::cerr << "\t" << i << ":\t" << b[i] << "\n";
		}
#endif

#ifdef TIME_PARTS
		gettimeofday(&st, NULL);
#endif

		//vector<FX> solution = solve_linsys_FX(A, b, t+1);
		vector<FX> solution = solve_linsys_FX(A, b);

#ifdef TIME_PARTS
		gettimeofday(&et, NULL);
		unsigned long long elapsedus = ((unsigned long long)(et.tv_sec -
				st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
		tsout << "solving linear system," << elapsedus << ",";
#endif

		if (solution.size() == 0) {
			return vector<RecoveryPolyMulti<FX> >();  // Fail
		}

#ifdef VERBOSE_CH_MULTI
		std::cerr << "solution =\n";
		for (unsigned int i = 0; i < m; ++i) {
			std::cerr << "\t" << i << ":\t" << solution[i] << "\n";
		}
#endif

		// For each polynomial in the solution, check how many input points it agrees
		// with.  If it's at least h, add it to the list of polys to return.
		vector< RecoveryPolyMulti<FX> > polys;
		vector<unsigned short>::iterator gooditer;
		vector<unsigned short> vecagree = goodservers;
		for (unsigned int i = 0; i < m; ++i) {
			for (gooditer = vecagree.begin(); gooditer != vecagree.end(); ) {
				F phival;
				eval(phival, solution[i], indices[*gooditer]);
				if (phival == shares[i][*gooditer]) {
					++gooditer;
				} else {
					gooditer = vecagree.erase(gooditer);
				}
			}
		}
		if (vecagree.size() >= h) {
			RecoveryPolyMulti<FX> n(vecagree, solution);
			polys.push_back(n);
		}

		// Change t,k to 1
		t = 1;
		k = 1;

#ifdef VERBOSE_CH_MULTI
		std::cerr << "polys = \n";
		for (unsigned int i = 0; i < polys.size(); ++i) {
			for (unsigned int j = 0; j < m; ++j) {
				std::cerr << "\t" << j << "\t" << polys[i].phis[j] << "\n";
			}
			std::cerr << "\n";
		}
#endif

#ifdef TIME_PARTS
		tsout << "\n";

		char * testOutfile = getenv("TIME_PARTS_OUTFILE");
		if (testOutfile) {
			std::ofstream timefile;
			timefile.open (testOutfile, ios::app);
			timefile << tsout.str();
			timefile.close();
		} else {
			std::cerr << "Time Parts: " << tsout.str();
		}
#endif

		return polys;
	}
}
#endif

template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<RecoveryPolyMulti<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_ch_tk1(
		unsigned int ell, unsigned int h, unsigned int m,
	    const vector<unsigned short> &goodservers,
		const vec_F &indices, const vector<vec_F> &shares) {

#ifdef VERBOSE_CH_TK1
	std::cerr << "Starting Cohn-Heninger Linear Multi-Polynomial algorithm...\n";
#endif

    // Get number of responses
	unsigned int n = goodservers.size();

#ifdef VERBOSE_CH_TK1
	std::cerr << "n = " << n << "\n";
	std::cerr << "ell = " << ell << "\n";
	std::cerr << "h = " << h << "\n";
	std::cerr << "m = " << m << "\n";
#endif

#ifdef TIME_PARTS
	std::stringstream tsout;
	tsout << n << "," << ell << "," << h << ",";
#endif

	// Create N(z)
	FX N;
	SetCoeff(N, 0, 1);
	for (unsigned int i = 0; i < n; ++i) {
		FX newTerm;
		SetCoeff(newTerm, 0, -(indices[goodservers[i]]));
		SetCoeff(newTerm, 1);
		N *= newTerm;
	}

#ifdef VERBOSE_CH_TK1
	std::cerr << "N(z) = " << N << " (" << deg(N) << ")\n";
#endif

	// Create L[i](z), the interpolation of points in shares[i]
	vector<FX> L;
	vec_F goodIndices(INIT_SIZE, n);
	for (unsigned int j = 0; j < n; ++j) {
		goodIndices[j] = indices[goodservers[j]];
	}
	for (unsigned int i = 0; i < m; ++i) {
		vec_F goodShares(INIT_SIZE, n);
		for (unsigned int j = 0; j < n; ++j) {
			goodShares[j] = shares[i][goodservers[j]];
		}
		FX L_i;
		interpolate(L_i, goodIndices, goodShares);
		L.push_back(L_i);
	}

#ifdef VERBOSE_CH_TK1
	for (unsigned int i = 0; i < m; ++i) {
		std::cerr << "L[" << i << "](z) = " << L[i] << " (" << deg(L[i]) << ")\n";
	}
#endif

	// Create lattice basis
	vector<vector<FX> > lattice(m+1, vector<FX>(m+1));
	lattice[0][0] = N;
	FX z_toell(ell, 1);
	for (unsigned int i = 0; i < m; ++i) {
		lattice[i+1][0] = -(L[i]);
		lattice[i+1][i+1] = z_toell;
	}

#ifdef VERBOSE_CH_TK1
	std::cerr << "lattice =\n";
	for (unsigned int i = 0; i <= m; ++i) {
		for (unsigned int j = 0; j <= m; ++j) {
			std::cerr << "\t" << i << "," << j << ":\t" << lattice[i][j] << "\n";
		}
	}
#endif

	vector<vector<FX> > aux (m+1, vector<FX>());

#ifdef TIME_PARTS
	struct timeval st, et;
	gettimeofday(&st, NULL);
#endif

	vector<vector<FX> > reducedLattice = reduce_lattice_MS(lattice, aux);

#ifdef TIME_PARTS
	gettimeofday(&et, NULL);
	unsigned long long elapsedus = ((unsigned long long)(et.tv_sec -
			st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
	tsout << "lattice reduction," << elapsedus << ",";
#endif

#ifdef VERBOSE_CH_TK1
	std::cerr << "reducedLattice =\n";
	for (unsigned int i = 0; i < m+1; ++i) {
		for (unsigned int j = 0; j < m+1; ++j) {
			std::cerr << "\t" << i << "," << j << ":\t" << reducedLattice[i][j] << "\n";
		}
	}
	std::cerr << "reducedLattice degrees:\n";
	for (unsigned int i = 0; i < m+1; ++i) {
		std::cerr << "\t";
		for (unsigned int j = 0; j < m+1; ++j) {
			std::cerr << deg(reducedLattice[i][j]) << "\t";
		}
		std::cerr << "\n";
	}
#endif

	vector<unsigned int> smallRows;
	for (unsigned int i = 0; i < m+1; ++i) {
		int maxdeg = -2;
		for (unsigned int j = 0; j < m+1; ++j) {
			if (deg(reducedLattice[i][j]) > maxdeg) {
				maxdeg = deg(reducedLattice[i][j]);
			}
 		}
//		std::cerr << "i = " << i << ", maxdeg = " << maxdeg << ", h = " << h << "\n";
		if (maxdeg < (int)(h)) {
			smallRows.push_back(i);
		}
	}

#ifdef VERBOSE_CH_TK1
	std::cerr << "smallRows =";
	vector<unsigned int>::iterator sriter;
	for (sriter = smallRows.begin(); sriter != smallRows.end(); ++sriter) {
		std::cerr << " " << *sriter;
	}
	std::cerr << "\n";
#endif

	if (smallRows.size() < m) {
#ifdef VERBOSE_CH_TK1
        std::cerr << "FAIL: Not enough small rows\n";
#endif
        return vector<RecoveryPolyMulti<FX> >();  // FAIL
    }

    vector<vector<FX> > A;
    vector<FX> b;
    vector<unsigned int>::iterator smallIter = smallRows.begin();
    for (unsigned int i = 0; smallIter != smallRows.end() && i < m; ++smallIter, ++i) {
        vector<FX> currRow;
        for (unsigned int j = 1; j <= m; ++j) {
            currRow.push_back(RightShift(reducedLattice[*smallIter][j], ell));
        }
        A.push_back(currRow);
        b.push_back(-(reducedLattice[*smallIter][0]));
    }

#ifdef VERBOSE_CH_TK1
    std::cerr << "A =\n";
    for (unsigned int i = 0; i < m; ++i) {
        for (unsigned int j = 0; j < m; ++j) {
            std::cerr << "\t" << i << "," << j << ":\t" << A[i][j] << "\n";
        }
    }
    std::cerr << "b =\n";
    for (unsigned int i = 0; i < m; ++i) {
        std::cerr << "\t" << i << ":\t" << b[i] << "\n";
    }
#endif

#ifdef TIME_PARTS
    gettimeofday(&st, NULL);
#endif

    //vector<FX> solution = solve_linsys_FX(A, b, t+1);
    vector<FX> solution = solve_linsys_FX(A, b);

#ifdef TIME_PARTS
    gettimeofday(&et, NULL);
    unsigned long long elapsedus = ((unsigned long long)(et.tv_sec -
            st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
    tsout << "solving linear system," << elapsedus << ",";
#endif

    if (solution.size() == 0) {
#ifdef VERBOSE_CH_TK1
        std::cerr << "FAIL: solve_linsys_FX failed\n";
#endif
        return vector<RecoveryPolyMulti<FX> >();  // Fail
    }

#ifdef VERBOSE_CH_TK1
    std::cerr << "solution =\n";
    for (unsigned int i = 0; i < m; ++i) {
        std::cerr << "\t" << i << ":\t" << solution[i] << "\n";
    }
#endif

    // For each polynomial in the solution, check how many input points it agrees
    // with.  If it's at least h, add it to the list of polys to return.
    vector< RecoveryPolyMulti<FX> > polys;
    vector<unsigned short>::iterator gooditer;
    vector<unsigned short> vecagree = goodservers;
    for (unsigned int i = 0; i < m; ++i) {
        for (gooditer = vecagree.begin(); gooditer != vecagree.end(); ) {
            F phival;
            eval(phival, solution[i], indices[*gooditer]);
            if (phival == shares[i][*gooditer]) {
                ++gooditer;
            } else {
                gooditer = vecagree.erase(gooditer);
            }
        }
    }
    if (vecagree.size() >= h) {
        RecoveryPolyMulti<FX> n(vecagree, solution);
        polys.push_back(n);
    }

#ifdef VERBOSE_CH_TK1
    std::cerr << "polys = \n";
    for (unsigned int i = 0; i < polys.size(); ++i) {
        for (unsigned int j = 0; j < m; ++j) {
            std::cerr << "\t" << j << "\t" << polys[i].phis[j] << "\n";
        }
        std::cerr << "\n";
    }
#endif

#ifdef TIME_PARTS
    tsout << "\n";

    char * testOutfile = getenv("TIME_PARTS_OUTFILE");
    if (testOutfile) {
        std::ofstream timefile;
        timefile.open (testOutfile, ios::app);
        timefile << tsout.str();
        timefile.close();
    } else {
        std::cerr << "Time Parts: " << tsout.str();
    }
#endif

    return polys;
}

/*
// Comparison functions used to let polynomials be keys for maps
bool cmp (const GF2EX& a, const GF2EX& b) { 
    if (deg(a) != deg(b)) {
        return deg(a) < deg(b);
    }    
    for (int i = deg(a); i >= 0; --i) {
        unsigned char byte_a, byte_b;
        BytesFromGF2X(&byte_a, rep(coeff(a, i)), 1); 
        BytesFromGF2X(&byte_b, rep(coeff(b, i)), 1); 
        if (byte_a != byte_b) {
            return byte_a < byte_b;
        }    
    }    
    return false;
}

bool cmp (const ZZ_pX& a, const ZZ_pX& b) { 
    if (deg(a) != deg(b)) {
        return deg(a) < deg(b);
    }    
    for (int i = deg(a); i >= 0; --i) {
        ZZ rep_a = rep(coeff(a, i)); 
        ZZ rep_b = rep(coeff(b, i)); 
        if (rep_a != rep_b) {
            return rep_a < rep_b;
        }    
    }    
    return false;
}
*/

template<class F, class vec_F, class FX, class FXY, class mat_F>
class cmpPolys {
public:
    bool operator() (const FX& a, const FX& b) {
        if (a == b) {
            return false;
        }
        return cmp(a, b);
    }
};

template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<RecoveryPoly<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_dp (
        unsigned int ell, unsigned int h,
        const vector<unsigned short> &goodservers,
        const vec_F &indices, const vec_F &shares,
        const DPType dptype, unsigned int gord)
{
    unsigned int n = goodservers.size();

    map<FX, std::set<unsigned short>, cmpPolys<F,vec_F,FX,FXY,mat_F> > soln_map;
    vector<unsigned short> assume_servers, sub_goodservers;

#ifdef VERBOSE_DP
    std::cerr << "\tindices = ";
    for (unsigned int i = 0; i < n; ++i) {
        std::cerr << indices[goodservers[i]] << " ";
    }
    std::cerr << "\n";
    std::cerr << "\tshares = ";
    for (unsigned int i = 0; i < n; ++i) {
        std::cerr << shares[goodservers[i]] << " ";
    }
    std::cerr << "\n";
#endif
    if (dptype == ASSUME_CORRECT) {
#ifdef VERBOSE_DP
        std::cerr << "ASSUME CORRECT: g = " << gord << "\n";
#endif
        assume_servers.insert(assume_servers.begin(), goodservers.begin(), goodservers.begin() + n - h + gord);
        subset_iterator subIter = subset_iterator(assume_servers, gord);
        while( ! subIter.atend() ) {
#ifdef VERBOSE_DP
            std::cerr << "\t*subIter = ";
            for (unsigned int i = 0; i < gord; ++i) {
                std::cerr << (*subIter)[i] << " ";
            }
            std::cerr << "\n";
#endif
            vector<unsigned short>::const_iterator siIter = (*subIter).begin();
            vector<unsigned short>::const_iterator gsIter = goodservers.begin();
            for (; gsIter != goodservers.end(); ++gsIter) {
                if (siIter != (*subIter).end() && *gsIter == *siIter) {
                    ++siIter;
                } else {
                    sub_goodservers.push_back(*gsIter);
                }
            }
#ifdef VERBOSE_DP
            std::cerr << "\tsub_goodservers = ";
            for (unsigned int i = 0; i < n - gord; ++i) {
                std::cerr << sub_goodservers[i] << " ";
            }
            std::cerr << "\n";
#endif

            FX L, B(0, 1);
            unsigned short numagree, numdisagree;
            vector<unsigned short> vecagree;
            test_interpolate(gord - 1, shares, indices, *subIter, goodservers, numagree, numdisagree, vecagree, L);
            for (unsigned int i = 0; i < gord; ++i) {
                FX mult(0, -(indices[(*subIter)[i]]));
                SetCoeff(mult, 1, 1);
                B *= mult;
            }

#ifdef VERBOSE_DP
            std::cerr << "\tL(z) = " << L << "\n";
            std::cerr << "\tB(z) = " << B << "\n";
#endif

            vec_F new_shares;
            new_shares.SetLength(shares.length());
            for (unsigned int i = 0; i < sub_goodservers.size(); ++i) {
                new_shares[sub_goodservers[i]] = (shares[sub_goodservers[i]] - eval(L, indices[sub_goodservers[i]])) / eval(B, indices[sub_goodservers[i]]);
            }

#ifdef VERBOSE_DP
            std::cerr << "\tnew_shares =";
            for (int i = 0; i < new_shares.length(); ++i) {
                std::cerr << " " << new_shares[i];
            }
            std::cerr << "\n";
#endif

            TestType testType = BEST;
            vector<RecoveryPoly<FX> > sub_polys = findpolys(ell - gord, h - gord, sub_goodservers, indices, new_shares, testType);

#ifdef VERBOSE_DP
            std::cerr << "\tsub_polys =\n";
            for (unsigned int j = 0; j < sub_polys.size(); ++j) {
                std::cerr << "\t\t" << sub_polys[j].phi << ", { ";
                for (unsigned int i = 0; i < sub_polys[j].G.size(); ++i) {
                    std::cerr << sub_polys[j].G[i] << " ";
                }
                std::cerr << "}\n";
            }
            std::cerr << "\tpolys +=\n";
#endif
            typename vector<RecoveryPoly<FX> >::iterator rpIter;
            for (rpIter = sub_polys.begin(); rpIter != sub_polys.end(); ++rpIter) {
                soln_map[rpIter->phi * B + L].insert(rpIter->G.begin(), rpIter->G.end());
            }

            sub_goodservers.clear();
            ++subIter;
        }

    } else if (dptype == ASSUME_WRONG) {
#ifdef VERBOSE_DP
        std::cerr << "ASSUME WRONG: g = " << gord << "\n";
#endif
        assume_servers.insert(assume_servers.begin(), goodservers.begin(), goodservers.begin() + h + gord);
        subset_iterator subIter = subset_iterator(assume_servers, gord);
        while( ! subIter.atend() ) {
#ifdef VERBOSE_DP
            std::cerr << "\t*subIter = ";
            for (unsigned int i = 0; i < gord; ++i) {
                std::cerr << (*subIter)[i] << " ";
            }
            std::cerr << "\n";
#endif
            vector<unsigned short>::const_iterator siIter = (*subIter).begin();
            vector<unsigned short>::const_iterator gsIter = goodservers.begin();
            for (; gsIter != goodservers.end(); ++gsIter) {
                if (siIter != (*subIter).end() && *gsIter == *siIter) {
                    ++siIter;
                } else {
                    sub_goodservers.push_back(*gsIter);
                }
            }
#ifdef VERBOSE_DP
            std::cerr << "\tsub_goodservers = ";
            for (unsigned int i = 0; i < n - gord; ++i) {
                std::cerr << sub_goodservers[i] << " ";
            }
            std::cerr << "\n";
#endif

            TestType testType = BEST;
            vector<RecoveryPoly<FX> > sub_polys = findpolys(ell, h, sub_goodservers, indices, shares, testType);

#ifdef VERBOSE_DP
            std::cerr << "\tsub_polys =\n";
            for (unsigned int j = 0; j < sub_polys.size(); ++j) {
                std::cerr << "\t\t" << sub_polys[j].phi << ", { ";
                for (unsigned int i = 0; i < sub_polys[j].G.size(); ++i) {
                    std::cerr << sub_polys[j].G[i] << " ";
                }
                std::cerr << "}\n";
            }
            std::cerr << "\tpolys +=\n";
#endif
            typename vector<RecoveryPoly<FX> >::iterator rpIter;
            for (rpIter = sub_polys.begin(); rpIter != sub_polys.end(); ++rpIter) {
                soln_map[rpIter->phi].insert(rpIter->G.begin(), rpIter->G.end());
            }

            sub_goodservers.clear();
            ++subIter;
        }

    } else if (dptype == ASSUME_SHARES) {
#ifdef VERBOSE_DP
        std::cerr << "ASSUME SHARES: d = " << gord << "\n";
#endif
        assume_servers.insert(assume_servers.begin(), goodservers.begin(), goodservers.begin() + gord);
        sub_goodservers.insert(sub_goodservers.begin(), goodservers.begin() + gord, goodservers.end());

        for (unsigned int r = 0; r <= gord && r <= ell + 1; ++r) {
            subset_iterator subIter = subset_iterator(assume_servers, r);
            while( ! subIter.atend() ) {
#ifdef VERBOSE_DP
                std::cerr << "\t*subIter = ";
                for (unsigned int i = 0; i < r; ++i) {
                    std::cerr << (*subIter)[i] << " ";
                }
                std::cerr << "\n";
#endif

                FX L, B(0, 1);
                if (r != 0) {
                    unsigned short numagree, numdisagree;
                    vector<unsigned short> vecagree;
                    test_interpolate(r - 1, shares, indices, *subIter, goodservers, numagree, numdisagree, vecagree, L);
                    for (unsigned int i = 0; i < r; ++i) {
                        FX mult(0, -(indices[(*subIter)[i]]));
                        SetCoeff(mult, 1, 1);
                        B *= mult;
                    }
                } else {
                    L = FX::zero();
                }

#ifdef VERBOSE_DP
                std::cerr << "\tL(z) = " << L << "\n";
                std::cerr << "\tB(z) = " << B << "\n";
#endif

                vec_F new_shares;
                new_shares.SetLength(shares.length());
                for (unsigned int i = 0; i < sub_goodservers.size(); ++i) {
                    new_shares[sub_goodservers[i]] = (shares[sub_goodservers[i]] - eval(L, indices[sub_goodservers[i]])) / eval(B, indices[sub_goodservers[i]]);
                }

#ifdef VERBOSE_DP
                std::cerr << "\tnew_shares =";
                for (int i = 0; i < new_shares.length(); ++i) {
                    std::cerr << " " << new_shares[i];
                }
                std::cerr << "\n";
#endif

                TestType testType = BEST;
                vector<RecoveryPoly<FX> > sub_polys = findpolys(ell - r, h - r, sub_goodservers, indices, new_shares, testType);

#ifdef VERBOSE_DP
                std::cerr << "\tsub_polys =\n";
                for (unsigned int j = 0; j < sub_polys.size(); ++j) {
                    std::cerr << "\t\t" << sub_polys[j].phi << ", { ";
                    for (unsigned int i = 0; i < sub_polys[j].G.size(); ++i) {
                        std::cerr << sub_polys[j].G[i] << " ";
                    }
                    std::cerr << "}\n";
                }
                std::cerr << "\tpolys +=\n";
#endif
                typename vector<RecoveryPoly<FX> >::iterator rpIter;
                for (rpIter = sub_polys.begin(); rpIter != sub_polys.end(); ++rpIter) {
                    soln_map[rpIter->phi * B + L].insert(rpIter->G.begin(), rpIter->G.end());
                }

                ++subIter;
            }
        }
        sub_goodservers.clear();

    } else {
#ifdef VERBOSE_DP
        std::cerr << "FAIL: Invalid DPType\n";
#endif
        return vector<RecoveryPoly<FX> >();  // FAIL
    }

    vector<RecoveryPoly<FX> > polys;
    typename map<FX, std::set<unsigned short>, cmpPolys<F,vec_F,FX,FXY,mat_F> >::iterator soln_iter;
    for (soln_iter = soln_map.begin(); soln_iter != soln_map.end(); ++soln_iter) {
        vector<unsigned short> new_G (soln_iter->second.begin(), soln_iter->second.end());
        polys.push_back(RecoveryPoly<FX>(new_G, soln_iter->first));
    }
    return polys;
}


template<class F, class vec_F, class FX, class FXY, class mat_F>
vector< RecoveryPoly<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_best (
    unsigned int n, int ell,
	unsigned int h, const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares, TestType &testType, DPType &dpType,
    int &gord)
{
/*
    char buffer[64];
    char testTypeStr[16], dpTypeStr[16];

    std::stringstream cmd;
    cmd << "grep --color=never '^" << n << "," << ell << "," << h << ",' best.csv | cut --output-delimiter=' ' -d',' -f4,5,6";
#ifdef VERBOSE_FINDPOLYS
    std::cerr << "cmd = " << cmd.str() << "\n";
#endif

    FILE * f = popen(cmd.str().c_str(), "r");
    if (f == 0) {
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: Could not read from best.csv\n";
#endif
        return vector<RecoveryPoly<FX> >();  // FAIL
    }
    if (fgets(buffer, sizeof(buffer), f) != buffer) {
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: Could not get data from command\n";
#endif
        fclose(f);
        return vector<RecoveryPoly<FX> >();  // FAIL
    }
    fclose(f);

    int ret = 0;
    ret = sscanf(buffer, "%s %s %d", testTypeStr, dpTypeStr, &gord); 
    if (ret != 1 && ret != 3) {
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: Could not parse data from best.csv\n";
#endif
        return vector<RecoveryPoly<FX> >();  // FAIL
    }
    for (unsigned int i = 2; i < MAX_TESTTYPE; ++i) {
        if (testTypeStr == testTypeStrings[i]) {
            testType = (TestType)(i);
            break;
        }
    }
    
    if (testType == DP) {
        for (unsigned int i = 1; i < MAX_DPTYPE; ++i) {
            if (dpTypeStr == DPTypeStrings[i]) {
                dpType = (DPType)(i);
                break;
            }
        }
    } else {
        dpType = UNDEFINED_DPTYPE;
        gord = 1;
    }
    */

    portfolioChoice(n, ell, h, testType, dpType, gord);

#ifdef VERBOSE_FINDPOLYS
    std::cerr << "The best algoritm is:\n";
    std::cerr << "    testType = " << testTypeStrings[testType] << "\n";
    std::cerr << "    dpType = " << DPTypeStrings[dpType] << "\n";
    std::cerr << "    gord = " << gord << "\n";
#endif

    DPType usedDPType = dpType;
    vector<RecoveryPoly<FX> > polys = findpolys(ell, h, goodservers, indices, shares, testType, usedDPType, gord);

    return polys;
}


template<class F, class vec_F, class FX, class FXY, class mat_F>
vector< RecoveryPoly<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys (int k,
	unsigned int t, const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares, TestType &testType, DPType &dpType, int &gord)
{
    unsigned int n = goodservers.size();

#ifdef VERBOSE_DP
    std::cerr << "RUNNING findpolys WITH (" << testTypeStrings[testType] << ", " << n << ", " << k << ", " << t << ")\n";
#endif
    
    vector<RecoveryPoly<FX> > polys;
    switch(testType) {
    case BEST:
        polys = findpolys_best(n, k, t, goodservers, indices, shares, testType, dpType, gord);
        break;
    case BRUTE:
        polys = findpolys_brute(k, t, goodservers, indices, shares);
        break;
    case KOTTER:
    case CH_MS:
        polys = findpolys_gs(k, t, goodservers, indices, shares, testType);
        break;
    case BW:
        polys = findpolys_bw(k, t, goodservers, indices, shares);
        break;
    case DP:
        polys = findpolys_dp(k, t, goodservers, indices, shares, dpType, gord);
        break;
    case CH_MULTI:
    case CH_TK1:
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: Invalid single poly test type\n";
#endif
        return vector<RecoveryPoly<FX> >();  // FAIL
        break;
    case UNDEFINED:
    case UNKNOWN:
    case MAX_TESTTYPE:
    default:
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: Not a valid test type\n";
#endif
        return vector<RecoveryPoly<FX> >();  // FAIL
        break;
    }

    (void)n;
    return polys;
}

template<class F, class vec_F, class FX, class FXY, class mat_F>
vector< RecoveryPolyMulti<FX> > RSDecoder<F,vec_F,FX,FXY,mat_F>::findpolys_multi (unsigned int k,
	unsigned int t, const vector<unsigned short> &goodservers,
	const vec_F &indices, const vector<vec_F> &shares, TestType testType)
{
	unsigned int n = goodservers.size();
	unsigned int m = shares.size();

    vector<RecoveryPolyMulti<FX> > multipolys;
    vector<RecoveryPoly<FX> > polys;
    switch(testType) {

    case BEST:
    case BRUTE:
    case KOTTER:
    case CH_MS:
    case BW:
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: Invalid multi poly test type\n";
#endif
        return vector<RecoveryPolyMulti<FX> >();  // FAIL
        break;

    case CH_MULTI:
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "FAIL: The test type CH_MULTI is not implemented\n";
#endif
        return vector<RecoveryPolyMulti<FX> >();  // FAIL
        break;

    case CH_TK1:
        multipolys = findpolys_ch_tk1(k, t, m, goodservers, indices, shares);
        break;

    case UNDEFINED:
    case UNKNOWN:
    case MAX_TESTTYPE:
    default:
#ifdef VERBOSE_FINDPOLYS
        std::cerr << "Invalid test type!\n";
#endif
		return vector<RecoveryPolyMulti<FX> >();  // Fail
    }

    (void)n;
	return multipolys;
}


template <class F, class FX> 
void addResult (vector<DecoderResult<F> > &results, nservers_t h,
	dbsize_t word_number, vector<RecoveryPoly<FX> > polys)
{
    typename vector<RecoveryPoly<FX> >::iterator poly_iter;
    if (results.empty()) {
        for (poly_iter = polys.begin(); poly_iter != polys.end(); ++poly_iter) {
            if (poly_iter->G.size() >= h) {
                F wz; 
                eval(wz, poly_iter->phi, F::zero());
                map<dbsize_t, F> new_recovered;
                new_recovered[word_number] = wz;
                results.push_back(DecoderResult<F>(poly_iter->G, new_recovered));
            }
        }
        return;
    }

    if (results[0].recovered.find(word_number) != results[0].recovered.end()) {
#ifdef VERBOSE_RECOVER
        std::cerr << "This one is already done!\n";
#endif
        return;
    }

    ssize_t result_size = results.size();
    for (ssize_t idx=0; idx < result_size; ++idx) {
	vector<nservers_t> orig_G = results[idx].G;

        ssize_t found = -1;
        for (poly_iter = polys.begin(); poly_iter != polys.end(); ++poly_iter) {
            vector<nservers_t> intersection = intersect(poly_iter->G, orig_G);
            if (intersection.size() >= h) {
                F wz; 
                eval(wz, poly_iter->phi, F::zero());
                if (found == -1) {
		    results[idx].G = intersection;
                    results[idx].recovered[word_number] = wz; 
                    found = idx;
                } else {
                    results.push_back(DecoderResult<F>(intersection,
					results[idx].recovered));
                    results[results.size()-1].recovered[word_number] = wz;
                }
            }   
        }
        if (found == -1) {
	    results.erase(results.begin()+idx);
            --idx;
            --result_size;
        }
    }   
}


template<class F, class vec_F, class FX, class FXY, class mat_F>
bool RSDecoder<F,vec_F,FX,FXY,mat_F>::Recover (dbsize_t bytes_per_word,
	nservers_t t, nservers_t h, const vector<nservers_t> &goodservers,
	const vector<vector<vec_F> > &values, const vec_F &indices,
	vector<vector<DecoderResult<F> > > &results,
	vector<std::set<dbsize_t> > &decoded,
	const nqueries_t multi_only)
{
#ifdef VERBOSE_RECOVER
    std::cerr << "h = " << h << "\n";
#endif
    TestType ttype;
    nservers_t k = goodservers.size();

    // Keep track of undecoded words by query
    map<nqueries_t, std::set<dbsize_t> > undecoded;

    nqueries_t num_queries = values.size();
    for (nqueries_t q = 0; q < num_queries; ++q) {
	undecoded[q] = std::set<dbsize_t>();
	dbsize_t num_words = values[q].size();
	for (dbsize_t i = 0; i < num_words; ++i) {
#ifdef VERBOSE_RECOVER
	    std::cerr << "Decoding (" << q << ", " << i << "):\n";
	    std::cerr << "    (k, t, h) = (" << k << ", " << t << ", " << h <<
		")\n";
#endif
	    if (decoded[q].find(i) != decoded[q].end()) {
#ifdef VERBOSE_RECOVER
		std::cerr << "    Already decoded.\n";
#endif
		continue;
	    }

	    if (q >= multi_only) {
		// Try EasyRecover
		bool easy_result = EasyRecover(bytes_per_word, t, h, results[q], i,
			values[q][i], indices);
		if (easy_result) {
#ifdef VERBOSE_RECOVER
		    std::cerr << "    EASY: PASS\n";
#endif
		    decoded[q].insert(i);
		    continue;
		}
#ifdef VERBOSE_RECOVER
		std::cerr << "    EASY: FAIL\n";
#endif

		// Try Berlekamp-Welch
		ttype = BW;
		vector<RecoveryPoly<FX> > bw_result = findpolys(t, h, goodservers, indices, values[q][i], ttype);
		if (bw_result.size() > 0) {
#ifdef VERBOSE_RECOVER
		    std::cerr << "    BW: PASS\n";
#endif
		    addResult<F,FX>(results[q], h, i, bw_result);
		    decoded[q].insert(i);
		    continue;
		}
#ifdef VERBOSE_RECOVER
		std::cerr << "    BW: FAIL\n";
#endif

		// Try findpolys_best
		ttype = BEST;
		// Only run best in single-polynomial decoding realm or when tk1
		// can't give an answer (h == t+1)
		if (h > sqrt(k*t) || h == t+1) {
		    vector<RecoveryPoly<FX> > best_result = findpolys(t, h, goodservers, indices, values[q][i], ttype);
		    if (best_result.size() > 0) {
#ifdef VERBOSE_RECOVER
			std::cerr << "    BEST: PASS (#results: " << best_result.size() << ")\n";
			for (size_t iw=0;iw<best_result.size(); ++iw) {
			    std::cerr << "        " << best_result[iw].phi << "\n";
			}
#endif
			addResult<F,FX>(results[q], h, i, best_result);
			decoded[q].insert(i);
			continue;
		    }
#ifdef VERBOSE_RECOVER
		    std::cerr << "    BEST: FAIL\n";
#endif
		}
	    }

	    // Save for multi
#ifdef VERBOSE_RECOVER
	    std::cerr << "    Saving for multi...\n";
#endif
	    undecoded[q].insert(i);
	}
	if (undecoded[q].empty()) {
#ifdef VERBOSE_RECOVER
	    std::cerr << "Query " << q << " is done\n";
#endif
	    undecoded.erase(q);
	}
    }

    // Check if done
    if (undecoded.empty()) {
#ifdef VERBOSE_RECOVER
        std::cerr << "All queries are done\n";
#endif
        return true;
    }

    // If there are enough for tk1, do tk1
    double queries_needed = (double)(k - h) / (double)(h - t - 1);
#ifdef VERBOSE_RECOVER
    std::cerr << "queries_needed = " << queries_needed << "\n";
#endif

    // Do multi while we have enough queries
    vector<nqueries_t> ud_queries;
    vector<typename std::set<dbsize_t>::iterator> ud_iters;
    map<nqueries_t, std::set<dbsize_t> >::iterator ud_iter;
    for (ud_iter = undecoded.begin(); ud_iter != undecoded.end(); ++ud_iter) {
	ud_queries.push_back(ud_iter->first);
	ud_iters.push_back(ud_iter->second.begin());
    }
    while (ud_queries.size() >= queries_needed) {
        vector<vec_F> multi_values;
	for (nqueries_t i = 0; i < ud_queries.size(); ++i) {
	//for (ud_iter = undecoded.begin(); ud_iter != undecoded.end(); ++ud_iter) {
	    nqueries_t q = ud_queries[i];
	    dbsize_t w = *(ud_iters[i]);
            multi_values.push_back(values[q][w]);
        }
        ttype = CH_TK1;
        vector<RecoveryPolyMulti<FX> > multi_result = findpolys_multi(t, h, goodservers, indices, multi_values, ttype);

        if (multi_result.size() == 0) {
#ifdef VERBOSE_RECOVER
            std::cerr << "    TK1: FAIL\n";
#endif
	    for (nqueries_t i = 0; i < ud_queries.size(); ++i) {
		++(ud_iters[i]);
	    }
#endif
        } else {
#ifdef VERBOSE_RECOVER
	    std::cerr << "    TK1: PASS\n";
#endif

	    for (nqueries_t i = 0; i < ud_queries.size(); ++i) {
	    //for (ud_iter = undecoded.begin(), i = 0; ud_iter != undecoded.end(); ++ud_iter, ++i) {
		nqueries_t q = ud_queries[i];
		dbsize_t w = *(ud_iters[i]);
		++(ud_iters[i]);
		vector<RecoveryPoly<FX> > ith_result;
		typename vector<RecoveryPolyMulti<FX> >::iterator iter;
		for (iter = multi_result.begin(); iter != multi_result.end();
			++iter) {
		    ith_result.push_back(RecoveryPoly<FX>(iter->G, iter->phis[i]));
		}
		addResult<F,FX>(results[q], h, w, ith_result);
		decoded[q].insert(w);
		undecoded[q].erase(w);
	    }
	}

	for (nqueries_t i = 0; i < ud_queries.size(); ++i) {
	//for (ud_iter = undecoded.begin(); ud_iter != undecoded.end(); ++ud_iter) {
	    nqueries_t q = ud_queries[i];
	    if (ud_iters[i] == undecoded[q].end()) {
	    //if (w >= ud_iter->second.size()) {
#ifdef VERBOSE_RECOVER
		if (undecoded[q].empty()) {
		    undecoded.erase(q);
		    std::cerr << "Query " << q << " is done\n";
		} else {
		    std::cerr << "Query " << q << " was not fully decoded\n";
		}
#endif
		ud_queries.erase(ud_queries.begin() + i);
		ud_iters.erase(ud_iters.begin() + i);
		--i;
	    }
	}
    }

    // Check if done
    if (undecoded.empty()) {
#ifdef VERBOSE_RECOVER
        std::cerr << "All queries are done\n";
#endif
        return true;
    }
#ifdef VERBOSE_RECOVER
    std::cerr << "NOT all queries were decoded\n";
#endif
    return false;
}


// Depth-first search on the tree of coefficients; used by the
// Roth-Ruckenstein algorithm.
template<class F, class vec_F, class FX, class FXY, class mat_F>
void RSDecoder<F,vec_F,FX,FXY,mat_F>::dfs(vector<FX> &res,
		int u, 
		vector<int> &pi, 
		vector<int> &Deg,
		vector<F> &Coeff, 
		vector<FXY> &Q, 
		int &t, 
		int degreebound)
{
#ifdef TEST_RR
    cout << "\nVertex " << u << ": pi[" << u << "] = " << pi[u] <<
	", deg[" << u << "] = " << Deg[u] << ", Coeff[" << u <<
	"] = " << Coeff[u] << "\n";
    cout << "Q[" << u << "] = " << Q[u] << "\n";
#endif
    // if Q_u(x,0) == 0 then output y-root
    if ( IsZero(Q[u]) || IsZero(Q[u].rep[0])) { 
	// Output f_u[x]

	FX fux;
	while (Deg[u] >= 0) {
	    SetCoeff(fux, Deg[u], Coeff[u]);
	    u = pi[u];
	}
#ifdef TEST_RR
	cout << "Outputting " << fux << "\n";
#endif
	res.push_back(fux);
    } else if (degreebound < 0 || Deg[u] < degreebound) {
	// Construct Q_u(0,y)
	FX Qu0y;
	int degqu = deg(Q[u]);
	for (int d = 0; d <= degqu; ++d) {
	    SetCoeff(Qu0y, d, coeff(Q[u].rep[d], 0));
	}
      
#ifdef TEST_RR
	cout << "Q[" << u << "](0,y) = " << Qu0y << "\n";
#endif
	// Find its roots
	vec_F rootlist =
	    findroots_FX<F,vec_F,FX,
		typename DT::vec_FX,typename DT::vec_pair_FX_long>(Qu0y);
#ifdef TEST_RR
	cout << "Rootlist = " << rootlist << "\n";
#endif
	int numroots = rootlist.length();
	for (int r=0; r<numroots; ++r) {
	    // For each root a, recurse on <<Q[u](x, x*y + a)>>
	    // where <<Q>> is Q/x^k for the maximal k such that the
	    // division goes evenly.
	    int v = t;
	    ++t;
	    pi.push_back(u);
	    Deg.push_back(Deg[u]+1);
	    Coeff.push_back(rootlist[r]);
	    // mapit(Q(x,y), a) computes <<Q(x, x*y + a)>>
	    Q.push_back(mapit(Q[u], rootlist[r]));
	    dfs(res, v, pi, Deg, Coeff, Q, t, degreebound);
	}
    }
}

// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  This routine only works over fields F.
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<FX> RSDecoder<F,vec_F,FX,FXY,mat_F>::rr_findroots(
	const FXY &P, int degreebound)
{
    vector<int> pi;
    vector<int> Deg;
    vector<F> Coeff;
    vector<FXY> Q;
    vector<FX> res;
    int t = 1;

#ifdef TEST_RR
    std::cerr << "Now in findroots(" << P << ", " << degreebound << ")\n";
#endif

    FXY P0 = backShiftX(P, minX(P));
    pi.push_back(-1);
    Deg.push_back(-1);
    Coeff.push_back(F::zero());
    Q.push_back(P0);
    int u = 0; 

    if (degreebound < 0) {
	degreebound = degX(P0);
    }
    dfs(res, u, pi, Deg, Coeff, Q, t, degreebound);

    return res;
}


// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  This routine handles the case where F is the
// integers mod p1*p2 as well as fields.  (But it handles the ring case
// in a specialization of this templated method.) This routine may alsofromsage[*rowIter][i];
// return some spurious values.
template<class F, class vec_F, class FX, class FXY, class mat_F>
vector<FX> RSDecoder<F,vec_F,FX,FXY,mat_F>::findroots(const FXY &P, int degreebound)
{
//	std::cerr << "P = " << P << "\n";
    vector<FX> result = rr_findroots(P, degreebound);
//    for (unsigned int i = 0; i < result.size(); ++i) {
//    	std::cerr << "result[" << i << "] = " << result[i] << "\n";
//    }
    return result;
}

