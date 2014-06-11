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
#include <vector>
#include <sstream>
#include <string.h>
#include <NTL/ZZ.h>
#include "rsdecoder.h"

NTL_CLIENT


uint64_t hasseop = 0, kotter_usec = 0;

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
template<>
string RSDecoder_ZZ_p::append(const string &s, const ZZ_p &wz,
	unsigned int bytes_per_word)
{
    unsigned char * w = new unsigned char[bytes_per_word];
    BytesFromZZ(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    delete[] w;
    return r;
}

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
template<>
string RSDecoder_GF2E::append(const string &s, const GF2E &wz,
	unsigned int bytes_per_word)
{
    unsigned char * w = new unsigned char[bytes_per_word];
    BytesFromGF2X(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    delete[] w;
    return r;
}

// Project a polynomial created with a ZZ_pContext of k*p down to the
// current ZZ_pContext of p.
static ZZ_pXY project_down(const ZZ_pXY &P)
{
    ZZ_pXY newP;
    stringstream ss (stringstream::in | stringstream::out);
    ss << P;
    ss >> newP;
    return newP;
}

// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  The global ZZ_pContext should already be set to
// p1 * p2, where p2 is prime, and p1 is either prime or 1.
// This routine may also return some spurious values.
template<>
vector<ZZ_pX> RSDecoder_ZZ_p::findroots(
	const ZZ_pXY &P, int degreebound)
{
    // If we're already working mod a prime, just go ahead
    if (p1 == 1) {
	return rr_findroots(P, degreebound);
    }

    // We have to find the roots mod each prime separately, and combine
    // the results with the CRT if p1 > 1.
    ZZ_pBak pbak;
    pbak.save();

    vector<ZZ_pX> roots_p1, roots_p2, roots;

    ZZ_p::init(p1);
    ZZ_pXY P1 = project_down(P);
    roots_p1 = rr_findroots(P1, degreebound);

    ZZ_p::init(p2);
    ZZ_pXY P2 = project_down(P);
    roots_p2 = rr_findroots(P2, degreebound);

    pbak.restore();

    // Calcuate a1 and a2 s.t. a1 = (0, 1) mod (p1,p2) and
    // a2 = (1, 0) mod (p1, p2).
    ZZ_p a1 = to_ZZ_p(p1);
    ZZ_p a2 = to_ZZ_p(p2);
    // cerr << "p1 = " << p1 << "\np2 = " << p2 << "\np1 * p2 = " << p1 * p2 << "\n";
    a1 *= to_ZZ_p(InvMod(AddMod(p1, 0, p2), p2));
    // cerr << "a1 = " << a1 << "\n";
    a2 *= to_ZZ_p(InvMod(AddMod(p2, 0, p1), p1));
    // cerr << "a2 = " << a2 << "\n";
    
    // For each pair, use the CRT to combine them
    unsigned short num_p1 = roots_p1.size();
    unsigned short num_p2 = roots_p2.size();
    for (unsigned short i=0; i<num_p1; ++i) {
	for (unsigned short j=0; j<num_p2; ++j) {
	    ZZ_pX comb = roots_p1[i] * a2 + roots_p2[j] * a1;
	    roots.push_back(comb);
	}
    }

    return roots;
}


#ifdef TEST_FINDPOLYS

#include <sstream>
int main()
{
    ZZ modulus, one;
    modulus = 10007;
    one = 1;
    ZZ_p::init(modulus);
    RSDecoder_ZZ_p decoder(one, modulus);

    stringstream ss(stringstream::in | stringstream::out);

    vector<unsigned short> goodservers;
    goodservers.push_back(0);
    goodservers.push_back(1);
    goodservers.push_back(2);
    goodservers.push_back(3);
    goodservers.push_back(4);
    goodservers.push_back(5);
    vec_ZZ_p indices, shares;
    ss << "[ 1 2 3 4 5 6 ]";
    ss >> indices;
    ss << "[ 4 3 7507 1 2504 7 ]";
    ss >> shares;
    vector<RecoveryPoly<ZZ_pX> > ret = decoder.findpolys(3, 5, goodservers, indices, shares);
    for (vector<RecoveryPoly<ZZ_pX> >::const_iterator iter = ret.begin(); iter != ret.end(); ++iter) {
	cout << "{ " << iter->phi << ", [ ";
	for (vector<unsigned short>::const_iterator gi = iter->G.begin(); gi != iter->G.end(); ++gi) {
	    cout << *gi << " ";
	}
	cout << "] }\n";
    }
}

#endif

#ifdef TEST_RR

// Test the algorithm on the Example 15 from the McElice paper.
int main()
{
    ZZ modulus;
    modulus = 19;
    ZZ p1;
    p1 = 1;
    ZZ_p::init(modulus);
    RSDecoder_ZZ_p decoder(p1, modulus);

    stringstream ss(stringstream::in | stringstream::out);
    ZZ_pXY P;

    ss << "[[4 12 5 11 8 13] [14 14 9 16 8] [14 13 1] [2 11 1] [17]] ";
    ss >> P;

    vector<ZZ_pX> roots = decoder.findroots(P, 1);
    cout << "\nRoots found:\n";
    for (unsigned int i=0; i<roots.size(); ++i) {
	cout << roots[i] << endl;

    }

    cout << "\nExpected output (in some order):\n";
    cout << "[18 14]\n[14 16]\n[8 8]\n";
    return 0;
}
#endif

#ifdef TIME_FINDPOLYS

#include <iostream>
#include <sstream>
#include <math.h>
#include <sys/time.h>
#include <time.h>

template <class F, class vec_F, class FX, class FXY, class mat_F>
static RSDecoder<F,vec_F,FX,FXY,mat_F> do_init()
{
    RSDecoder<F,vec_F,FX,FXY,mat_F> decoder;
    return decoder;
}

template <>
RSDecoder_ZZ_p do_init()
{
    ZZ modulus, one;

    stringstream ss(stringstream::in | stringstream::out);
    // 128-bit modulus
#ifdef USE_W8
    ss << "257 ";
#elif defined USE_W16
    ss << "65537 ";
#elif defined USE_W32
    ss << "4294967311 ";
#else
    ss << "340282366920938463463374607431768211507 ";
#endif
    ss >> modulus;
    one = 1;

    ZZ_p::init(modulus);
    return RSDecoder_ZZ_p(one, modulus);
}

template <>
RSDecoder_GF2E do_init()
{
    // Initialize the GF2E modulus to the one used by AES
    GF2X AES_P;
#ifdef USE_GF24
    SetCoeff(AES_P, 4, 1);
    SetCoeff(AES_P, 1, 1);
    SetCoeff(AES_P, 0, 1);
#else
    SetCoeff(AES_P, 8, 1);
    SetCoeff(AES_P, 4, 1);
    SetCoeff(AES_P, 3, 1);
    SetCoeff(AES_P, 1, 1);
    SetCoeff(AES_P, 0, 1);
#if 0  // GF(2^16)
    SetCoeff(AES_P, 16, 1);
    SetCoeff(AES_P, 5, 1);
    SetCoeff(AES_P, 3, 1);
    SetCoeff(AES_P, 2, 1);
    SetCoeff(AES_P, 0, 1);
#endif
#endif

    GF2E::init(AES_P);
    GF2X::HexOutput = 1;

    return RSDecoder_GF2E();
}

ZZ binomial (ZZ n, ZZ k) {
	if (n -k > k) {
		k = n - k;
	}
	ZZ result(ZZ::zero() + 1);
	for (unsigned int i = 0; i < k; ++i) {
		result *= n - i;
		result /= i + 1;
	}
	return result;
}

// Keeps track of seed value for testing purposes.
ZZ seed_g;

template<class F, class vec_F, class FX, class FXY, class mat_F>
static void time_findpolys(TestType testType, int k, int t, int h, int s, int multi_t, int multi_k,
        DPType dpType, int gord)
{
    const char *brute = getenv("PIRC_BRUTE");
    if (brute == NULL) brute = "0";

    TestType origTestType = testType;

	std::stringstream ssout;

#ifdef USE_GF28
	ssout << "gf28,";
#elif defined USE_GF24
	ssout << "gf24,";
#elif defined USE_W8
	ssout << "w8,";
#elif defined USE_W16
	ssout << "w16,";
#elif defined USE_W32
	ssout << "w32,";
#else
	ssout << "w128,";
#endif
	ssout << testTypeStrings[testType] << "," << k << "," << t << "," << h << ",";

	// Test preprocessing
	bool dnr = false;
	int dnr_low_denominator = 9;
	unsigned int ch_multi_mins = 1, ch_multi_t = 1, ch_multi_k = 1,
			curr_t, curr_k, curr_m, max_m;
	const unsigned int add_to_mins = 10;
	ZZ min_runtime(ZZ::zero() - 1), curr_runtime, ineq;
	switch (testType) {
	case KOTTER:
		// For now, do not run kotter when t == 0

		if (t <= 0) {
			dnr = true;
			ssout << "infty,dnr_kotter_low_t,";
			break;
		}

        if (h <= sqrt(k * t)) {
            dnr = true;
            ssout << "infty,dnr_gs_inequality,";
            break;
        }

        break;
	case CH_MS:
		// Do not run kotter, ch_ms when h^2-kt <= dnr_denom
		if (h * h - k * t <= dnr_low_denominator) {
			dnr = true;
			ssout << "infty,dnr_low_denominator,";
			break;
		}

		// Do not run kotter, ch_ms when h <= sqrt(k*t)
		if (h <= sqrt(k * t)) {
			dnr = true;
			ssout << "infty,dnr_gs_inequality,";
			break;
		}

		// Do not run kotter, ch_ms when t < 0;
		if (t < 0) {
			dnr = true;
			ssout << "infty,dnr_invalid_arg,";
			break;
		}
		break;
	case BW:
		// Do not run bw when 2h <= k+t
		if (2 * h <= k + t) {
			dnr = true;
			ssout << "infty,dnr_bw_inequality,";
			break;
		}

		// Do not run bw when t < 0;
		if (t < 0) {
			dnr = true;
			ssout << "infty,dnr_invalid_arg,";
			break;
		}
		break;
	case CH_MULTI:
		// Do not run ch_multi when t < 0;
		if (t < 0) {
			dnr = true;
			ssout << "infty,dnr_invalid_arg,";
			break;
		}

		if (multi_t != 0 && multi_k != 0) {
			ch_multi_t = (unsigned int)multi_t;
			ch_multi_k = (unsigned int)multi_k;
			break;
		}
		// Choose m,t,k
		ch_multi_mins = (unsigned int)(log(double(k)/double(h)) / log(double(h)/double(t))) + 1;
		max_m = ch_multi_mins + add_to_mins;
		for (curr_t = 0;; ++curr_t) {
			if (max_m < ch_multi_mins) {
				break;
			}
			for (curr_m = ch_multi_mins; curr_m <= (unsigned int)max_m; ++curr_m) {
				// ceil of k approx.
				curr_k = (unsigned int)(pow(double(h)/double(k), 1/double(curr_m)) * curr_t) + 1;
				ineq = curr_m*t*binomial(to_ZZ(curr_m+curr_t),to_ZZ(curr_m+1))
						+ k*binomial(to_ZZ(curr_m+curr_k),to_ZZ(curr_m+1))
						- h*curr_k*binomial(to_ZZ(curr_m+curr_t),to_ZZ(curr_m))
						+ h*curr_k*(curr_m-1);
				curr_runtime = power(binomial(to_ZZ(curr_m+curr_t),to_ZZ(curr_m)),3)*((k-1)*curr_t)*((k-1)*curr_t-curr_k*h);
				if (ineq >= 0) {
					// sub one (pseudo-floor)
					--curr_k;
					ineq = curr_m*t*binomial(to_ZZ(curr_m+curr_t),to_ZZ(curr_m+1))
							+ k*binomial(to_ZZ(curr_m+curr_k),to_ZZ(curr_m+1))
							- h*curr_k*binomial(to_ZZ(curr_m+curr_t),to_ZZ(curr_m))
							+ h*curr_k*(curr_m-1);
					curr_runtime = power(binomial(to_ZZ(curr_m+curr_t),to_ZZ(curr_m)),3)*((k-1)*curr_t)*((k-1)*curr_t-curr_k*h);
					if (min_runtime != -1 && curr_runtime > min_runtime) {
						max_m = curr_m - 1;
						continue;
					}
					if (ineq >= 0) {
						continue;
					}
				}
				max_m = curr_m - 1;
				if (min_runtime == -1 || curr_runtime < min_runtime) {
					min_runtime = curr_runtime;
					s = curr_m;
					ch_multi_t = curr_t;
					ch_multi_k = curr_k;
				}
			}
		}
//		std::cerr << "Best case: m = " << s << ", t = " << ch_multi_t << ", k = " << ch_multi_k << "\n";
//		return;
		break;
	case CH_TK1:
		// Set s
        if (getenv("PIRC_M")) {
            s = atoi(getenv("PIRC_M"));
        } else {
    		s = ceil((double)(k - h) / (h - t - 1));
        }

		// Do not run ch_multi when t < 0;
		if (t < 0) {
			dnr = true;
			ssout << "infty,dnr_invalid_arg,";
			break;
		}
		break;
	case BRUTE:
		// Do not run brute when t < -1;
		if (t < -1) {
			dnr = true;
			ssout << "infty,dnr_invalid_arg,";
			break;
		}

		// Do not run if n > 20 and 1/4*n < h < 3/4*n
		if (k > 20 && h > ((double)k / 4) && h < (3 * (double)k / 4)) {
			dnr = true;
			ssout << "infty,dnr_big_brute,";
			break;
		}
		break;
    case DP:
        break;
	case BEST:
        // Do not run for h <= t + 1
        if (h <= t + 1) {
            dnr = true;
            ssout << "infty,dnr_best_invalid_case,";
            break;
        }
        break;
	case UNDEFINED:
	case UNKNOWN:
	case MAX_TESTTYPE:
	default:
		break;
	}

	if (!dnr) {
		RSDecoder<F, vec_F, FX, FXY, mat_F> decoder =
		do_init<F,vec_F,FX,FXY,mat_F>();

		// Construct s random polynomials of degree t
		vector<FX> randpolys;
		for (int i = 0; i < s; ++i) {
#if defined(USE_GF28) || defined(USE_GF24)
			randpolys.push_back(random_GF2EX(t+1));
#else
			randpolys.push_back(random_ZZ_pX(t+1));
#endif

#ifdef VERBOSE_FINDPOLYS
			std::cerr << "Original " << i << ": " << randpolys.back() << "\n";
#endif
		}

		vec_F indices;
		vector<vec_F> shares(s);

#if 0
		{
		F r = random_F();
		struct timeval st, et;
		gettimeofday(&st, NULL);
		for (int i = 0; i<1000000; ++i) {
		r += r;
		}
		gettimeofday(&et, NULL);
		unsigned long long elapsedus = ((unsigned long long)(et.tv_sec -
			st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
		cerr << "+: " << elapsedus << endl;
		gettimeofday(&st, NULL);
		for (int i = 0; i<1000000; ++i) {
		r *= r;
		}
		gettimeofday(&et, NULL);
		elapsedus = ((unsigned long long)(et.tv_sec -
			st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
		cerr << "*: " << elapsedus << endl;
		gettimeofday(&et, NULL);
		for (int i = 0; i<10; ++i) {
		r = power(r,rep(r));
		}
		gettimeofday(&et, NULL);
		elapsedus = ((unsigned long long)(et.tv_sec -
			st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);
		cerr << "^: " << elapsedus/10 << endl;
		exit(0);
		}
#endif

		indices.SetLength(k);
		for (int i = 0; i < s; ++i) {
			shares[i].SetLength(k);
		}

		// Construct the indices and shares
		for (int i=0; i<k; ++i) {
#if defined(USE_GF28) || defined(USE_GF24)
			unsigned char b = i+1;
			conv(indices[i], GF2XFromBytes(&b, 1));
#else
			indices[i] = i+1;
#endif
			for (int j = 0; j < s; ++j) {
				eval(shares[j][i], randpolys[j], indices[i]);
			}
		}

		// Pick a random subset of them to be wrong
		vector<unsigned short> allservers, wrongservers;

		for (int i=0; i<k; ++i) {
		allservers.push_back(i);
		}

		random_subset(allservers, wrongservers, k-h);

		bool rand_when = getenv("USE_RAND") && strchr(getenv("USE_RAND"), 'w');
		bool rand_amount = getenv("USE_RAND") && strchr(getenv("USE_RAND"), 'a');
		bool rand_srv = getenv("USE_RAND") && strchr(getenv("USE_RAND"), 's');

#ifdef TEST_TK1_ADD1
        map<unsigned short, F> srvvals;
#endif

		for(vector<unsigned short>::iterator iter = wrongservers.begin();
			iter != wrongservers.end(); ++iter) {
			F srvval;
			if (rand_srv) {
				random(srvval);
#ifdef TEST_TK1_ADD1
                srvvals[*iter] = srvval;
#endif
			}
			for (int i = 0; i < s; ++i) {
				if (!rand_when || (rand() % 2)) {
					F randval;
					if (rand_amount) {
						random(randval);
						shares[i][*iter] += randval;
					} else if (rand_srv) {
                shares[i][*iter] += srvval * (i ? 1 : *iter);
					} else {
						shares[i][*iter] += (i ? 1 : *iter);
					}
				}
			}
		}

#ifdef TEST_TK1_ADD1
        std::cerr << "shares (First Time) =\n";
        for (unsigned int i = 0; i < s; ++i) {
            std::cerr << "\t";
            for (unsigned int j = 0; j < k; ++j) {
                std::cerr << shares[i][j] << " ";
            }
            std::cerr << "\n";
        }
#endif

		vector<RecoveryPoly<FX> > recs;
		vector<RecoveryPolyMulti<FX> > recs_multi;

		struct timeval st, et;
		vec_F sharesVec;
		switch (testType) {
		case CH_MULTI:
		case CH_TK1:
			gettimeofday(&st, NULL);
			recs_multi = decoder.findpolys_multi(t, h, allservers, indices, shares, testType);
			gettimeofday(&et, NULL);
			break;
        case DP:
		case BRUTE:
		case BW:
		case KOTTER:
		case CH_MS:
		case UNDEFINED:
		case UNKNOWN:
		case MAX_TESTTYPE:
		default:
			sharesVec = shares[0];
			gettimeofday(&st, NULL);
			recs = decoder.findpolys(t, h, allservers, indices, sharesVec, testType, dpType, gord);
			gettimeofday(&et, NULL);
			break;
		}
		uint64_t elapsedus = ((uint64_t)(et.tv_sec -
			st.tv_sec)) * 1000000 + (et.tv_usec - st.tv_usec);

	    int numres = 0, numcorrect = 0;
	    switch (testType) {
	    case CH_MULTI:
	    case CH_TK1:
			for(typename vector<RecoveryPolyMulti<FX> >::const_iterator iter = recs_multi.begin();
				iter != recs_multi.end(); ++iter) {
				int numcorrect_multi = 0;
				for (int i = 0; i < s; ++i) {
					if (iter->phis[i] == randpolys[i]) {
						++numcorrect_multi;
					}
				}
				if (numcorrect_multi == s) {
					std::cout << "Correct:";
					for (int i = 0; i < s; ++i) {
						std::cout << " " << iter->phis[i];
					}
					std::cout << "\n";
					++numcorrect;
				}
				++numres;
			}

#if TEST_TK1_ADD1
            while (numcorrect == 0) {
                std::cerr << "TK1 FAILED!  TRYING AGAIN WITH ONE MORE POLY ...\n";
#if defined(USE_GF28) || defined(USE_GF24)
                randpolys.push_back(random_GF2EX(t+1));
#else
                randpolys.push_back(random_ZZ_pX(t+1));
#endif
                shares.push_back(vec_F());
                shares[s].SetLength(k);

                for (int i=0; i<k; ++i) {
#if defined(USE_GF28) || defined(USE_GF24)
                    unsigned char b = i+1;
                    conv(indices[i], GF2XFromBytes(&b, 1));
#else
                    indices[i] = i+1;
#endif
                    eval(shares[s][i], randpolys[s], indices[i]);
                }
                
                for(vector<unsigned short>::iterator iter = wrongservers.begin();
                        iter != wrongservers.end(); ++iter) {
                    if (!rand_when || (rand() % 2)) {
                        F randval;
                        if (rand_amount) {
                            random(randval);
                            shares[s][*iter] += randval;
                        } else if (rand_srv) {
                            shares[s][*iter] += srvvals[*iter] * (s ? 1 : *iter);
                        } else {
                            shares[s][*iter] += (s ? 1 : *iter);
                        }
                    }
                }

                std::cerr << "shares (Second Time) =\n";
                for (unsigned int i = 0; i <= s; ++i) {
                    std::cerr << "\t";
                    for (unsigned int j = 0; j < k; ++j) {
                        std::cerr << shares[i][j] << " ";
                    }
                    std::cerr << "\n";
                }

                recs_multi = decoder.findpolys_multi(t, h, allservers, indices, shares, testType);
                for(typename vector<RecoveryPolyMulti<FX> >::const_iterator iter = recs_multi.begin();
                    iter != recs_multi.end(); ++iter) {
                    int numcorrect_multi = 0;
                    for (int i = 0; i < s; ++i) {
                        if (iter->phis[i] == randpolys[i]) {
                            ++numcorrect_multi;
                        }   
                    }   
                    if (numcorrect_multi == s) {
                        std::cout << "Correct:";
                        for (int i = 0; i < s; ++i) {
                            std::cout << " " << iter->phis[i];
                        }   
                        std::cout << "\n";
                        ++numcorrect;
                    }   
                    ++numres;
                }
                if (numcorrect == 0) {
                    std::cerr << "STILL FAILED ...\n";
                }

#else
	    	// If failed, more output
	    	if (numcorrect == 0) {
                std::ofstream tk1file;
	    		tk1file.open ("testing/tk1fail.txt", ios::app);
	    		tk1file << "SEED: " << seed_g << "\n";
	    		tk1file << "n = " << k << "\n";
	    		tk1file << "ell = " << t << "\n";
	    		tk1file << "h = " << h << "\n";
	    		tk1file << "m = " << s << "\n";
	    		tk1file << "t = k = 1" << "\n";
	    		tk1file << "Original Polys:\n";
	    		for (int i = 0; i < s; ++i) {
	    			tk1file << "\t" << i << ":\t" << randpolys[i] << "\n";
	    		}
	    		tk1file << "indices:\n";
	    		for (int i = 0; i < k; ++i) {
	    			tk1file << "\t" << i << ":\t" << indices[i] << "\n";
	    		}
	    		tk1file << "shares (poly #, index):\n";
	    		for (int i = 0; i < s; ++i) {
	    			for (int j = 0; j < k; ++j) {
	    				tk1file << "\t" << i << ",\t" << j << ":\t" << shares[i][j] << "\n";
	    			}
	    		}
	    		tk1file << "----------------------------------------------------------\n";
	    		tk1file.close();
#endif
	    	}
	    	break;
	    default:
			for(typename vector<RecoveryPoly<FX> >::const_iterator iter = recs.begin();
				iter != recs.end(); ++iter) {
				if (iter->phi == randpolys[0]) {
					std::cout << "Correct: " << iter->phi << "\n";
					++numcorrect;
				}
				++numres;
			
            }
			break;
	    }
	    cout << numres << " result" << (numres == 1 ? "" : "s") << "\n";

		ssout << elapsedus << "," << (numcorrect==0?"fail,":"pass,");
	}

    if (origTestType == BEST) {
        std::string oldssout = ssout.str();
        ssout.str("");
        size_t beststart = oldssout.find("best");
        ssout << oldssout.substr(0, beststart);
        ssout << testTypeStrings[testType];
        ssout << oldssout.substr(beststart + 4);
    }

    unsigned int kotter_m, kotter_L, ch_ms_m, ch_ms_k, bw_numcols;
    switch (testType) {
    case BRUTE:
    	break;
    case KOTTER:
    	if (!dnr) {
    		kotter_m = 1 + (unsigned int)(floor( t*k / (h*h-t*k)));
        	kotter_L = (kotter_m*h - 1)/t;
			if (getenv("PIRC_L")) kotter_L = atoi(getenv("PIRC_L"));
			if (getenv("PIRC_m")) kotter_m = atoi(getenv("PIRC_m"));
			ssout << "m=" << kotter_m << ",L=" << kotter_L << ",";
    	}
        break;
    case CH_MS:
    	if (!dnr) {
			ch_ms_m = ((h - t) * k) / (h * h - t * k) + 1;
			ch_ms_k = (h * ch_ms_m) / k;
			ssout << "m=" << ch_ms_m << ",k=" << ch_ms_k << ",";
    	}
        break;
    case BW:
    	bw_numcols = k + k - h - h + t + 2;
    	ssout << "numcols=" << bw_numcols << ",";
    	break;
    case CH_MULTI:
    	ssout << "m=" << s << ",t=" << ch_multi_t << ",k=" << ch_multi_k << ",";
    case CH_TK1:
    	ssout << "m=" << s << ",";
    	break;
    case DP:
        ssout << "dpType=" << DPTypeStrings[dpType] << ",gord=" << gord << ",";
        break;
    default:
    	break;
    }
    ssout << std::endl;

	char * testOutfile = getenv("TEST_OUTFILE");
	if (testOutfile) {
		std::ofstream timefile;
		timefile.open (testOutfile, ios::app);
		timefile << ssout.str();
		timefile.close();
	} else {
		std::cerr << "Test Info: " << ssout.str();
	}
}

void usage (int argc, char **argv) {
    std::cerr << "Usage:\n\n";
    std::cerr << "    " << argv[0] << " --help\n";
    std::cerr << "    " << argv[0] << " ( best | brute | bw | ch_ms | kotter ) k t h\n";
    std::cerr << "    " << argv[0] << " ( ch_multi | ch_tk1 ) k t h m [ multi_t multi_k ]\n";
    std::cerr << "    " << argv[0] << " dp k t h ( assume_correct | assume_wrong | assume_shares ) [ g | d ]\n\n";
    std::cerr << "    where:\n";
    std::cerr << "        k        The number of servers that respond\n";
    std::cerr << "        t        The number of servers that can collude\n";
    std::cerr << "        h        The number of honest servers\n";
    std::cerr << "        m        For the multi-polynomial cases, the number of polynomials\n";
    std::cerr << "        multi_t  For ch_multi\n";
    std::cerr << "        multi_k  For ch_multi\n";
    std::cerr << "        g        For dp with assume_correct or assume_wrong, the number of\n";
    std::cerr << "                 shares to assume correct/wrong\n";
    std::cerr << "        d        For dp with assume_shares, the number of share to make\n";
    std::cerr << "                 assumption of correct/wrong about\n\n";
}

int main(int argc, char **argv)
{
    if (strcmp(argv[1], "--help") == 0) {
        usage(argc, argv);
        return 0;
    }

	char * givenSeed = getenv("USE_SEED");
	if (givenSeed) {
		ZZ seedzz;
		std::stringstream seedstr;
		seedstr << givenSeed;
		seedstr >> seedzz;
#ifdef VERBOSE_FINDPOLYS
		std::cerr << "SEED: " << seedzz << "\n";
#endif
		SetSeed(seedzz);
		seed_g = seedzz;
	} else {
		unsigned char randbuf[128];
		ifstream urand("/dev/urandom");
		urand.read((char *)randbuf, sizeof(randbuf));
		urand.close();
		ZZ randzz = ZZFromBytes(randbuf, sizeof(randbuf));
#ifdef VERBOSE_FINDPOLYS
		std::cerr << "SEED: " << randzz << "\n";
#endif
		SetSeed(randzz);
		seed_g = randzz;
	}

	std::string testTypeStr = argc > 1 ? argv[1] : "";
	TestType testType = UNDEFINED;
	unsigned int i;
	for (i = 2; i < MAX_TESTTYPE; ++i) {
		if (testTypeStr == testTypeStrings[i]) {
			testType = (TestType)(i);
			break;
		}
	}
	if (i == MAX_TESTTYPE) {
		testType = UNKNOWN;
	}

    int k = argc > 2 ? atoi(argv[2]) : 10;
    int t = argc > 3 ? atoi(argv[3]) : 5;
    int h = argc > 4 ? atoi(argv[4]) : int(sqrt(double(k*t))) + 1;

    int s = 1, gord = 1, multi_t = 0, multi_k = 0;
    DPType dpType = UNDEFINED_DPTYPE;
    std::string dpTypeStr;
    switch (testType) {
    case CH_MULTI:
    case CH_TK1:
        s = argc > 5 ? atoi(argv[5]) : 1;
        multi_t = argc > 7 ? atoi(argv[6]) : 0;
        multi_k = argc > 7 ? atoi(argv[7]) : 0;
        break;
    case DP:
        dpTypeStr = argc > 5 ? argv[5] : "";
        for (i = 1; i < MAX_DPTYPE; ++i) {
            if (dpTypeStr == DPTypeStrings[i]) {
                dpType = (DPType)(i);
                break;
            }
        }
        gord = argc > 6 ? atoi(argv[6]) : 1;
        break;
    default:
        break;
    }

#if defined(USE_GF28) || defined(USE_GF24)
    time_findpolys<GF2E,vec_GF2E,GF2EX,GF2EXY,mat_GF2E>(testType, k, t, h, s, multi_t, multi_k, dpType, gord);
#else
    time_findpolys<ZZ_p,vec_ZZ_p,ZZ_pX,ZZ_pXY,mat_ZZ_p>(testType, k, t, h, s, multi_t, multi_k, dpType, gord);
#endif

    return 0;
}

#endif

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


