// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

/*
 * Setup.h
 *
 */

#ifndef MATH_SETUP_H_
#define MATH_SETUP_H_

#include <iostream>
using namespace std;

#ifndef PREP_DIR
#define PREP_DIR "Player-Data/"
#endif

/*
 * Routines to create and read setup files for the finite fields
 */

// get main directory for prep. data
string get_prep_dir(int nparties, int lg2p, int gf2ndegree);

// Read online setup file for gfp and gf2n
//void read_setup(const string& dir_prefix);
//void read_setup(int nparties, int lg2p, int gf2ndegree);


#endif /* MATH_SETUP_H_ */
