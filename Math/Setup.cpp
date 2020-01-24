// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#include "Math/Setup.h"
#include "Math/gfp.h"
#include "Math/gf2n.h"

#include "Tools/mkpath.h"

#include <fstream>



string get_prep_dir(int nparties, int lg2p, int gf2ndegree)
{
  if (gf2ndegree == 0)
    gf2ndegree = gf2n::default_length();
  stringstream ss;
  ss << PREP_DIR << nparties << "-" << lg2p << "-" << gf2ndegree << "/";
  return ss.str();
}

//// Only read enough to initialize the fields (i.e. for OT offline or online phase only)
//void read_setup(const string& /*dir_prefix*/)
//{
//	/*
//  int lg2;
//  bigint p;
//
//  string filename = dir_prefix + "Params-Data";
//  cerr << "loading params from: " << filename << endl;
//
//  // backwards compatibility hack
//  if (dir_prefix.compare("") == 0)
//    filename = string(PREP_DIR "Params-Data");
//
//  ifstream inpf(filename.c_str());
//  if (inpf.fail()) { throw file_error(filename.c_str()); }
//  inpf >> p;
//  inpf >> lg2;
//
//  inpf.close();
//
//  gfp::init_field(p);
//  gf2n::init_field(lg2);
//  */
//}
//
//void read_setup(int /*nparties*/, int /*lg2p*/, int /*gf2ndegree*/)
//{
//  //string dir = get_prep_dir(nparties, lg2p, gf2ndegree);
//  //read_setup(dir);
//}
