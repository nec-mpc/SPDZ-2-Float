// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

#ifndef _gfp
#define _gfp

#include <iostream>
using namespace std;

#include "Math/gf2n.h"
#include "Math/field_types.h"
#include "Tools/random.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/format.hpp>

namespace mp = boost::multiprecision;

/* This is a wrapper class for the modp data type
 * It is used to be interface compatible with the gfp
 * type, which then allows us to template the Share
 * data type.
 *
 * So gfp is used ONLY for the stuff in the finite fields
 * we are going to be doing MPC over, not the modp stuff
 * for the FHE scheme
 */

#define GFP_SIZE 4

class gfp
{
	uint64_t x[GFP_SIZE];

  public:

  typedef gfp value_type;

  static int t()
    { return 0; }

  static DataFieldType field_type() { return DATA_MODP; }
  static char type_char() { return 'p'; }
  static string type_string() { return "gfp"; }

  static int size() { return GFP_SIZE*sizeof(uint64_t); }

  void assign(const gfp& g)     { }
  void assign_zero()            {  memset(x, 0, GFP_SIZE*sizeof(uint64_t));}
  void assign_one()             {  memset(x, 1, GFP_SIZE*sizeof(uint64_t));}
  void assign(word aa)          { x[0] = aa; for(int i=1; i<GFP_SIZE; ++i) x[i]=0;}
  void assign(long aa)          { x[0] = aa; for(int i=1; i<GFP_SIZE; ++i) x[i]=0;}
  void assign(int aa)           { x[0] = aa; for(int i=1; i<GFP_SIZE; ++i) x[i]=0;}
  void assign(mp::uint256_t aa) { uint64_t bit_mask = 0xFFFFFFFFFFFFFFFF;
	                               x[0] = (uint64_t) (aa & bit_mask);
                                  x[1] = (uint64_t) ((aa >> 64) & bit_mask);
                                  x[2] = (uint64_t) ((aa >> 128) & bit_mask);
                                  x[3] = (uint64_t) ((aa >> 192) & bit_mask);
  	  	  	  	  	  	  	  	  	 }
  void assign(const char*) {}

  gfp()              { }
  gfp(const gfp& g)  { memcpy(&x[0], &g.x[0], GFP_SIZE*sizeof(uint64_t)); }
  gfp(const __m128i& x) { *this=x; }
  gfp(const int128& x) { *this=x.a; }
  gfp(int x)         { assign(x); }
  ~gfp()             { ; }

  uint64_t get() const {return x[0];}

  gfp& operator=(const gfp& g)
    { if (&g!=this) {
    	for (int i=0; i<GFP_SIZE; ++i) x[i] = g.x[i];
    }
      return *this;
    }

  gfp& operator=(const __m128i /*other*/)
    {
      return *this;
    }

  void to_m128i(__m128i& /*ans*/)
    {
//      memcpy(&ans, a.x, sizeof(ans));
    }

  __m128i to_m128i()
    {
      __m128i nil;
      return nil;
    }


  bool is_zero() const            { return true;  }
  bool is_one()  const            { return false; }
  bool is_bit()  const            { return is_zero() or is_one(); }
  bool equal(const gfp& /*y*/) const  { return true; }
  bool operator==(const gfp& y) const { return equal(y); }
  bool operator!=(const gfp& y) const { return !equal(y); }

  // x+y
  template <int T>
  void add(const gfp& x,const gfp& y)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i]+y.x[i];
    }
  template <int T>
  void add(const gfp& x)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] += x.x[i];
	}
  template <int T>
  void add(void* /*x*/)
    { /*ZpD.Add<T>(a.x,a.x,(mp_limb_t*)x); */}
  template <int T>
  void add(octetStream& os)
    { add<T>(os.consume(size())); }
  void add(const gfp& x,const gfp& y)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i]+y.x[i];
	}
  void add(const gfp& x)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] += x.x[i];
    }
  void add(void* /*x*/)
    { /*ZpD.Add(a.x,a.x,(mp_limb_t*)x);*/ }
  void sub(const gfp& x,const gfp& y)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i]-y.x[i];
    }
  void sub(const gfp& x)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] -= x.x[i];
    }

  void mul(const gfp& x,const gfp& y)
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i]*y.x[i];
    }
  void mul(const gfp& x) 
    {
	  for(int i=0; i<GFP_SIZE; ++i) this->x[i] *= x.x[i];
    }

  gfp operator+(const gfp& x) { gfp res; res.add(*this, x); return res; }
  gfp operator-(const gfp& x) { gfp res; res.sub(*this, x); return res; }
  gfp operator*(const gfp& x) { gfp res; res.mul(*this, x); return res; }
  gfp& operator+=(const gfp& x) { add(x); return *this; }
  gfp& operator-=(const gfp& x) { sub(x); return *this; }
  gfp& operator*=(const gfp& x) { mul(x); return *this; }

  gfp operator-() { gfp res = *this; res.negate(); return res; }

  void square(const gfp& /*aa*/)
    { /*Sqr(a,aa.a,ZpD);*/ }
  void square()
    { /*Sqr(a,a,ZpD);*/ }
  void invert()
    { /*Inv(a,a,ZpD);*/ }
  void invert(const gfp& /*aa*/)
    { /*Inv(a,aa.a,ZpD);*/ }
  void negate() 
    { /*Negate(a,a,ZpD);*/ }
  void power(long i)
    { /*Power(a,a,i,ZpD);*/ }


  void randomize(PRNG& /*G*/)
    { /*a.randomize(G,ZpD);*/ }
  // faster randomization, see implementation for explanation
  void almost_randomize(PRNG& G);

  void output(ostream& s,bool /*human*/) const
    {
#ifndef MP_PRINT
	  s << x[0];
#else
	  mp::uint256_t res;
	  res = 0;
	  for(size_t i=0; i<4; ++i)
	  {
		  mp::uint256_t tmp = (mp::uint256_t) (x[i]);
		  res += (tmp << (64 * i));
	  }
	  s << res;
#endif
    }
  void input(istream& s,bool human)
    {
	  s >> x[0];
    }

  friend ostream& operator<<(ostream& s,const gfp& x)
    { x.output(s,true);
      return s;
    }
  friend istream& operator>>(istream& s,gfp& x)
    { x.input(s,true);
      return s;
    }

  /* Bitwise Ops 
   *   - Converts gfp args to bigints and then converts answer back to gfp
   */
  void AND(const gfp& x,const gfp& y);
  void XOR(const gfp& x,const gfp& y);
  void OR(const gfp& x,const gfp& y);
  void SHL(const gfp& x,int n);
  void SHR(const gfp& x,int n);

  gfp operator&(const gfp& x) { gfp res; res.AND(*this, x); return res; }
  gfp operator^(const gfp& x) { gfp res; res.XOR(*this, x); return res; }
  gfp operator|(const gfp& x) { gfp res; res.OR(*this, x); return res; }
  gfp operator<<(int i) { gfp res; res.SHL(*this, i); return res; }
  gfp operator>>(int i) { gfp res; res.SHR(*this, i); return res; }

  void pack(octetStream& os) const {  os.append((octet*) &x[0], GFP_SIZE*sizeof(uint64_t)); }
  void unpack(octetStream& os) { os.consume((octet*) &x[0], GFP_SIZE*sizeof(uint64_t)); }
};


#endif
