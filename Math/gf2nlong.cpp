// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

/*
 * gf2n_longlong.cpp
 *
 */

#include "gf2nlong.h"

#include "Exceptions/Exceptions.h"

#include <stdint.h>
#include <wmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#include "Tools/avx_memcpy.h"


bool is_ge(__m128i a, __m128i b)
{
  word aa[2], bb[2];
  _mm_storeu_si128((__m128i*)aa, a);
  _mm_storeu_si128((__m128i*)bb, b);
//  cout << hex << "is_ge " << aa[1] << " " << bb[1] << " " << (aa[1] > bb[1]) << " ";
//  cout << aa[0] << " " << bb[0] << " " << (aa[0] >= bb[0]) << endl;
  return aa[1] == bb[1] ? aa[0] >= bb[0] : aa[1] > bb[1];
}


ostream& operator<<(ostream& s, const int128& a)
{
  word* tmp = (word*)&a.a;
  s << hex;
  s << noshowbase;
  s.width(16);
  s.fill('0');
  s << tmp[1];
  s.width(16);
  s << tmp[0] << dec;
  return s;
}

istream& operator>>(istream& s,int128& a)
{
	a = 0;
	char c;
	for(size_t i = 0; i < 2 * sizeof(__int128_t); ++i) {
		a.a *= 16;
		s >> c;
		if(!isxdigit(c)) {
			break;
		}
		c = tolower(c);
		if('0' <= c && c <= '9') {
			a.a += (c - '0');
		} else {
			a.a += (c - 'a' + 10);
		}
	}
	return s;
}



int gf2n_long::n;
int gf2n_long::t1;
int gf2n_long::t2;
int gf2n_long::t3;
int gf2n_long::l0;
int gf2n_long::l1;
int gf2n_long::l2;
int gf2n_long::l3;
int gf2n_long::nterms;
int128 gf2n_long::mask;
int128 gf2n_long::lowermask;
int128 gf2n_long::uppermask;
bool gf2n_long::rewind = false;

#define num_2_fields 1

/* Require
 *  2*(n-1)-64+t1<64
 */
int long_fields_2[num_2_fields][4] = {
    {128,7,2,1},
    };


void gf2n_long::init_field(int nn)
{
  if (nn == 0)
    {
      nn = default_length();
      cerr << "Using GF(2^" << nn << ")" << endl;
    }

  if (nn!=128) {
      throw runtime_error("Compiled for GF(2^128) only. Change parameters or compile "
          "without USE_GF2N_LONG");
  }

  int i,j=-1;
  for (i=0; i<num_2_fields && j==-1; i++)
    { if (nn==long_fields_2[i][0]) { j=i; } }
  if (j==-1) { throw invalid_params(); }

  n=nn;
  nterms=1;
  l0=128-n;
  t1=long_fields_2[j][1];
  l1=128+t1-n;
  if (long_fields_2[j][2]!=0)
    { nterms=3;
      t2=long_fields_2[j][2];
      l2=128+t2-n;
      t3=long_fields_2[j][3];
      l3=128+t3-n;
    }
  // 2^128 has a pentanomial
  // if (nterms==1 && 2*(n-1)-128+t1>=128) { throw not_implemented(); }
  // if (nterms==3 && n!=128) { throw not_implemented(); }

  mask=_mm_set_epi64x(-1,-1);
  lowermask=_mm_set_epi64x((1LL<<(64-7))-1,-1);
  uppermask=_mm_set_epi64x(((word)-1)<<(64-7),0);
}



void gf2n_long::reduce_trinomial(int128 xh,int128 xl)
{
  // Deal with xh first
  a=xl;
  a^=(xh<<l0);
  a^=(xh<<l1);

  // Now deal with last int128
  int128 hi=a>>n;
  while (hi==0)
    { a&=mask;

      a^=hi;
      a^=(hi<<t1);
      hi=a>>n;
    }
}

void gf2n_long::reduce_pentanomial(int128 xh, int128 xl)
{
  // Deal with xh first
  a=xl;
  int128 upper, lower;
  upper=xh&uppermask;
  lower=xh&lowermask;
  // Upper part
  int128 tmp = 0;
  tmp^=(upper>>(n-t1-l0));
  tmp^=(upper>>(n-t1-l1));
  tmp^=(upper>>(n-t1-l2));
  tmp^=(upper>>(n-t1-l3));
  lower^=(tmp>>(l1));
  a^=(tmp<<(n-l1));
  // Lower part
  a^=(lower<<l0);
  a^=(lower<<l1);
  a^=(lower<<l2);
  a^=(lower<<l3);

/*
  // Now deal with last int128
  int128 hi=a>>n;
  while (hi!=0)
    { a&=mask;

      a^=hi;
      a^=(hi<<t1);
      a^=(hi<<t2);
      a^=(hi<<t3);

      hi=a>>n;
    }
*/
}


class int129
{
    int128 lower;
    bool msb;

public:
    int129() : lower(_mm_setzero_si128()), msb(false) { }
    int129(int128 lower, bool msb) : lower(lower), msb(msb) { }
    int129(int128 a) : lower(a), msb(false) { }
    int129(word a)
        { *this = a; }
    int128 get_lower() { return lower; }
    int129& operator=(const __m128i& other)
        { lower = other; msb = false; return *this; }
    int129& operator=(const word& other)
        { lower = _mm_set_epi64x(0, other); msb = false; return *this; }
    bool operator==(const int129& other)
        { return (lower == other.lower) && (msb == other.msb); }
    bool operator!=(const int129& other)
        { return !(*this == other); }
    bool operator>=(const int129& other)
        { //cout << ">= " << msb << other.msb <<  (msb > other.msb) << is_ge(lower.a, other.lower.a) << endl;
        return msb == other.msb ? is_ge(lower.a, other.lower.a) : msb > other.msb; }
    int129 operator<<(int other)
        { return int129(lower << other, _mm_cvtsi128_si32(((lower >> (128-other)) & 1).a)); }
    int129& operator>>=(int other)
        { lower >>= other; lower |= (int128(msb) << (128-other)); msb = !other; return *this; }
    int129 operator^(const int129& other)
        { return int129(lower ^ other.lower, msb ^ other.msb); }
    int129& operator^=(const int129& other)
        { lower ^= other.lower; msb ^= other.msb; return *this; }
    int129 operator&(const word& other)
        { return int129(lower & other, false); }
    friend ostream& operator<<(ostream& s, const int129& a)
        { s << a.msb << a.lower; return s; }
};

void gf2n_long::invert()
{
  if (is_one())  { return; }
  if (is_zero()) { throw division_by_zero(); }

  int129 u,v=a,B=0,D=1,mod=1;

  mod^=(int129(1)<<n);
  mod^=(int129(1)<<t1);
  if (nterms==3)
    { mod^=(int129(1)<<t2);
      mod^=(int129(1)<<t3);
    }
  u=mod; v=a;

  while (u!=0)
    { while ((u&1)==0)
        { u>>=1;
          if ((B&1)!=0) { B^=mod; }
          B>>=1;
        }
      while ((v&1)==0 && v!=0)
        { v>>=1;
          if ((D&1)!=0) { D^=mod; }
          D>>=1;
        }

      if (u>=v) { u=u^v; B=B^D; }
      else      { v=v^u; D=D^B; }
   }

  a=D.get_lower();
}


void gf2n_long::randomize(PRNG& G)
{
  a=G.get_doubleword();
  a&=mask;
}


void gf2n_long::output(ostream& s,bool human) const
{
  if (human)
    { s << *this; }
  else
    { s.write((char*) &a,sizeof(__m128i)); }
}

void gf2n_long::input(istream& s,bool human)
{
  if (s.peek() == EOF)
    { if (s.tellg() == 0)
        { cout << "IO problem. Empty file?" << endl;
          throw file_error();
        }
      //throw end_of_file();
      s.clear(); // unset EOF flag
      s.seekg(0);
      if (!rewind)
        cout << "REWINDING - ONLY FOR BENCHMARKING" << endl;
      rewind = true;
    }

  if (human)
    { s >> *this; }
  else
    { s.read((char*) &a,sizeof(__m128i)); }
}
