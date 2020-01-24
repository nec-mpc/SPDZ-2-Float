// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#include "Math/gfp.h"

#include "Exceptions/Exceptions.h"

void gfp::almost_randomize(PRNG& /*G*/)
{
//  G.get_octets((octet*)a.x,t()*sizeof(mp_limb_t));
//  a.x[t()-1]&=ZpD.mask;
}

void gfp::AND(const gfp& x,const gfp& y)
{
	for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i] & y.x[i];
}

void gfp::OR(const gfp& x,const gfp& y)
{
	for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i] | y.x[i];
}

void gfp::XOR(const gfp& x,const gfp& y)
{
	for(int i=0; i<GFP_SIZE; ++i) this->x[i] = x.x[i] ^ y.x[i];
}


void gfp::SHL(const gfp& x,int n)
{
	if (!x.is_zero()) {
		assign(x.x[0] << n);
	}
	else {
		assign_zero();
	}
}

void gfp::SHR(const gfp& x,int n)
{
	if (!x.is_zero()) {
		assign(x.x[0] >> n);
	}
	else {
		assign_zero();
	}
}


