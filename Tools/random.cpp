// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#include "Tools/random.h"
#include <stdio.h>

#include <iostream>
using namespace std;


PRNG::PRNG() : cnt(0)
{
  #ifdef USE_AES
    useC=(Check_CPU_support_AES()==0);
  #endif

    srand((int)time(NULL));
  
}

void PRNG::ReSeed()
{
	for(int i=0; i<SEED_SIZE; ++i)
		seed[i] = (octet)rand();
	InitSeed();
}


void PRNG::SetSeed(octet* inp)
{
  memcpy(seed,inp,SEED_SIZE*sizeof(octet));
  InitSeed();
}

void PRNG::SetSeed(PRNG& G)
{
  octet tmp[SEED_SIZE];
  G.get_octets(tmp, sizeof(tmp));
  SetSeed(tmp);
}

void PRNG::InitSeed()
{
  #ifdef USE_AES
     if (useC)
        { aes_schedule(KeyScheduleC,seed); }
     else
        { aes_schedule(KeySchedule,seed); }
     memset(state,0,RAND_SIZE*sizeof(octet));
     for (int i = 0; i < PIPELINES; i++)
         state[i*AES_BLK_SIZE] = i;
  #else
     memcpy(state,seed,SEED_SIZE*sizeof(octet));
  #endif
  next(); 
  //cout << "SetSeed : "; print_state(); cout << endl;
}


void PRNG::print_state() const
{
  int i;
  for (i=0; i<SEED_SIZE; i++)
    { if (seed[i]<10){ cout << "0"; }
      cout << hex << (int) seed[i]; 
    }
  cout << "\t";
  for (i=0; i<RAND_SIZE; i++)
    { if (random[i]<10) { cout << "0"; }
      cout << hex << (int) random[i]; 
    }
  cout << "\t";
  for (i=0; i<SEED_SIZE; i++)
    { if (state[i]<10) { cout << "0"; }
      cout << hex << (int) state[i];
    }
  cout << " " << dec << cnt << " : ";
}


void PRNG::hash()
{
  #ifndef USE_AES
    // Hash seed to get a random value
    blk_SHA_CTX ctx;
    blk_SHA1_Init(&ctx);
    blk_SHA1_Update(&ctx,state,SEED_SIZE);
    blk_SHA1_Final(random,&ctx);
  #else
    if (useC)
       { aes_encrypt(random,state,KeyScheduleC); }
    else
       { ecb_aes_128_encrypt<PIPELINES>((__m128i*)random,(__m128i*)state,KeySchedule); }
  #endif
  // This is a new random value so we have not used any of it yet
  cnt=0;
}



void PRNG::next()
{
  // Increment state
  for (int i = 0; i < PIPELINES; i++)
    {
      int64_t* s = (int64_t*)&state[i*AES_BLK_SIZE];
      s[0] += PIPELINES;
      if (s[0] == 0)
          s[1]++;
    }
  hash();
}


double PRNG::get_double()
{
  // We need four bytes of randomness
  if (cnt>RAND_SIZE-4) { next(); }
  unsigned int a0=random[cnt],a1=random[cnt+1],a2=random[cnt+2],a3=random[cnt+3];
  double ans=(a0+(a1<<8)+(a2<<16)+(a3<<24));
  cnt=cnt+4;
  unsigned int den=0xFFFFFFFF;
  ans=ans/den;
  //print_state(); cout << " DBLE " <<  ans << endl;
  return ans;
}


unsigned int PRNG::get_uint()
{
  // We need four bytes of randomness
  if (cnt>RAND_SIZE-4) { next(); }
  unsigned int a0=random[cnt],a1=random[cnt+1],a2=random[cnt+2],a3=random[cnt+3];
  cnt=cnt+4;
  unsigned int ans=(a0+(a1<<8)+(a2<<16)+(a3<<24));
  // print_state(); cout << " UINT " << ans << endl;
  return ans;
}



void PRNG::get_octetStream(octetStream& ans,int len)
{
  ans.resize(len);
  for (int i=0; i<len; i++)
    { ans.data[i]=get_uchar(); }
  ans.len=len;
  ans.ptr=0;
}


void PRNG::get(int& res, int n_bits, bool positive)
{
  res = get_uint();
  res &= (1 << n_bits) - 1;
  if (positive and get_bit())
    res = -res;
}
