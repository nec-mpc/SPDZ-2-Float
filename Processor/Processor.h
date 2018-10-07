// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#ifndef _Processor
#define _Processor

/* This is a representation of a processing element
 *   Consisting of 256 clear and 256 shared registers
 */

#include "Math/Share.h"
#include "Math/gf2n.h"
#include "Math/gfp.h"
#include "Math/Integer.h"
#include "Exceptions/Exceptions.h"
#include "Networking/Player.h"
#include "Auth/MAC_Check.h"
#include "Data_Files.h"
#include "Input.h"
#include "PrivateOutput.h"
#include "Machine.h"
#include "ExternalClients.h"
#include "Binary_File_IO.h"
#include "Instruction.h"

#include <stack>

class ProcessorBase
{
  // Stack
  stack<long> stacki;

protected:
  // Optional argument to tape
  int arg;

public:
  void pushi(long x) { stacki.push(x); }
  void popi(long& x) { x = stacki.top(); stacki.pop(); }

  int get_arg() const
    {
      return arg;
    }

  void set_arg(int new_arg)
    {
      arg=new_arg;
    }
};

class Processor : public ProcessorBase
{
  vector<gf2n>  C2;
  vector<gfp>   Cp;
  vector<Share<gf2n> > S2;
  vector<Share<gfp> >  Sp;
  vector<long> Ci;

  // This is the vector of partially opened values and shares we need to store
  // as the Open commands are split in two
  vector<gf2n> PO2;
  vector<gfp>  POp;
  vector<Share<gf2n> > Sh_PO2;
  vector<Share<gfp> >  Sh_POp;

  int reg_max2,reg_maxp,reg_maxi;
  int thread_num;

  // Data structure used for reading/writing data to/from a socket (i.e. an external party to SPDZ)
  octetStream socket_stream;

  #ifdef DEBUG
    vector<int> rw2;
    vector<int> rwp;
    vector<int> rwi;
  #endif

  template <class T>
  vector< Share<T> >& get_S();
  template <class T>
  vector<T>& get_C();

  template <class T>
  vector< Share<T> >& get_Sh_PO();
  template <class T>
  vector<T>& get_PO();

  public:
  Data_Files& DataF;
  Player& P;
  MAC_Check<gf2n>& MC2;
  MAC_Check<gfp>& MCp;
  Machine& machine;

  string private_input_filename;

  Input<gf2n> input2;
  Input<gfp> inputp;
  
  PrivateOutput<gf2n> privateOutput2;
  PrivateOutput<gfp>  privateOutputp;

  ifstream public_input;
  ifstream private_input;
  ofstream public_output;
  ofstream private_output;

  unsigned int PC;
  TempVars temp;
  PRNG prng;

  int sent, rounds;

  ExternalClients external_clients;
  Binary_File_IO binary_file_io;
  
  static const int reg_bytes = 4;
  
  void reset(const Program& program,int arg); // Reset the state of the processor
  string get_filename(const char* basename, bool use_number);

  Processor(int thread_num,Data_Files& DataF,Player& P,
          MAC_Check<gf2n>& MC2,MAC_Check<gfp>& MCp,Machine& machine,
          const Program& program);
  ~Processor();

  int get_thread_num()
    {
      return thread_num;
    }

  #ifdef DEBUG  
    const gf2n& read_C2(int i) const
      { if (rw2[i]==0)
	  { throw Processor_Error("Invalid read on clear register"); }
        return C2.at(i);
      }
    const Share<gf2n> & read_S2(int i) const
      { if (rw2[i+reg_max2]==0)
          { throw Processor_Error("Invalid read on shared register"); }
        return S2.at(i);
      }
    gf2n& get_C2_ref(int i)
      { rw2[i]=1;
        return C2.at(i);
      }
    Share<gf2n> & get_S2_ref(int i)
      { rw2[i+reg_max2]=1;
        return S2.at(i);
      }
    void write_C2(int i,const gf2n& x)
      { rw2[i]=1;
        C2.at(i)=x;
      }
    void write_S2(int i,const Share<gf2n> & x)
      { rw2[i+reg_max2]=1;
        S2.at(i)=x;
      }

    const gfp& read_Cp(int i) const
      { if (rwp[i]==0)
	  { throw Processor_Error("Invalid read on clear register"); }
        return Cp.at(i);
      }
    const Share<gfp> & read_Sp(int i) const
      { if (rwp[i+reg_maxp]==0)
          { throw Processor_Error("Invalid read on shared register"); }
        return Sp.at(i);
      }
    gfp& get_Cp_ref(int i)
      { rwp[i]=1;
        return Cp.at(i);
      }
    Share<gfp> & get_Sp_ref(int i)
      { rwp[i+reg_maxp]=1;
        return Sp.at(i);
      }
    void write_Cp(int i,const gfp& x)
      { rwp[i]=1;
        Cp.at(i)=x;
      }
    void write_Sp(int i,const Share<gfp> & x)
      { rwp[i+reg_maxp]=1;
        Sp.at(i)=x;
      }

    const long& read_Ci(int i) const
      { if (rwi[i]==0)
          { throw Processor_Error("Invalid read on integer register"); }
        return Ci.at(i);
      }
    long& get_Ci_ref(int i)
      { rwi[i]=1;
        return Ci.at(i);
      }
    void write_Ci(int i,const long& x)
      { rwi[i]=1;
        Ci.at(i)=x;
      }
 #else
    const gf2n& read_C2(int i) const
      { return C2[i]; }
    const Share<gf2n> & read_S2(int i) const
      { return S2[i]; }
    gf2n& get_C2_ref(int i)
      { return C2[i]; }
    Share<gf2n> & get_S2_ref(int i)
      { return S2[i]; }
    void write_C2(int i,const gf2n& x)
      { C2[i]=x; }
    void write_S2(int i,const Share<gf2n> & x)
      { S2[i]=x; }
  
    const gfp& read_Cp(int i) const
      { return Cp[i]; }
    const Share<gfp> & read_Sp(int i) const
      { return Sp[i]; }
    gfp& get_Cp_ref(int i)
      { return Cp[i]; }
    Share<gfp> & get_Sp_ref(int i)
      { return Sp[i]; }
    void write_Cp(int i,const gfp& x)
      { Cp[i]=x; }
    void write_Sp(int i,const Share<gfp> & x)
      { Sp[i]=x; }

    const long& read_Ci(int i) const
      { return Ci[i]; }
    long& get_Ci_ref(int i)
      { return Ci[i]; }
    void write_Ci(int i,const long& x)
      { Ci[i]=x; }
  #endif

  // Template-based access
  template<class T> Share<T>& get_S_ref(int i);
  template<class T> T& get_C_ref(int i);

  // Access to external client sockets for reading clear/shared data
  void read_socket_ints(int client_id, const vector<int>& registers);
  // Setup client public key
  void read_client_public_key(int client_id, const vector<int>& registers);
  void init_secure_socket(int client_id, const vector<int>& registers);
  void init_secure_socket_internal(int client_id, const vector<int>& registers);
  void resp_secure_socket(int client_id, const vector<int>& registers);
  void resp_secure_socket_internal(int client_id, const vector<int>& registers);
  
  void write_socket(const RegType reg_type, const SecrecyType secrecy_type, const bool send_macs,
                             int socket_id, int message_type, const vector<int>& registers);

  template <class T>
  void read_socket_vector(int client_id, const vector<int>& registers);
  template <class T>
  void read_socket_private(int client_id, const vector<int>& registers, bool send_macs);

  // Read and write secret numeric data to file (name hardcoded at present)
  template <class T>
  void read_shares_from_file(int start_file_pos, int end_file_pos_register, const vector<int>& data_registers);
  template <class T>
  void write_shares_to_file(const vector<int>& data_registers);
  
  // Access to PO (via calls to POpen start/stop)
  template <class T>
  void POpen_Start(const vector<int>& reg,const Player& P,MAC_Check<T>& MC,int size);

  template <class T>
  void prep_shares(const vector<int>& reg, vector< Share<T> >& shares, int size);

  template <class T>
  void POpen_Stop(const vector<int>& reg,const Player& P,MAC_Check<T>& MC,int size);

  template <class T>
  void POpen_Stop_prep_opens(const vector<int>& reg, vector<T>& PO, vector<T>& C, int size);

  template <class T>
  void POpen(const vector<int>& reg,const Player& P,MAC_Check<T>& MC,int size);

  // Print the processor state
  friend ostream& operator<<(ostream& s,const Processor& P);

  private:
    void maybe_decrypt_sequence(int client_id);
    void maybe_encrypt_sequence(int client_id);

#if defined(EXTENDED_SPDZ)
  public:

  void POpen_Ext(const vector<int>& reg,int size);
  void PTriple_Ext(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c);
  void PInput_Ext(Share<gfp>& input_value, const int input_party_id);
  void PMult_Ext(const vector<int>& reg, int size);
  void PMult_Stop_prep_products(const vector<int>& reg, int size);
  void PAddm_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c);
  void PSubml_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c);
  void PSubmr_Ext(gfp& a, Share<gfp>& b, Share<gfp>& c);
  void PLdsi_Ext(gfp& value, Share<gfp>& share);
  void PBit_Ext(Share<gfp>& share);
  void PInverse_Ext(Share<gfp>& share_value, Share<gfp>& share_inverse);
  void PMulm_Ext(Share<gfp>& sec_product, const Share<gfp> & sec_factor, const gfp & clr_factor);
  void PAdds_Ext(Share<gfp>& sum, const Share<gfp>& a, const Share<gfp>& b);
  void PSubs_Ext(Share<gfp>& diff, const Share<gfp>& a, const Share<gfp>& b);

  void PShares2mpz(const vector< Share<gfp> >& shares, mpz_t * share_values);
  void Pmpz2gfps(const mpz_t * mpz_values, vector<gfp>& gfps);
  void Pmpz2share(const mpz_t * mpzv, Share<gfp> & shv);

  void * spdz_gfp_ext_handle;
  mpz_t mpz_share_aux, mpz_arg_aux;

  mpz_t * po_shares, * po_opens;
  size_t po_size;
  void alloc_po_mpz(const size_t required_size)
  {
	  po_shares = new mpz_t[po_size = required_size];
	  po_opens = new mpz_t[po_size];
	  for(size_t i = 0; i < po_size; i++) { mpz_init(po_shares[i]); mpz_init(po_opens[i]); }
  }
  void free_po_mpz()
  {
	  if(po_size > 0)
	  {
		  for(size_t i = 0; i < po_size; i++) { mpz_clear(po_shares[i]); mpz_clear(po_opens[i]); }
		  delete po_shares; po_shares = NULL;
		  delete po_opens; po_opens = NULL;
		  po_size = 0;
	  }
  }

  mpz_t * pi_inputs;
  size_t pi_size;
  void alloc_pi_mpz(const size_t required_size)
  {
	  pi_inputs = new mpz_t[pi_size = required_size];
	  for(size_t i = 0; i < pi_size; i++) mpz_init(pi_inputs[i]);
  }
  void free_pi_mpz()
  {
	  for(size_t i = 0; i < pi_size; i++) mpz_clear(pi_inputs[i]);
	  delete pi_inputs; pi_inputs = NULL;
	  pi_size = 0;
  }

  mpz_t * pm_shares, * pm_products;
  size_t pm_size;
  void alloc_pm_mpz(const size_t required_size)
  {
	  assert(required_size%2 == 0);
	  pm_size = required_size;
	  pm_shares = new mpz_t[pm_size];
	  pm_products = new mpz_t[pm_size/2];
	  for(size_t i = 0; i < pm_size/2; i++)
	  {
		  mpz_init(pm_products[i]);
		  mpz_init(pm_shares[2*i]);
		  mpz_init(pm_shares[2*i+1]);
	  }
  }
  void free_pm_mpz()
  {
	  if(pm_size > 0)
	  {
		  for(size_t i = 0; i < pm_size/2; i++)
		  {
			  mpz_clear(pm_products[i]);
			  mpz_clear(pm_shares[2*i]);
			  mpz_clear(pm_shares[2*i+1]);
		  }
		  delete pm_shares; pm_shares = NULL;
		  delete pm_products; pm_products = NULL;
		  pm_size = 0;
	  }
  }

  mpz_t * go_shares, * go_opens;
  size_t go_size;
  void alloc_go_mpz(const size_t required_size)
  {
	  go_shares = new mpz_t[go_size = required_size];
	  go_opens = new mpz_t[go_size];
	  for(size_t i = 0; i < pi_size; i++) { mpz_init(go_shares[i]); mpz_init(go_opens[i]); }
  }
  void free_go_mpz()
  {
	  for(size_t i = 0; i < po_size; i++) { mpz_clear(go_shares[i]); mpz_clear(go_opens[i]); }
	  delete go_shares; go_shares = NULL;
	  delete go_opens; go_opens = NULL;
	  go_size = 0;
  }

  mpz_t * gi_inputs;
  size_t gi_size;
  void alloc_gi_mpz(const size_t required_size)
  {
	  gi_inputs = new mpz_t[gi_size = required_size];
	  for(size_t i = 0; i < gi_size; i++) mpz_init(gi_inputs[i]);
  }
  void free_gi_mpz()
  {
	  for(size_t i = 0; i < gi_size; i++) mpz_clear(gi_inputs[i]);
	  delete gi_inputs; gi_inputs = NULL;
	  gi_size = 0;
  }

  mpz_t * gm_shares, * gm_products;
  size_t gm_size;
  void alloc_gm_mpz(const size_t required_size)
  {
	  gm_shares = new mpz_t[gm_size = required_size];
	  for(size_t i = 0; i < gm_size; i++) mpz_init(gm_shares[i]);
	  gm_products = new mpz_t[gm_size/2];
	  for(size_t i = 0; i < gm_size/2; i++) mpz_init(gm_products[i]);
  }
  void free_gm_mpz()
  {
	  for(size_t i = 0; i < gm_size; i++) mpz_clear(gm_shares[i]);
	  for(size_t i = 0; i < gm_size/2; i++) mpz_clear(gm_products[i]);
	  delete gm_shares; gm_shares = NULL;
	  delete gm_products; gm_products = NULL;
	  gm_size = 0;
  }

#endif

};

class spdz_ext_ifc
{
public:
	spdz_ext_ifc();
	~spdz_ext_ifc();

	void * x_lib_handle;

	int (*x_init)(void ** handle, const int pid, const int num_of_parties, const int thread_id,
					const char * field, const int open_count, const int mult_count, const int bits_count);
    int (*x_term)(void * handle);

    int (*x_offline)(void * handle, const int offline_size);

    int (*x_opens)(void * handle, const size_t share_count, const mpz_t * shares, mpz_t * opens, int verify);

    int (*x_triple)(void * handle, mpz_t a, mpz_t b, mpz_t c);

	int (*x_verify)(void * handle, int * error);

	int (*x_input)(void * handle, const int input_of_pid, const size_t num_of_inputs, mpz_t * inputs);

	int (*x_mult)(void * handle, const size_t share_count, const mpz_t * shares, mpz_t * products, int verify);

    int (*x_mix_add)(void * handle, mpz_t share, const mpz_t scalar);

    int (*x_mix_sub_scalar)(void * handle, mpz_t share, const mpz_t scalar);

    int (*x_mix_sub_share)(void * handle, const mpz_t scalar, mpz_t share);

    int (*x_mix_mul)(void * handle, mpz_t share, const mpz_t scalar);

    int (*x_adds)(void * handle, mpz_t share1, const mpz_t share2);

    int (*x_subs)(void * handle, mpz_t share1, const mpz_t share2);

    int (*x_share_immediates)(void * handle, const int party_id, const size_t value_count, const mpz_t * values, mpz_t * shares);

    int (*x_bit)(void * handle, mp_limb_t * share);

    int (*x_inverse)(void * handle, mpz_t share_value, mpz_t share_inverse);

    static int load_extension_method(const char * method_name, void ** proc_addr, void * libhandle);
};

template<> inline Share<gf2n>& Processor::get_S_ref(int i) { return get_S2_ref(i); }
template<> inline gf2n& Processor::get_C_ref(int i)        { return get_C2_ref(i); }
template<> inline Share<gfp>& Processor::get_S_ref(int i)  { return get_Sp_ref(i); }
template<> inline gfp& Processor::get_C_ref(int i)         { return get_Cp_ref(i); }

template<> inline vector< Share<gf2n> >& Processor::get_S()       { return S2; }
template<> inline vector< Share<gfp> >& Processor::get_S()        { return Sp; }

template<> inline vector<gf2n>& Processor::get_C()                { return C2; }
template<> inline vector<gfp>& Processor::get_C()                 { return Cp; }

template<> inline vector< Share<gf2n> >& Processor::get_Sh_PO()   { return Sh_PO2; }
template<> inline vector<gf2n>& Processor::get_PO()               { return PO2; }
template<> inline vector< Share<gfp> >& Processor::get_Sh_PO()    { return Sh_POp; }
template<> inline vector<gfp>& Processor::get_PO()                { return POp; }

#endif

