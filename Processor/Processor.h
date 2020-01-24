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

#if defined(EXTENDED_SPDZ_GFP)
  public:

  void POpen_Ext(const vector<int>& reg,int size);
  void MPOpen_Ext(const vector<int>& reg,int size);
  void PTriple_Ext(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c);
  void PInput_Ext(Share<gfp>& input_value, const int input_party_id);
  void PMult_Ext(const vector<int>& reg, int size);
  void MPMult_Ext(const vector<int>& reg, int size);
  void PMult_Stop_prep_products(const vector<int>& reg, int size, const std::vector< Share<gfp> > & products);
  void PAddm_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c);
  void MPAddm_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c);
  void PSubml_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c);
  void MPSubml_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c);
  void PSubmr_Ext(gfp& a, Share<gfp>& b, Share<gfp>& c);
  void MPSubmr_Ext(gfp& a, Share<gfp>& b, Share<gfp>& c);
  void PLdsi_Ext(gfp& value, Share<gfp>& share);
  void MPLdsi_Ext(gfp& value, Share<gfp>& share);
  void PBit_Ext(Share<gfp>& share);
  void PInverse_Ext(Share<gfp>& share_value, Share<gfp>& share_inverse);
  void PMulm_Ext(Share<gfp>& sec_product, const Share<gfp> & sec_factor, const gfp & clr_factor);
  void MPMulm_Ext(Share<gfp>& sec_product, const Share<gfp> & sec_factor, const gfp & clr_factor);
  void PAdds_Ext(Share<gfp>& sum, const Share<gfp>& a, const Share<gfp>& b);
  void MPAdds_Ext(Share<gfp>& sum, const Share<gfp>& a, const Share<gfp>& b);
  void PSubs_Ext(Share<gfp>& diff, const Share<gfp>& a, const Share<gfp>& b);
  void MPSubs_Ext(Share<gfp>& diff, const Share<gfp>& a, const Share<gfp>& b);
  void PMovs_Ext(Share<gfp>& dest, const Share<gfp>& source);

  void * spdz_gfp_ext_handle;

#endif


#if defined(EXTENDED_SPDZ_GF2N)
  public:

  void GOpen_Ext(const vector<int>& reg,int size);
  void GTriple_Ext(Share<gf2n>& a, Share<gf2n>& b, Share<gf2n>& c);
  void GInput_Ext(Share<gf2n>& input_value, const int input_party_id);
  void GMult_Ext(const vector<int>& reg, int size);
  void GMult_Stop_prep_products(const vector<int>& reg, int size, const std::vector< Share<gf2n> > & products);
  void GAddm_Ext(Share<gf2n>& a, gf2n& b, Share<gf2n>& c);
  void GSubml_Ext(Share<gf2n>& a, gf2n& b, Share<gf2n>& c);
  void GSubmr_Ext(gf2n& a, Share<gf2n>& b, Share<gf2n>& c);
  void GLdsi_Ext(gf2n& value, Share<gf2n>& share);
  void GBit_Ext(Share<gf2n>& share);
  void GInverse_Ext(Share<gf2n>& share_value, Share<gf2n>& share_inverse);
  void GMulm_Ext(Share<gf2n>& sec_product, const Share<gf2n> & sec_factor, const gf2n & clr_factor);
  void GAdds_Ext(Share<gf2n>& sum, const Share<gf2n>& a, const Share<gf2n>& b);
  void GSubs_Ext(Share<gf2n>& diff, const Share<gf2n>& a, const Share<gf2n>& b);
  void GMovs_Ext(Share<gf2n>& dest, const Share<gf2n>& source);

  void * spdz_gf2n_ext_handle;

#endif

#if defined(EXTENDED_SPDZ_Z2N)
  void Skew_Bit_Decomp_Ext(const vector<int>& dest_reg, const Share<gfp>& src_ring, int size);
  void MP_Skew_Bit_Decomp_Ext(const vector<int>& dest_reg, const Share<gfp>& src_ring, int size);
  void Skew_Bit_Recomp_Ext(const vector<int>& dest_reg, const Share<gf2n>& src_bit, int size);
  void Skew_Bit_Inject_Ext(const vector<int>& dest_reg, const Share<gf2n>& src_bit, int size);
  void Skew_Ring_Recomp_Ext(Share<gfp>& dest_ring, const vector<int>& src_reg, int size);
  void MP_Skew_Ring_Recomp_Ext(Share<gfp>& dest_ring, const vector<int>& src_reg, int size);
  void MP_Skew_Bit_Inject_Ext(const vector<int>& dest_reg, const Share<gf2n>& src_bit, int size);
#endif

};

#if defined(EXTENDED_SPDZ_GFP) || defined(EXTENDED_SPDZ_GF2N)
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

    int (*x_opens)(void * handle, const size_t share_count, const uint64_t * shares, uint64_t * opens, int verify);

    int (*x_closes)(void * handle, const int party_id, const size_t value_count, const uint64_t * values, uint64_t * shares);

    int (*x_triple)(void * handle, uint64_t * a, uint64_t * b, uint64_t * c);

	int (*x_verify)(void * handle, int * error);

	int (*x_input)(void * handle, const int input_of_pid, const size_t num_of_inputs, uint64_t * inputs);

	int (*x_mult)(void * handle, const size_t share_count, const uint64_t * xshares, const uint64_t * yshares, uint64_t * products, int verify);

    int (*x_mix_add)(void * handle, const uint64_t * share, const uint64_t * scalar, uint64_t * sum);

    int (*x_mix_sub_scalar)(void * handle, const uint64_t * share, const uint64_t * scalar, uint64_t * diff);

    int (*x_mix_sub_share)(void * handle, const uint64_t * scalar, const uint64_t * share, uint64_t * diff);

    int (*x_mix_mul)(void * handle, const uint64_t * share, const uint64_t * scalar, uint64_t * product);

    int (*x_adds)(void * handle, const uint64_t * share1, const uint64_t * share2, uint64_t * sum);

    int (*x_subs)(void * handle, const uint64_t * share1, const uint64_t * share2, uint64_t * diff);

    int (*x_bit)(void * handle, uint64_t * share);

    int (*x_inverse)(void * handle, uint64_t * share_value, uint64_t * share_inverse);
#if defined(EXTENDED_SPDZ_Z2N)
    int (*x_skew_decomp)(void * handle, const size_t bits_count, const uint64_t * ring_shares, uint64_t * bit_shares);

    int (*x_skew_recomp)(void * handle, const size_t bits_count, const uint64_t * bits_shares, uint64_t * ring_shares);

	int (*x_skew_inject)(void * handle, const uint64_t * bit_shares, uint64_t * ring_shares);

	int (*x_mp_closes)(void * handle, const int party_id, const size_t value_count, const uint64_t * values, uint64_t * shares);

	int (*x_mp_opens)(void * handle, const size_t share_count, const uint64_t * shares, uint64_t * opens, int verify);

	int (*x_mp_adds)(void * handle, const uint64_t * share1, const uint64_t * share2, uint64_t * sum);

	int (*x_mp_mix_add)(void * handle, const uint64_t * share, const uint64_t * scalar, uint64_t * sum);

	int (*x_mp_subs)(void * handle, const uint64_t * share1, const uint64_t * share2, uint64_t * diff);

	int (*x_mp_mix_sub_share)(void * handle, const uint64_t * scalar, const uint64_t * share, uint64_t * diff);

	int (*x_mp_mix_sub_scalar)(void * handle, const uint64_t * share, const uint64_t * scalar, uint64_t * diff);

	int (*x_mp_mix_mul)(void * handle, const uint64_t * share, const uint64_t * scalar, uint64_t * product);

	int (*x_mp_mult)(void * handle, const size_t share_count, const uint64_t * xshares, const uint64_t * yshares, uint64_t * products, int verify);

	int (*x_mp_skew_decomp)(void * handle, const size_t bits_count, const uint64_t * ring_shares, uint64_t * bit_shares);

	int (*x_mp_skew_recomp)(void * handle, const size_t bits_count, const uint64_t * bits_shares, uint64_t * ring_shares);

	int (*x_mp_skew_inject)(void * handle, const uint64_t * bit_shares, uint64_t * ring_shares);
#endif

    static int load_extension_method(const char * method_name, void ** proc_addr, void * libhandle);
};
#endif

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

