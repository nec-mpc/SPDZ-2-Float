// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#include "Processor/Processor.h"
#include "Auth/MAC_Check.h"

#include <string>

#if defined(EXTENDED_SPDZ_GFP) || defined(EXTENDED_SPDZ_GF2N)
#include <sys/stat.h>
#include <dlfcn.h>
#include <list>

spdz_ext_ifc the_ext_lib;
#endif

Processor::Processor(int thread_num,Data_Files& DataF,Player& P,
        MAC_Check<gf2n>& MC2,MAC_Check<gfp>& MCp,Machine& machine,
        const Program& program)
: thread_num(thread_num),DataF(DataF),P(P),MC2(MC2),MCp(MCp),machine(machine),
  private_input_filename(get_filename(PREP_DIR "Private-Input-",true)),
  input2(*this,MC2),inputp(*this,MCp),privateOutput2(*this),privateOutputp(*this),sent(0),rounds(0),
  external_clients(ExternalClients(P.my_num(), DataF.prep_data_dir)),binary_file_io(Binary_File_IO())
{
  reset(program,0);

  public_input.open(get_filename("Programs/Public-Input/",false).c_str());
  private_input.open(private_input_filename.c_str());
  public_output.open(get_filename(PREP_DIR "Public-Output-",true).c_str(), ios_base::out);
  private_output.open(get_filename(PREP_DIR "Private-Output-",true).c_str(), ios_base::out);

#if defined(EXTENDED_SPDZ_GFP)
    spdz_gfp_ext_handle = NULL;
	cout << "Processor " << thread_num << " SPDZ GFP extension library initializing." << endl;
//	if(0 != (*the_ext_lib.x_init)(&spdz_gfp_ext_handle, P.my_num(), P.num_players(), thread_num, "gfp127", 100, 100, 100))
	if(0 != (*the_ext_lib.x_init)(&spdz_gfp_ext_handle, P.my_num(), P.num_players(), thread_num, "Z2n_Ring", 100, 100, 100))
	{
		cerr << "SPDZ GFP extension library initialization failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
	cout << "SPDZ GFP extension library initialized." << endl;
#endif

#if defined(EXTENDED_SPDZ_GF2N)
	spdz_gf2n_ext_handle = NULL;
	cout << "Processor " << thread_num << " SPDZ GF2N extension library initializing." << endl;
//	if(0 != (*the_ext_lib.x_init)(&spdz_gf2n_ext_handle, P.my_num(), P.num_players(), thread_num, "gf2n64", 0, 0, 100))
	if(0 != (*the_ext_lib.x_init)(&spdz_gf2n_ext_handle, P.my_num(), P.num_players(), thread_num, "Z2_Bool", 0, 0, 100))
	{
		cerr << "SPDZ GF2N extension library initialization failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
	cout << "SPDZ GF2N extension library initialized." << endl;
#endif
}

Processor::~Processor()
{
  cerr << "Sent " << sent << " elements in " << rounds << " rounds" << endl;

#if defined(EXTENDED_SPDZ_GFP)
	(*the_ext_lib.x_term)(spdz_gfp_ext_handle);
#endif

#if defined(EXTENDED_SPDZ_GF2N)
	(*the_ext_lib.x_term)(spdz_gf2n_ext_handle);
#endif

#if defined(EXTENDED_SPDZ_GFP) || defined(EXTENDED_SPDZ_GF2N)
	dlclose(the_ext_lib.x_lib_handle);
#endif
}

string Processor::get_filename(const char* prefix, bool use_number)
{
  stringstream filename;
  filename << prefix;
  if (!use_number)
    filename << machine.progname;
  if (use_number)
    filename << P.my_num();
  if (thread_num > 0)
    filename << "-" << thread_num;
  cerr << "Opening file " << filename.str() << endl;
  return filename.str();
}


void Processor::reset(const Program& program,int arg)
{
  reg_max2 = program.num_reg(GF2N);
  reg_maxp = program.num_reg(MODP);
  reg_maxi = program.num_reg(INT);
  C2.resize(reg_max2); Cp.resize(reg_maxp);
  S2.resize(reg_max2); Sp.resize(reg_maxp);
  Ci.resize(reg_maxi);
  this->arg = arg;

  #ifdef DEBUG
    rw2.resize(2*reg_max2);
    for (int i=0; i<2*reg_max2; i++) { rw2[i]=0; }
    rwp.resize(2*reg_maxp);
    for (int i=0; i<2*reg_maxp; i++) { rwp[i]=0; }
    rwi.resize(2*reg_maxi);
    for (int i=0; i<2*reg_maxi; i++) { rwi[i]=0; }
  #endif
}

#include "Networking/sockets.h"
#include "Math/Setup.h"

// Write socket (typically SPDZ engine -> external client), for different register types.
// RegType and SecrecyType determines how registers are read and the socket stream is packed.
// If message_type is > 0, send message_type in bytes 0 - 3, to allow an external client to
//  determine the data structure being sent in a message.
void Processor::write_socket(const RegType reg_type, const SecrecyType secrecy_type, const bool send_macs,
                             int socket_id, int message_type, const vector<int>& registers)
{
  if (socket_id >= (int)external_clients.external_client_sockets.size())
  {
    cerr << "No socket connection exists for client id " << socket_id << endl;
    return;  
  }
  int m = registers.size();
  socket_stream.reset_write_head();

  //First 4 bytes is message_type (unless indicate not needed)
  if (message_type != 0) {
    socket_stream.store(message_type);
  }

  for (int i = 0; i < m; i++)
  {
    if (reg_type == MODP && secrecy_type == SECRET) {
      // Send vector of secret shares and optionally macs
      get_S_ref<gfp>(registers[i]).get_share().pack(socket_stream);
      if (send_macs)
        get_S_ref<gfp>(registers[i]).get_mac().pack(socket_stream);
    }
    else if (reg_type == MODP && secrecy_type == CLEAR) {
      // Send vector of clear public field elements
      get_C_ref<gfp>(registers[i]).pack(socket_stream);
    }
    else if (reg_type == INT && secrecy_type == CLEAR) {
      // Send vector of 32-bit clear ints
      socket_stream.store((int&)get_Ci_ref(registers[i]));
    } 
    else {
      stringstream ss;
      ss << "Write socket instruction with unknown reg type " << reg_type << 
        " and secrecy type " << secrecy_type << "." << endl;      
      throw Processor_Error(ss.str());
    }
  }

  try {
    socket_stream.Send(external_clients.external_client_sockets[socket_id]);
  }
    catch (bad_value& e) {
    cerr << "Send error thrown when writing " << m << " values of type " << reg_type << " to socket id " 
      << socket_id << "." << endl;
  }
}


// Receive vector of 32-bit clear ints
void Processor::read_socket_ints(int client_id, const vector<int>& registers)
{
  if (client_id >= (int)external_clients.external_client_sockets.size())
  {
    cerr << "No socket connection exists for client id " << client_id << endl; 
    return; 
  }

  int m = registers.size();
  socket_stream.reset_write_head();
  socket_stream.Receive(external_clients.external_client_sockets[client_id]);
  for (int i = 0; i < m; i++)
  {
    int val;
    socket_stream.get(val);
    write_Ci(registers[i], (long)val);
  }
}

// Receive vector of public field elements
template <class T>
void Processor::read_socket_vector(int client_id, const vector<int>& registers)
{
  if (client_id >= (int)external_clients.external_client_sockets.size())
  {
    cerr << "No socket connection exists for client id " << client_id << endl;
    return;  
  }

  int m = registers.size();
  socket_stream.reset_write_head();
  socket_stream.Receive(external_clients.external_client_sockets[client_id]);
  for (int i = 0; i < m; i++)
  {
    get_C_ref<T>(registers[i]).unpack(socket_stream);
  }
}

// Receive vector of field element shares over private channel
template <class T>
void Processor::read_socket_private(int client_id, const vector<int>& registers, bool read_macs)
{
  if (client_id >= (int)external_clients.external_client_sockets.size())
  {
    cerr << "No socket connection exists for client id " << client_id << endl;
    return;  
  }
  int m = registers.size();
  socket_stream.reset_write_head();
  socket_stream.Receive(external_clients.external_client_sockets[client_id]);

  for (int i = 0; i < m; i++)
  {
    temp.ansp.unpack(socket_stream);
    get_Sp_ref(registers[i]).set_share(temp.ansp);
    if (read_macs)
    {
      temp.ansp.unpack(socket_stream);
      get_Sp_ref(registers[i]).set_mac(temp.ansp);
    }
  }
}

// Read socket for client public key as 8 ints, calculate session key for client.
void Processor::read_client_public_key(int client_id, const vector<int>& registers) {

  read_socket_ints(client_id, registers);

  // After read into registers, need to extract values
  vector<int> client_public_key (registers.size(), 0);
  for(unsigned int i = 0; i < registers.size(); i++) {
    client_public_key[i] = (int&)get_Ci_ref(registers[i]);
  }
}


// Read share data from a file starting at file_pos until registers filled.
// file_pos_register is written with new file position (-1 is eof).
// Tolerent to no file if no shares yet persisted.
template <class T> 
void Processor::read_shares_from_file(int start_file_posn, int end_file_pos_register, const vector<int>& data_registers) {
  string filename;
  filename = "Persistence/Transactions-P" + to_string(P.my_num()) + ".data";

  unsigned int size = data_registers.size();

  vector< Share<T> > outbuf(size);

  int end_file_posn = start_file_posn;

  try {
    binary_file_io.read_from_file<T>(filename, outbuf, start_file_posn, end_file_posn);

    for (unsigned int i = 0; i < size; i++)
    {
      get_Sp_ref(data_registers[i]).set_share(outbuf[i].get_share());
      get_Sp_ref(data_registers[i]).set_mac(outbuf[i].get_mac());
    }

    write_Ci(end_file_pos_register, (long)end_file_posn);    
  }
  catch (file_missing& e) {
    cerr << "Got file missing error, will return -2. " << e.what() << endl;
    write_Ci(end_file_pos_register, (long)-2);
  }
}

// Append share data in data_registers to end of file. Expects Persistence directory to exist.
template <class T>
void Processor::write_shares_to_file(const vector<int>& data_registers) {
  string filename;
  filename = "Persistence/Transactions-P" + to_string(P.my_num()) + ".data";

  unsigned int size = data_registers.size();

  vector< Share<T> > inpbuf (size);

  for (unsigned int i = 0; i < size; i++)
  {
    inpbuf[i] = get_S_ref<T>(data_registers[i]);
  }

  binary_file_io.write_to_file<T>(filename, inpbuf);
}

template <class T>
void Processor::prep_shares(const vector<int>& reg, vector< Share<T> >& shares, int size)
{
	if (size>1)
	{
		for (typename vector<int>::const_iterator reg_it=reg.begin(); reg_it!=reg.end(); reg_it++)
		{
			typename vector<Share<T> >::iterator begin=get_S<T>().begin()+*reg_it;
			shares.insert(shares.end(),begin,begin+size);
		}
	}
	else
	{
		int sz=reg.size();
		for (int i=0; i<sz; i++)
		{
			shares.push_back(get_S_ref<T>(reg[i]));
		}
	}
}

template <class T>
void Processor::POpen_Stop_prep_opens(const vector<int>& reg, vector<T>& PO, vector<T>& C, int size)
{
	if (size>1)
	{
		typename vector<T>::iterator PO_it=PO.begin();
		for (typename vector<int>::const_iterator reg_it=reg.begin(); reg_it!=reg.end(); reg_it++)
		{
			for (typename vector<T>::iterator C_it=C.begin()+*reg_it; C_it!=C.begin()+*reg_it+size; C_it++)
			{
			  *C_it=*PO_it;
			  PO_it++;
			}
		}
	}
	else
	{
		for (unsigned int i=0; i<reg.size(); i++)
		{
			get_C_ref<T>(reg[i]) = PO[i];
		}
	}
}

template <class T>
void Processor::POpen_Start(const vector<int>& reg,const Player& P,MAC_Check<T>& MC,int size)
{
	int sz=reg.size();

	vector< Share<T> >& Sh_PO = get_Sh_PO<T>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(reg, Sh_PO, size);

	vector<T>& PO = get_PO<T>();
	PO.resize(sz*size);

	MC.POpen_Begin(PO,Sh_PO,P);
}


template <class T>
void Processor::POpen_Stop(const vector<int>& reg,const Player& P,MAC_Check<T>& MC,int size)
{
	vector< Share<T> >& Sh_PO = get_Sh_PO<T>();
	vector<T>& PO = get_PO<T>();
	vector<T>& C = get_C<T>();
	int sz=reg.size();
	PO.resize(sz*size);
	MC.POpen_End(PO,Sh_PO,P);

	POpen_Stop_prep_opens(reg, PO, C, size);

	sent += reg.size() * size;
	rounds++;
}

void unzip_open(vector<int>& dest, vector<int>& source, const vector<int>& reg)
{
	int n = reg.size() / 2;
	source.resize(n);
	dest.resize(n);
	for (int i = 0; i < n; i++)
	{
		source[i] = reg[2 * i + 1];
		dest[i] = reg[2 * i];
	}
}

template<class T>
void Processor::POpen(const vector<int>& reg, const Player& P, MAC_Check<T>& MC,
		int size)
{
	vector<int> source, dest;
	unzip_open(dest, source, reg);
	POpen_Start(source, P, MC, size);
	POpen_Stop(dest, P, MC, size);
}

ostream& operator<<(ostream& s,const Processor& P)
{
  s << "Processor State" << endl;
  s << "Char 2 Registers" << endl;
  s << "Val\tClearReg\tSharedReg" << endl;
  for (int i=0; i<P.reg_max2; i++)
    { s << i << "\t";
      P.read_C2(i).output(s,true);
      s << "\t";
      P.read_S2(i).output(s,true);
      s << endl;
    }
  s << "Char p Registers" << endl;
  s << "Val\tClearReg\tSharedReg" << endl;
  for (int i=0; i<P.reg_maxp; i++)
    { s << i << "\t";
      P.read_Cp(i).output(s,true);
      s << "\t";
      P.read_Sp(i).output(s,true);
      s << endl;
    }

  return s;
}


#if defined(EXTENDED_SPDZ_GFP)

void Processor::POpen_Ext(const vector<int>& reg, int size)
{
	vector<int> dest, source;
	unzip_open(dest, source, reg);

	int sz=source.size();
	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sz*size);
	vector<gfp>& C = get_C<gfp>();
	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(source, Sh_PO, size);

	//the extension library is given the shares' values and returns opens' values
	if(0 != (*the_ext_lib.x_opens)(spdz_gfp_ext_handle, Sh_PO.size(), (const uint64_t*)Sh_PO.data(), (uint64_t*)PO.data(), 1))
	{
		cerr << "Processor::POpen_Ext extension library open failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	POpen_Stop_prep_opens(dest, PO, C, size);

	sent += dest.size() * size;
	rounds++;
}

void Processor::MPOpen_Ext(const vector<int>& reg, int size)
{
	vector<int> dest, source;
	unzip_open(dest, source, reg);

	int sz=source.size();
	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sz*size);
	vector<gfp>& C = get_C<gfp>();
	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(source, Sh_PO, size);

	//the extension library is given the shares' values and returns opens' values
	if(0 != (*the_ext_lib.x_mp_opens)(spdz_gfp_ext_handle, Sh_PO.size(), (const uint64_t*)Sh_PO.data(), (uint64_t*)PO.data(), 1))
	{
		cerr << "Processor::MPOpen_Ext extension library open failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	POpen_Stop_prep_opens(dest, PO, C, size);

	sent += dest.size() * size;
	rounds++;
}


void Processor::PTriple_Ext(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_triple)(spdz_gfp_ext_handle, (uint64_t*)&a, (uint64_t*)&b, (uint64_t*)&c))
	{
		cerr << "Processor::PTriple_Ext extension library triple failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PInput_Ext(Share<gfp>& input_value, const int input_party_id)
{
	if(0 != (*the_ext_lib.x_input)(spdz_gfp_ext_handle, input_party_id, 1, (uint64_t*)&input_value))
	{
		cerr << "Processor::PInput_Ext extension library input failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PMult_Ext(const vector<int>& reg, int size)
{
	vector<int> xsources, ysources, dest;
	int n = reg.size() / 3;
	xsources.reserve(n);
	ysources.reserve(n);
	dest.reserve(n);

	for (int i = 0; i < n; i++)
	{
		dest.push_back(reg[3 * i]);
		xsources.push_back(reg[3 * i + 1]);
		ysources.push_back(reg[3 * i + 2]);
	}

	std::vector< Share<gfp> > xsh, ysh, products(n*size);
	xsh.reserve(n*size);
	ysh.reserve(n*size);

	prep_shares(xsources, xsh, size);
	prep_shares(ysources, ysh, size);

	if(0 != (*the_ext_lib.x_mult)(spdz_gfp_ext_handle, n*size, (const uint64_t*)xsh.data(), (const uint64_t*)ysh.data(), (uint64_t*)products.data(), 1))
	{
		cerr << "Processor::PMult_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	PMult_Stop_prep_products(dest, size, products);

	sent += dest.size() * size;
	rounds++;
}

void Processor::MPMult_Ext(const vector<int>& reg, int size)
{
	vector<int> xsources, ysources, dest;
	int n = reg.size() / 3;
	xsources.reserve(n);
	ysources.reserve(n);
	dest.reserve(n);

	for (int i = 0; i < n; i++)
	{
		dest.push_back(reg[3 * i]);
		xsources.push_back(reg[3 * i + 1]);
		ysources.push_back(reg[3 * i + 2]);
	}

	// size = 1
	std::vector< Share<gfp> > xsh, ysh, products(n*size);
	xsh.reserve(n*size);
	ysh.reserve(n*size);

	prep_shares(xsources, xsh, size);
	prep_shares(ysources, ysh, size);
	if(0 != (*the_ext_lib.x_mp_mult)(spdz_gfp_ext_handle, n*size, (const uint64_t*)xsh.data(), (const uint64_t*)ysh.data(), (uint64_t*)products.data(), 1))
	{
		cerr << "Processor::MPMult_Ext extension library mp_mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
	PMult_Stop_prep_products(dest, size, products);
	sent += dest.size() * size;
	rounds++;
}


void Processor::PMult_Stop_prep_products(const vector<int>& reg, int size, const std::vector< Share<gfp> > & products)
{
	std::vector< Share<gfp> >::const_iterator sitr = products.begin();
	if (size>1)
	{
		for (typename vector<int>::const_iterator reg_it=reg.begin(); reg_it!=reg.end(); reg_it++)
		{
			vector< Share<gfp> >::iterator insert_point=get_S<gfp>().begin()+*reg_it;
			for(int i = 0; i < size; ++i)
			{
				*(insert_point + i) = *sitr++;
			}
		}
	}
	else
	{
		int sz=reg.size();
		for(int i = 0; i < sz; ++i)
		{
			get_S_ref<gfp>(reg[i]) = *sitr++;
//			get_S_ref<gfp>(reg[i]) = *(sitr + i);
		}
	}
}

void Processor::PAddm_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_mix_add)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::PAddm_Ext extension library mix_add failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPAddm_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_mp_mix_add)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::MPAddm_Ext extension library mix_add failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PSubml_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_mix_sub_scalar)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::PSubml_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPSubml_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_mp_mix_sub_scalar)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::MPSubml_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PSubmr_Ext(gfp& a, Share<gfp>& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_mix_sub_share)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::PSubmr_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPSubmr_Ext(gfp& a, Share<gfp>& b, Share<gfp>& c)
{
	if(0 != (*the_ext_lib.x_mp_mix_sub_share)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::MPSubmr_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PLdsi_Ext(gfp& value, Share<gfp>& share)
{
	if(0 != (*the_ext_lib.x_closes)(spdz_gfp_ext_handle, 0, 1, (const uint64_t *)&value, (uint64_t *)&share))
	{
		cerr << "Processor::PLdsi_Ext extension library closes failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPLdsi_Ext(gfp& value, Share<gfp>& share)
{
	if(0 != (*the_ext_lib.x_mp_closes)(spdz_gfp_ext_handle, 0, 1, (const uint64_t *)&value, (uint64_t *)&share))
	{
		cerr << "Processor::PLdsi_Ext extension library closes failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PBit_Ext(Share<gfp>& share)
{
	if(0 != (*the_ext_lib.x_bit)(spdz_gfp_ext_handle, (uint64_t *)&share))
	{
		cerr << "Processor::PBit_Ext extension library bit failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PInverse_Ext(Share<gfp>& share_value, Share<gfp>& share_inverse)
{
	if(0 != (*the_ext_lib.x_inverse)(spdz_gfp_ext_handle, (uint64_t *)&share_value, (uint64_t *)&share_inverse))
	{
		cerr << "Processor::PInverse_Ext extension library inverse failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PMulm_Ext(Share<gfp>& product, const Share<gfp>& share, const gfp & scalar)
{

	if(0 != (*the_ext_lib.x_mix_mul)(spdz_gfp_ext_handle, (const uint64_t *)&share, (const uint64_t *)&scalar, (uint64_t *)&product))
	{
		cerr << "Processor::PMulm_Ext extension library mix_mul failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPMulm_Ext(Share<gfp>& product, const Share<gfp>& share, const gfp & scalar)
{

	if(0 != (*the_ext_lib.x_mp_mix_mul)(spdz_gfp_ext_handle, (const uint64_t *)&share, (const uint64_t *)&scalar, (uint64_t *)&product))
	{
		cerr << "Processor::MPMulm_Ext extension library mp_mix_mul failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PAdds_Ext(Share<gfp>& sum, const Share<gfp>& a, const Share<gfp>& b)
{
	if(0 != (*the_ext_lib.x_adds)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&sum))
	{
		cerr << "Processor::PAdds_Ext extension library x_adds failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPAdds_Ext(Share<gfp>& sum, const Share<gfp>& a, const Share<gfp>& b)
{
	if(0 != (*the_ext_lib.x_mp_adds)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&sum))
	{
		cerr << "Processor::MPAdds_Ext extension library x_mp_adds failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PSubs_Ext(Share<gfp>& diff, const Share<gfp>& a, const Share<gfp>& b)
{
	if(0 != (*the_ext_lib.x_subs)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&diff))
	{
		cerr << "Processor::PSubs_Ext extension library x_subs failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MPSubs_Ext(Share<gfp>& diff, const Share<gfp>& a, const Share<gfp>& b)
{
	if(0 != (*the_ext_lib.x_mp_subs)(spdz_gfp_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&diff))
	{
		cerr << "Processor::MPSubs_Ext extension library x_mp_subs failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PMovs_Ext(Share<gfp>& dest, const Share<gfp>& source)
{
	memcpy(&dest, &source, 4 * sizeof(uint64_t));
}

#endif

#if defined(EXTENDED_SPDZ_GF2N)

void Processor::GOpen_Ext(const vector<int>& reg,int size)
{
	vector<int> dest, source;
	unzip_open(dest, source, reg);

	int sz=source.size();
	vector<gf2n>& PO = get_PO<gf2n>();
	PO.resize(sz*size);
	vector<gf2n>& C = get_C<gf2n>();
	vector< Share<gf2n> >& Sh_PO = get_Sh_PO<gf2n>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(source, Sh_PO, size);

	//the extension library is given the shares' values and returns opens' values
	if(0 != (*the_ext_lib.x_opens)(spdz_gf2n_ext_handle, Sh_PO.size(), (const uint64_t*)Sh_PO.data(), (uint64_t*)PO.data(), 1))
	{
		cerr << "Processor::GOpen_Ext extension library open failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	POpen_Stop_prep_opens(dest, PO, C, size);

	sent += dest.size() * size;
	rounds++;
}

void Processor::GTriple_Ext(Share<gf2n>& a, Share<gf2n>& b, Share<gf2n>& c)
{
	if(0 != (*the_ext_lib.x_triple)(spdz_gf2n_ext_handle, (uint64_t*)&a, (uint64_t*)&b, (uint64_t*)&c))
	{
		cerr << "Processor::GTriple_Ext extension library triple failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GInput_Ext(Share<gf2n>& input_value, const int input_party_id)
{
	if(0 != (*the_ext_lib.x_input)(spdz_gf2n_ext_handle, input_party_id, 1, (uint64_t*)&input_value))
	{
		cerr << "Processor::GInput_Ext extension library input failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GMult_Ext(const vector<int>& reg, int size)
{
	vector<int> xsources, ysources, dest;
	int n = reg.size() / 3;
	xsources.reserve(n);
	ysources.reserve(n);
	dest.reserve(n);

	for (int i = 0; i < n; i++)
	{
		dest.push_back(reg[3 * i]);
		xsources.push_back(reg[3 * i + 1]);
		ysources.push_back(reg[3 * i + 2]);
	}

	std::vector< Share<gf2n> > xsh, ysh, products(n*size);
	xsh.reserve(n*size);
	ysh.reserve(n*size);

	prep_shares(xsources, xsh, size);
	prep_shares(ysources, ysh, size);

	if(0 != (*the_ext_lib.x_mult)(spdz_gf2n_ext_handle, n*size, (const uint64_t*)xsh.data(), (const uint64_t*)ysh.data(), (uint64_t*)products.data(), 1))
	{
		cerr << "Processor::GMult_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	GMult_Stop_prep_products(dest, size, products);

	sent += dest.size() * size;
	rounds++;

}

void Processor::GMult_Stop_prep_products(const vector<int>& reg, int size, const std::vector< Share<gf2n> > & products)
{
	std::vector< Share<gf2n> >::const_iterator sitr = products.begin();
	if (size>1)
	{
		for (typename vector<int>::const_iterator reg_it=reg.begin(); reg_it!=reg.end(); reg_it++)
		{
			vector< Share<gf2n> >::iterator insert_point=get_S<gf2n>().begin()+*reg_it;
			for(int i = 0; i < size; ++i)
			{
				*(insert_point + i) = *sitr++;
			}
		}
	}
	else
	{
		int sz=reg.size();
		for(int i = 0; i < sz; ++i)
		{
			get_S_ref<gf2n>(reg[i]) = *sitr++;
		}
	}
}

void Processor::GAddm_Ext(Share<gf2n>& a, gf2n& b, Share<gf2n>& c)
{
	if(0 != (*the_ext_lib.x_mix_add)(spdz_gf2n_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::GAddm_Ext extension library mix_add failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GSubml_Ext(Share<gf2n>& a, gf2n& b, Share<gf2n>& c)
{
	if(0 != (*the_ext_lib.x_mix_sub_scalar)(spdz_gf2n_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::GSubml_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GSubmr_Ext(gf2n& a, Share<gf2n>& b, Share<gf2n>& c)
{
	if(0 != (*the_ext_lib.x_mix_sub_share)(spdz_gf2n_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&c))
	{
		cerr << "Processor::GSubmr_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GLdsi_Ext(gf2n& value, Share<gf2n>& share)
{
	if(0 != (*the_ext_lib.x_closes)(spdz_gf2n_ext_handle, 0, 1, (const uint64_t *)&value, (uint64_t *)&share))
	{
		cerr << "Processor::GLdsi_Ext extension library closes failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GBit_Ext(Share<gf2n>& share)
{
	if(0 != (*the_ext_lib.x_bit)(spdz_gf2n_ext_handle, (uint64_t *)&share))
	{
		cerr << "Processor::GBit_Ext extension library bit failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GInverse_Ext(Share<gf2n>& share_value, Share<gf2n>& share_inverse)
{
	if(0 != (*the_ext_lib.x_inverse)(spdz_gf2n_ext_handle, (uint64_t *)&share_value, (uint64_t *)&share_inverse))
	{
		cerr << "Processor::GInverse_Ext extension library inverse failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GMulm_Ext(Share<gf2n>& sec_product, const Share<gf2n> & sec_factor, const gf2n & clr_factor)
{
	if(0 != (*the_ext_lib.x_mix_mul)(spdz_gf2n_ext_handle, (const uint64_t *)&sec_factor, (const uint64_t *)&clr_factor, (uint64_t *)&sec_product))
	{
		cerr << "Processor::PMulm_Ext extension library mix_mul failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GAdds_Ext(Share<gf2n>& sum, const Share<gf2n>& a, const Share<gf2n>& b)
{
	if(0 != (*the_ext_lib.x_adds)(spdz_gf2n_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&sum))
	{
		cerr << "Processor::GAdds_Ext extension library x_adds failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GSubs_Ext(Share<gf2n>& diff, const Share<gf2n>& a, const Share<gf2n>& b)
{
	if(0 != (*the_ext_lib.x_subs)(spdz_gf2n_ext_handle, (const uint64_t *)&a, (const uint64_t *)&b, (uint64_t *)&diff))
	{
		cerr << "Processor::GSubs_Ext extension library x_subs failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::GMovs_Ext(Share<gf2n>& dest, const Share<gf2n>& source)
{
#ifdef USE_GF2N_LONG
	memcpy(&dest, &source, 2 * sizeof(uint64_t));
#else
	memcpy(&dest, &source, sizeof(uint64_t));
#endif
}

#endif

#if defined(EXTENDED_SPDZ_Z2N)
void Processor::Skew_Bit_Decomp_Ext(const vector<int>& dest_reg, const Share<gfp>& src_ring, int size)
{
	size_t n = dest_reg.size() / 3;

	std::vector< Share<gf2n> > dest_bits(dest_reg.size());

	if(0 != (*the_ext_lib.x_skew_decomp)(spdz_gfp_ext_handle, n, (const uint64_t *)&src_ring, (uint64_t *)dest_bits.data()))
	{
		cerr << "Processor::Skew_Bit_Decomp_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	GMult_Stop_prep_products(dest_reg, size, dest_bits);
}

void Processor::MP_Skew_Bit_Decomp_Ext(const vector<int>& dest_reg, const Share<gfp>& src_ring, int size)
{
	size_t n = dest_reg.size() / 3;

	std::vector< Share<gf2n> > dest_bits(dest_reg.size());

	if(0 != (*the_ext_lib.x_mp_skew_decomp)(spdz_gfp_ext_handle, n, (const uint64_t *)&src_ring, (uint64_t *)dest_bits.data()))
	{
		cerr << "Processor::MP_Skew_Bit_Decomp_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	GMult_Stop_prep_products(dest_reg, size, dest_bits);
}


void Processor::Skew_Bit_Recomp_Ext(const vector<int>& dest_reg, const Share<gf2n>& src_bit, int size)
{
	size_t n = dest_reg.size() / 3;

	std::vector< Share<gf2n> > dest_bits(dest_reg.size());

	if(0 != (*the_ext_lib.x_skew_decomp)(spdz_gf2n_ext_handle, n, (const uint64_t*)&src_bit, (uint64_t *)dest_bits.data()))
	{
		cerr << "Processor::Skew_Bit_Recomp_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	GMult_Stop_prep_products(dest_reg, size, dest_bits);
}

void Processor::Skew_Bit_Inject_Ext(const vector<int>& dest, const Share<gf2n>& src_bit, int size)
{
	std::vector< Share<gfp> > dest_ring(dest.size());

	if(0 != (*the_ext_lib.x_skew_inject)(spdz_gf2n_ext_handle, (const uint64_t*)&src_bit, (uint64_t *)dest_ring.data()))
	{
		cerr << "Processor::Skew_Bit_Inject_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	PMult_Stop_prep_products(dest, size, dest_ring);
}

void Processor::MP_Skew_Bit_Inject_Ext(const vector<int>& dest, const Share<gf2n>& src_bit, int size)
{
	std::vector< Share<gfp> > dest_ring(dest.size());

	if(0 != (*the_ext_lib.x_mp_skew_inject)(spdz_gf2n_ext_handle, (const uint64_t*)&src_bit, (uint64_t *)dest_ring.data()))
	{
		cerr << "Processor::MP_Skew_Bit_Inject_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	PMult_Stop_prep_products(dest, size, dest_ring);
}

void Processor::Skew_Ring_Recomp_Ext(Share<gfp>& dest_ring, const vector<int>& src_reg, int size)
{
	size_t n = src_reg.size();

	std::vector< Share<gf2n> > src_bits;
	src_bits.reserve(src_reg.size());

	prep_shares(src_reg, src_bits, size);

	if(0 != (*the_ext_lib.x_skew_recomp)(spdz_gf2n_ext_handle, n, (const uint64_t*)src_bits.data(), (uint64_t *)&dest_ring))
	{
		cerr << "Processor::Skew_Bit_Recomp_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::MP_Skew_Ring_Recomp_Ext(Share<gfp>& dest_ring, const vector<int>& src_reg, int size)
{
	size_t n = src_reg.size();

	std::vector< Share<gf2n> > src_bits;
	src_bits.reserve(src_reg.size());

	prep_shares(src_reg, src_bits, size);

	if(0 != (*the_ext_lib.x_mp_skew_recomp)(spdz_gf2n_ext_handle, n, (const uint64_t*)src_bits.data(), (uint64_t *)&dest_ring))
	{
		cerr << "Processor::MP_Skew_Bit_Recomp_Ext extension library mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}
#endif

#if defined(EXTENDED_SPDZ_GFP) || defined(EXTENDED_SPDZ_GF2N)

#define LOAD_LIB_METHOD(Name,Proc)	\
if(0 != load_extension_method(Name, (void**)(&Proc), x_lib_handle)) { dlclose(x_lib_handle); abort(); }

spdz_ext_ifc::spdz_ext_ifc()
{
	x_lib_handle = NULL;
	*(void**)(&x_init) = NULL;
	*(void**)(&x_term) = NULL;
	*(void**)(&x_offline) = NULL;
	*(void**)(&x_opens) = NULL;
	*(void**)(&x_triple) = NULL;
	*(void**)(&x_verify) = NULL;
	*(void**)(&x_input) = NULL;
	*(void**)(&x_mult) = NULL;
	*(void**)(&x_mix_add) = NULL;
	*(void**)(&x_mix_sub_scalar) = NULL;
	*(void**)(&x_mix_sub_share) = NULL;
	*(void**)(&x_closes) = NULL;
	*(void**)(&x_bit) = NULL;
	*(void**)(&x_inverse) = NULL;
#if defined(EXTENDED_SPDZ_Z2N)
	*(void**)(&x_skew_decomp) = NULL;
	*(void**)(&x_skew_recomp) = NULL;
	*(void**)(&x_skew_inject) = NULL;
	*(void**)(&x_mp_closes) = NULL;
	*(void**)(&x_mp_opens) = NULL;
	*(void**)(&x_mp_mix_add) = NULL;
	*(void**)(&x_mp_mix_sub_share) = NULL;
	*(void**)(&x_mp_mix_sub_scalar) = NULL;
	*(void**)(&x_mp_mult) = NULL;
	*(void**)(&x_mp_skew_decomp) = NULL;
	*(void**)(&x_mp_skew_recomp) = NULL;
	*(void**)(&x_mp_skew_inject) = NULL;
#endif

	//get the SPDZ-2 extension library for env-var
	const char * spdz_ext_lib = getenv("SPDZ_EXT_LIB");
	if(NULL == spdz_ext_lib)
	{
		cerr << "SPDZ extension library defined not set" << endl;
		abort();
	}
	cout << "set extension library " << spdz_ext_lib << endl;

	//verify the SPDZ-2 extension library exists
	struct stat st;
	if(0 != stat(spdz_ext_lib, &st))
	{
		cerr << "failed to find extension library " << spdz_ext_lib << endl;
		abort();
	}
	cout << "found extension library " << spdz_ext_lib << endl;

	//load the SPDZ-2 extension library
	x_lib_handle = dlopen(spdz_ext_lib, RTLD_NOW);
	if(NULL == x_lib_handle)
	{
		const char * dlopen_err_msg = dlerror();
		cerr << "failed to load extension library [" << ((NULL != dlopen_err_msg)? dlopen_err_msg: "") << "]" << endl;
		abort();
	}

	//loading the SPDZ-2 extension library methods
	LOAD_LIB_METHOD("init", x_init)
	LOAD_LIB_METHOD("term", x_term)
	LOAD_LIB_METHOD("offline", x_offline)
	LOAD_LIB_METHOD("opens", x_opens)
	LOAD_LIB_METHOD("triple", x_triple)
	LOAD_LIB_METHOD("verify", x_verify)
	LOAD_LIB_METHOD("input", x_input)
	LOAD_LIB_METHOD("mult", x_mult)
	LOAD_LIB_METHOD("mix_add", x_mix_add)
	LOAD_LIB_METHOD("mix_sub_scalar", x_mix_sub_scalar)
	LOAD_LIB_METHOD("mix_sub_share", x_mix_sub_share)
	LOAD_LIB_METHOD("mix_mul", x_mix_mul)
	LOAD_LIB_METHOD("adds", x_adds)
	LOAD_LIB_METHOD("subs", x_subs)
	LOAD_LIB_METHOD("closes", x_closes)
	LOAD_LIB_METHOD("bit", x_bit)
	LOAD_LIB_METHOD("inverse", x_inverse)
#if defined(EXTENDED_SPDZ_Z2N)
	LOAD_LIB_METHOD("skew_decomp", x_skew_decomp)
	LOAD_LIB_METHOD("skew_recomp", x_skew_recomp)
	LOAD_LIB_METHOD("skew_inject", x_skew_inject)
	LOAD_LIB_METHOD("mp_closes", x_mp_closes)
	LOAD_LIB_METHOD("mp_opens", x_mp_opens)
	LOAD_LIB_METHOD("mp_adds", x_mp_adds)
	LOAD_LIB_METHOD("mp_mix_add", x_mp_mix_add)
	LOAD_LIB_METHOD("mp_subs", x_mp_subs)
	LOAD_LIB_METHOD("mp_mix_sub_share", x_mp_mix_sub_share)
	LOAD_LIB_METHOD("mp_mix_sub_scalar", x_mp_mix_sub_scalar)
	LOAD_LIB_METHOD("mp_mix_mul", x_mp_mix_mul)
	LOAD_LIB_METHOD("mp_mult", x_mp_mult)
	LOAD_LIB_METHOD("mp_skew_decomp", x_mp_skew_decomp)
	LOAD_LIB_METHOD("mp_skew_recomp", x_mp_skew_recomp)
	LOAD_LIB_METHOD("mp_skew_inject", x_mp_skew_inject)
#endif
}

spdz_ext_ifc::~spdz_ext_ifc()
{
	dlclose(x_lib_handle);
}

int spdz_ext_ifc::load_extension_method(const char * method_name, void ** proc_addr, void * libhandle)
{
	*proc_addr = dlsym(libhandle, method_name);
	const char * dlsym_error = dlerror();
	if(NULL != dlsym_error || NULL == *proc_addr)
	{
		cerr << "failed to load " << method_name << " extension [" << ((NULL != dlsym_error)? dlsym_error: "") << "]" << endl;
		return -1;
	}
	return 0;
}

#endif

template void Processor::POpen_Start(const vector<int>& reg,const Player& P,MAC_Check<gf2n>& MC,int size);
template void Processor::POpen_Start(const vector<int>& reg,const Player& P,MAC_Check<gfp>& MC,int size);
template void Processor::POpen_Stop(const vector<int>& reg,const Player& P,MAC_Check<gf2n>& MC,int size);
template void Processor::POpen_Stop(const vector<int>& reg,const Player& P,MAC_Check<gfp>& MC,int size);
template void Processor::POpen(const vector<int>& reg,const Player& P,MAC_Check<gf2n>& MC,int size);
template void Processor::POpen(const vector<int>& reg,const Player& P,MAC_Check<gfp>& MC,int size);
template void Processor::read_socket_private<gfp>(int client_id, const vector<int>& registers, bool send_macs);
template void Processor::read_socket_vector<gfp>(int client_id, const vector<int>& registers);
template void Processor::read_shares_from_file<gfp>(int start_file_pos, int end_file_pos_register, const vector<int>& data_registers);
template void Processor::write_shares_to_file<gfp>(const vector<int>& data_registers);
