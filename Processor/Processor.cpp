// (C) 2017 University of Bristol. See License.txt


#include "Processor/Processor.h"
#include "Networking/STS.h"
#include "Auth/MAC_Check.h"

#include "Auth/fake-stuff.h"
#include <sodium.h>
#include <string>

#if defined(EXTENDED_SPDZ)
#include <sys/stat.h>
#include <dlfcn.h>
#include <list>

spdz_ext_ifc the_ext_lib;
#endif

Processor::Processor(int thread_num,Data_Files& DataF,Player& P,
        MAC_Check<gf2n>& MC2,MAC_Check<gfp>& MCp,Machine& machine,
        const Program& program)
: thread_num(thread_num),DataF(DataF),P(P),MC2(MC2),MCp(MCp),machine(machine),
  input2(*this,MC2),inputp(*this,MCp),privateOutput2(*this),privateOutputp(*this),sent(0),rounds(0),
  external_clients(ExternalClients(P.my_num(), DataF.prep_data_dir)),binary_file_io(Binary_File_IO())
#if defined(EXTENDED_SPDZ)
  , po_shares(NULL), po_opens(NULL), po_size(0)
  , pi_inputs(NULL), pi_size(0)
  , pm_shares(NULL), pm_products(NULL), pm_size(0)
  , go_shares(NULL), go_opens(NULL), go_size(0)
  , gi_inputs(NULL), gi_size(0)
  , gm_shares(NULL), gm_products(NULL), gm_size(0)
#endif
{
  reset(program,0);

  public_input.open(get_filename("Programs/Public-Input/",false).c_str());
  private_input.open(get_filename("Player-Data/Private-Input-",true).c_str());
  public_output.open(get_filename("Player-Data/Public-Output-",true).c_str(), ios_base::out);
  private_output.open(get_filename("Player-Data/Private-Output-",true).c_str(), ios_base::out);

#if defined(EXTENDED_SPDZ)
    spdz_gfp_ext_handle = NULL;
	cout << "Processor " << thread_num << " SPDZ GFP extension library initializing." << endl;
	if(0 != (*the_ext_lib.ext_init)(&spdz_gfp_ext_handle, P.my_num(), P.num_players(), "gfp61", 100000, 1000, 100000))
	{
		cerr << "SPDZ GFP extension library initialization failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	cout << "SPDZ GFP extension library initialized." << endl;

    spdz_gf2n_ext_handle = NULL;
	cout << "SPDZ GF2N extension library initializing." << endl;
	if(0 != (*the_ext_lib.ext_init)(&spdz_gf2n_ext_handle, P.my_num(), P.num_players(), "gf2n40", 200, 200, 200))
	{
		cerr << "SPDZ GF2N extension library initialization failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	cout << "SPDZ GF2N extension library initialized." << endl;

	alloc_po_mpz(5000);
#endif
}

Processor::~Processor()
{
  cerr << "Sent " << sent << " elements in " << rounds << " rounds" << endl;
#if defined(EXTENDED_SPDZ)
	(*the_ext_lib.ext_term)(spdz_gfp_ext_handle);
	(*the_ext_lib.ext_term)(spdz_gf2n_ext_handle);
	dlclose(the_ext_lib.ext_lib_handle);

	free_po_mpz();
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
// Encryption is enabled if key material (for DH Auth Encryption and/or STS protocol) has been already setup.
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

  // Apply DH Auth encryption if session keys have been created.
  map<int,octet*>::iterator it = external_clients.symmetric_client_keys.find(socket_id);
  if (it != external_clients.symmetric_client_keys.end()) {
    socket_stream.encrypt(it->second);
  }

  // Apply STS commsec encryption if session keys have been created.
  try {
    maybe_encrypt_sequence(socket_id);
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
  maybe_decrypt_sequence(client_id);
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
  maybe_decrypt_sequence(client_id);
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
  maybe_decrypt_sequence(client_id);

  map<int,octet*>::iterator it = external_clients.symmetric_client_keys.find(client_id);
  if (it != external_clients.symmetric_client_keys.end())
  {
    socket_stream.decrypt(it->second);
  }
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

  external_clients.generate_session_key_for_client(client_id, client_public_key);  
}

void Processor::init_secure_socket_internal(int client_id, const vector<int>& registers) {
  external_clients.symmetric_client_commsec_send_keys.erase(client_id);
  external_clients.symmetric_client_commsec_recv_keys.erase(client_id);
  unsigned char client_public_bytes[crypto_sign_PUBLICKEYBYTES];
  sts_msg1_t m1;
  sts_msg2_t m2;
  sts_msg3_t m3;

  external_clients.load_server_keys_once();
  external_clients.require_ed25519_keys();

  // Validate inputs and state
  if(registers.size() != 8) {
      throw "Invalid call to init_secure_socket.";
  }
  if (client_id >= (int)external_clients.external_client_sockets.size())
  {
    cerr << "No socket connection exists for client id " << client_id << endl;
    throw "No socket connection exists for client";
  }

  // Extract client long term public key into bytes
  vector<int> client_public_key (registers.size(), 0);
  for(unsigned int i = 0; i < registers.size(); i++) {
    client_public_key[i] = (int&)get_Ci_ref(registers[i]);
  }
  external_clients.curve25519_ints_to_bytes(client_public_bytes,  client_public_key);

  // Start Station to Station Protocol
  STS ke(client_public_bytes, external_clients.server_publickey_ed25519, external_clients.server_secretkey_ed25519);
  m1 = ke.send_msg1();
  socket_stream.reset_write_head();
  socket_stream.append(m1.bytes, sizeof m1.bytes);
  socket_stream.Send(external_clients.external_client_sockets[client_id]);
  socket_stream.ReceiveExpected(external_clients.external_client_sockets[client_id],
                                96);
  socket_stream.consume(m2.pubkey, sizeof m2.pubkey);
  socket_stream.consume(m2.sig, sizeof m2.sig);
  m3 = ke.recv_msg2(m2);
  socket_stream.reset_write_head();
  socket_stream.append(m3.bytes, sizeof m3.bytes);
  socket_stream.Send(external_clients.external_client_sockets[client_id]);

  // Use results of STS to generate send and receive keys.
  vector<unsigned char> sendKey = ke.derive_secret(crypto_secretbox_KEYBYTES);
  vector<unsigned char> recvKey = ke.derive_secret(crypto_secretbox_KEYBYTES);
  external_clients.symmetric_client_commsec_send_keys[client_id] = make_pair(sendKey,0);
  external_clients.symmetric_client_commsec_recv_keys[client_id] = make_pair(recvKey,0);
}

void Processor::init_secure_socket(int client_id, const vector<int>& registers) {

  try {
      init_secure_socket_internal(client_id, registers);
  } catch (char const *e) {
      cerr << "STS initiator role failed with: " << e << endl;
      throw Processor_Error("STS initiator failed");
  }
}

void Processor::resp_secure_socket(int client_id, const vector<int>& registers) {
  try {
      resp_secure_socket_internal(client_id, registers);
  } catch (char const *e) {
      cerr << "STS responder role failed with: " << e << endl;
      throw Processor_Error("STS responder failed");
  }
}

void Processor::resp_secure_socket_internal(int client_id, const vector<int>& registers) {
  external_clients.symmetric_client_commsec_send_keys.erase(client_id);
  external_clients.symmetric_client_commsec_recv_keys.erase(client_id);
  unsigned char client_public_bytes[crypto_sign_PUBLICKEYBYTES];
  sts_msg1_t m1;
  sts_msg2_t m2;
  sts_msg3_t m3;

  external_clients.load_server_keys_once();
  external_clients.require_ed25519_keys();

  // Validate inputs and state
  if(registers.size() != 8) {
      throw "Invalid call to init_secure_socket.";
  }
  if (client_id >= (int)external_clients.external_client_sockets.size())
  {
    cerr << "No socket connection exists for client id " << client_id << endl;
    throw "No socket connection exists for client";
  }
  vector<int> client_public_key (registers.size(), 0);
  for(unsigned int i = 0; i < registers.size(); i++) {
    client_public_key[i] = (int&)get_Ci_ref(registers[i]);
  }
  external_clients.curve25519_ints_to_bytes(client_public_bytes,  client_public_key);

  // Start Station to Station Protocol for the responder
  STS ke(client_public_bytes, external_clients.server_publickey_ed25519, external_clients.server_secretkey_ed25519);
  socket_stream.reset_read_head();
  socket_stream.ReceiveExpected(external_clients.external_client_sockets[client_id],
                                32);
  socket_stream.consume(m1.bytes, sizeof m1.bytes);
  m2 = ke.recv_msg1(m1);
  socket_stream.reset_write_head();
  socket_stream.append(m2.pubkey, sizeof m2.pubkey);
  socket_stream.append(m2.sig, sizeof m2.sig);
  socket_stream.Send(external_clients.external_client_sockets[client_id]);

  socket_stream.ReceiveExpected(external_clients.external_client_sockets[client_id],
                                64);
  socket_stream.consume(m3.bytes, sizeof m3.bytes);
  ke.recv_msg3(m3);

  // Use results of STS to generate send and receive keys.
  vector<unsigned char> recvKey = ke.derive_secret(crypto_secretbox_KEYBYTES);
  vector<unsigned char> sendKey = ke.derive_secret(crypto_secretbox_KEYBYTES);
  external_clients.symmetric_client_commsec_recv_keys[client_id] = make_pair(recvKey,0);
  external_clients.symmetric_client_commsec_send_keys[client_id] = make_pair(sendKey,0);
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

void Processor::maybe_decrypt_sequence(int client_id)
{
  map<int, pair<vector<octet>,uint64_t> >::iterator it_cs = external_clients.symmetric_client_commsec_recv_keys.find(client_id);
  if (it_cs != external_clients.symmetric_client_commsec_recv_keys.end())
  {
    socket_stream.decrypt_sequence(&it_cs->second.first[0], it_cs->second.second);
    it_cs->second.second++;
  }
}

void Processor::maybe_encrypt_sequence(int client_id)
{
  map<int, pair<vector<octet>,uint64_t> >::iterator it_cs = external_clients.symmetric_client_commsec_send_keys.find(client_id);
  if (it_cs != external_clients.symmetric_client_commsec_send_keys.end())
  {
    socket_stream.encrypt_sequence(&it_cs->second.first[0], it_cs->second.second);
    it_cs->second.second++;
  }
}

#if defined(EXTENDED_SPDZ)

void Processor::POpen_Start_Ext_64(const vector<int>& reg, int size)
{
	int sz=reg.size();

	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(reg, Sh_PO, size);

	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sz*size);

	//the share values are saved as mpz
	if(Sh_PO.size() > po_size)
	{
		free_po_mpz();
		alloc_po_mpz(Sh_PO.size());
	}
	PShares2mpz(Sh_PO, po_shares);

	//the extension library is given the shares' values and returns opens' values
	if(0 != (*the_ext_lib.ext_start_open)(spdz_gfp_ext_handle, Sh_PO.size(), po_shares, po_opens, 1))
	{
		cerr << "Processor::POpen_Start_Ext_64 extension library start_open failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}

void Processor::POpen_Stop_Ext_64(const vector<int>& reg,int size)
{
	vector<gfp>& PO = get_PO<gfp>();
	vector<gfp>& C = get_C<gfp>();
	int sz=reg.size();
	PO.resize(sz*size);

	if(0 != (*the_ext_lib.ext_stop_open)(spdz_gfp_ext_handle))
	{
		cerr << "Processor::POpen_Stop_Ext_64 extension library stop_open failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	Pmpz2gfps(po_opens, PO);
	POpen_Stop_prep_opens(reg, PO, C, size);

	sent += reg.size() * size;
	rounds++;
}

void Processor::PTriple_Ext_64(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c)
{
	mpz_t ma, mb, mc;

	mpz_init(ma);
	mpz_init(mb);
	mpz_init(mc);

	if(0 != (*the_ext_lib.ext_triple)(spdz_gfp_ext_handle, &ma, &mb, &mc))
	{
		cerr << "Processor::PTriple_Ext_64 extension library triple failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	Pmpz2share(&ma, a);
	Pmpz2share(&mb, b);
	Pmpz2share(&mc, c);

	mpz_clear(ma);
	mpz_clear(mb);
	mpz_clear(mc);
}

void Processor::PInput_Ext_64(Share<gfp>& input_value, const int input_party_id)
{
	mpz_t mpz_input_value;
	mpz_init(mpz_input_value);
	if(0 != (*the_ext_lib.ext_input)(spdz_gfp_ext_handle, input_party_id, &mpz_input_value))
	{
		cerr << "Processor::PInput_Ext_64 extension library input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	Pmpz2share(&mpz_input_value, input_value);
	mpz_clear(mpz_input_value);
}

void Processor::PInput_Start_Ext_64(int player, int n_inputs)
{
	alloc_pi_mpz(n_inputs);
	if(0 != (*the_ext_lib.ext_start_input)(spdz_gfp_ext_handle, player, n_inputs, pi_inputs))
	{
		cerr << "Processor::PInput_Start_Ext_64 extension library start input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}

void Processor::PInput_Stop_Ext_64(int /*player*/, vector<int> targets)
{
	if(0 != (*the_ext_lib.ext_stop_input)(spdz_gfp_ext_handle))
	{
		cerr << "Processor::PInput_Stop_Ext_64 extension library stop input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(targets.size() == pi_size)
	{
		for(size_t i = 0; i < pi_size; ++i)
		{
			Share<gfp>& share = get_S_ref<gfp>(targets[i]);
			Pmpz2share(pi_inputs + i, share);
		}
	}
	else
	{
		cerr << "Processor::PInput_Stop_Ext_64 extension library stop input mismatched number of inputs " << targets.size() << "/" << pi_size << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	free_pi_mpz();
}

void Processor::PMult_Start_Ext_64(const vector<int>& reg, int size)
{
	int sz=reg.size();

	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(reg, Sh_PO, size);

	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sz*size);

	//the share values are saved as mpz
	alloc_pm_mpz(Sh_PO.size());
	PShares2mpz(Sh_PO, pm_shares);

	if(0 != (*the_ext_lib.ext_start_mult)(spdz_gfp_ext_handle, pm_size, pm_shares, pm_products, 1))
	{
		cerr << "Processor::PMult_Start_Ext_64 extension library start_mult failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	else
	{
		cout << "Processor::PMult_Start_Ext_64 extension library start_mult launched." << endl;
	}
}

void Processor::PMult_Stop_Ext_64(const vector<int>& reg, int size)
{
	if(0 != (*the_ext_lib.ext_stop_mult)(spdz_gfp_ext_handle))
	{
		cerr << "Processor::PMult_Stop_Ext_64 extension library stop_mult failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	PMult_Stop_prep_products(reg, size);
	free_pm_mpz();

	sent += reg.size() * size;
	rounds++;
}

void Processor::PMult_Stop_prep_products(const vector<int>& reg, int size)
{
	if (size>1)
	{
		size_t product_idx = 0;
		for (typename vector<int>::const_iterator reg_it=reg.begin(); reg_it!=reg.end(); reg_it++)
		{
			vector<Share<gfp> >::iterator insert_point=get_S<gfp>().begin()+*reg_it;
			for(int i = 0; i < size; ++i)
			{
				Pmpz2share(pm_products + (product_idx++), *(insert_point + i));
			}
		}
	}
	else
	{
		int sz=reg.size();
		for(int i = 0; i < sz; ++i)
		{
			Pmpz2share(pm_products + i, get_S_ref<gfp>(reg[i]));
		}
	}
}

void Processor::PAddm_Ext_64(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	mpz_t share_value, arg;
	mpz_init(share_value);
	mpz_init(arg);
	to_bigint(*((mpz_class*)(&share_value)), a.get_share());
	to_bigint(*((mpz_class*)(&arg)), b);
	if(0 == (*the_ext_lib.ext_mix_add)(spdz_gfp_ext_handle, &share_value, &arg))
	{
		Pmpz2share(&share_value, c);
	}
	else
	{
		cerr << "Processor::PAddm_Ext_64 extension library mix_add failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(share_value);
	mpz_clear(arg);
}

void Processor::PSubml_Ext_64(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	mpz_t share_value, arg;
	mpz_init(share_value);
	mpz_init(arg);
	to_bigint(*((mpz_class*)(&share_value)), a.get_share());
	to_bigint(*((mpz_class*)(&arg)), b);
	if(0 == (*the_ext_lib.ext_mix_sub_scalar)(spdz_gfp_ext_handle, &share_value, &arg))
	{
		Pmpz2share(&share_value, c);
	}
	else
	{
		cerr << "Processor::PSubml_Ext_64 extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(share_value);
	mpz_clear(arg);
}

void Processor::PSubmr_Ext_64(gfp& a, Share<gfp>& b, Share<gfp>& c)
{
	mpz_t share_value, arg;
	mpz_init(share_value);
	mpz_init(arg);
	to_bigint(*((mpz_class*)(&share_value)), b.get_share());
	to_bigint(*((mpz_class*)(&arg)), a);
	if(0 == (*the_ext_lib.ext_mix_sub_share)(spdz_gfp_ext_handle, &arg, &share_value))
	{
		Pmpz2share(&share_value, c);
	}
	else
	{
		cerr << "Processor::PSubmr_Ext_64 extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(share_value);
	mpz_clear(arg);
}

void Processor::PLdsi_Ext_64(gfp& value, Share<gfp>& share)
{
	mpz_t mpz_value, mpz_share;
	mpz_init(mpz_value);
	mpz_init(mpz_share);
	to_bigint(*((mpz_class*)(&mpz_value)), value);
	if(0 == (*the_ext_lib.ext_share_immediate)(spdz_gfp_ext_handle, &mpz_value, &mpz_share))
	{
		Pmpz2share(&mpz_share, share);
	}
	else
	{
		cerr << "Processor::PLdsi_Ext_64 extension library share_immediate failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(mpz_value);
	mpz_clear(mpz_share);
}

void Processor::PBit_Ext_64(Share<gfp>& share)
{
	mpz_t mpz_share;
	mpz_init(mpz_share);
	if(0 == (*the_ext_lib.ext_bit)(spdz_gfp_ext_handle, &mpz_share))
	{
		Pmpz2share(&mpz_share, share);
	}
	else
	{
		cerr << "Processor::PBit_Ext_64 extension library bit failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(mpz_share);
}

void Processor::PInverse_Ext_64(Share<gfp>& share_value, Share<gfp>& share_inverse)
{
	mpz_t mpz_share_value, mpz_share_inverse;
	mpz_init(mpz_share_value);
	mpz_init(mpz_share_inverse);
	if(0 == (*the_ext_lib.ext_inverse)(spdz_gfp_ext_handle, &mpz_share_value, &mpz_share_inverse))
	{
		Pmpz2share(&mpz_share_value, share_value);
		Pmpz2share(&mpz_share_inverse, share_inverse);
	}
	else
	{
		cerr << "Processor::PBit_Ext_64 extension library inverse failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(mpz_share_value);
	mpz_clear(mpz_share_inverse);
}

void Processor::PShares2mpz(const vector< Share<gfp> >& shares, mpz_t * share_values)
{
	size_t count = shares.size();
	for(size_t i = 0; i < count; i++)
	{
		to_bigint(*((mpz_class*)(share_values + i)), shares[i].get_share());
	}
}

void Processor::Pmpz2gfps(const mpz_t * mpz_values, vector<gfp>& gfps)
{
	size_t count = gfps.size();
	for(size_t i = 0; i < count; i++)
	{
		to_gfp(gfps[i], *((mpz_class*)(mpz_values + i)));
	}
}

void Processor::Pmpz2share(const mpz_t * mpzv, Share<gfp> & shv)
{
	gfp mac, value;
	to_gfp(value, *((mpz_class*)mpzv));
	mac.mul(MCp.get_alphai(), value);
	shv.set_share(value);
	shv.set_mac(mac);
}

void Processor::GOpen_Start_Ext_64(const vector<int>& reg,int size)
{
	int sz=reg.size();

	vector< Share<gf2n> >& Sh_PO = get_Sh_PO<gf2n>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(reg, Sh_PO, size);

	vector<gf2n>& PO = get_PO<gf2n>();
	PO.resize(sz*size);

	alloc_go_mpz(Sh_PO.size());
	GShares2mpz(Sh_PO, go_shares);

	if(0 != (*the_ext_lib.ext_start_open)(spdz_gf2n_ext_handle, go_size, go_shares, go_opens, 1))
	{
		cerr << "Processor::GOpen_Start_Ext_64 extension library start_open failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	else
	{
		cout << "Processor::GOpen_Start_Ext_64 extension start open launched." << endl;
	}
}

void Processor::GOpen_Stop_Ext_64(const vector<int>& reg,int size)
{
	vector<gf2n>& PO = get_PO<gf2n>();
	vector<gf2n>& C = get_C<gf2n>();
	int sz=reg.size();
	PO.resize(sz*size);

	if(0 != (*the_ext_lib.ext_stop_open)(spdz_gf2n_ext_handle))
	{
		cerr << "Processor::GOpen_Stop_Ext_64 extension library stop_open failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	Gmpz2gf2ns(go_opens, PO);
	free_go_mpz();
	POpen_Stop_prep_opens(reg, PO, C, size);

	sent += reg.size() * size;
	rounds++;
}

void Processor::GTriple_Ext_64(Share<gf2n>& a, Share<gf2n>& b, Share<gf2n>& c)
{
	mpz_t mpza, mpzb, mpzc;
	mpz_init(mpza);
	mpz_init(mpzb);
	mpz_init(mpzc);
	if(0 != (*the_ext_lib.ext_triple)(spdz_gf2n_ext_handle, &mpza, &mpzb, &mpzc))
	{
		cerr << "Processor::GTriple_Ext_64 extension library triple failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	Gmpz2share(&mpza, a);
	Gmpz2share(&mpzb, b);
	Gmpz2share(&mpzc, c);
	mpz_clear(mpza);
	mpz_clear(mpzb);
	mpz_clear(mpzc);
}

void Processor::GInput_Ext_64(Share<gf2n>& input_value, const int input_party_id)
{
	mpz_t mpzv;
	mpz_init(mpzv);
	if(0 != (*the_ext_lib.ext_input)(spdz_gf2n_ext_handle, input_party_id, &mpzv))
	{
		cerr << "Processor::GInput_Ext_64 extension library input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	Gmpz2share(&mpzv, input_value);
	mpz_clear(mpzv);
}

void Processor::GInput_Start_Ext_64(int player, int n_inputs)
{
	alloc_gi_mpz(n_inputs);
	if(0 != (*the_ext_lib.ext_start_input)(spdz_gf2n_ext_handle, player, n_inputs, &gi_inputs[0]))
	{
		cerr << "Processor::GInput_Start_Ext_64 extension library start input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}

void Processor::GInput_Stop_Ext_64(int /*player*/, vector<int> targets)
{
	if(0 != (*the_ext_lib.ext_stop_input)(spdz_gf2n_ext_handle))
	{
		cerr << "Processor::GInput_Stop_Ext_64 extension library stop input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(targets.size() == gi_size)
	{
		for(size_t i = 0; i < gi_size; ++i)
		{
			Share<gf2n>& share = get_S_ref<gf2n>(targets[i]);
			Gmpz2share(gi_inputs + i, share);
		}
	}
	else
	{
		cerr << "Processor::GInput_Stop_Ext_64 extension library stop input returned mismatched number of inputs " << targets.size() << "/" << gi_size << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	free_gi_mpz();
}

void Processor::GMult_Start_Ext_64(const vector<int>& reg, int size)
{
	int sz=reg.size();

	vector< Share<gf2n> >& Sh_PO = get_Sh_PO<gf2n>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	prep_shares(reg, Sh_PO, size);

	vector<gf2n>& PO = get_PO<gf2n>();
	PO.resize(sz*size);

	//the share values are saved as mpz
	alloc_gm_mpz(Sh_PO.size());
	GShares2mpz(Sh_PO, gm_shares);

	if(0 != (*the_ext_lib.ext_start_mult)(spdz_gf2n_ext_handle, gm_size, gm_shares, gm_products, 1))
	{
		cerr << "Processor::GMult_Start_Ext_64 extension library start_mult failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	else
	{
		cout << "Processor::GMult_Start_Ext_64 extension library start_mult launched." << endl;
	}
}

void Processor::GMult_Stop_Ext_64(const vector<int>& reg, int size)
{
	if(0 != (*the_ext_lib.ext_stop_mult)(spdz_gf2n_ext_handle))
	{
		cerr << "Processor::GMult_Stop_Ext_64 extension library stop_mult failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	GMult_Stop_prep_products(reg, size);
	free_gm_mpz();

	sent += reg.size() * size;
	rounds++;
}

void Processor::GMult_Stop_prep_products(const vector<int>& reg, int size)
{
	if (size>1)
	{
		size_t product_idx = 0;
		for (typename vector<int>::const_iterator reg_it=reg.begin(); reg_it!=reg.end(); reg_it++)
		{
			vector<Share<gf2n> >::iterator insert_point=get_S<gf2n>().begin()+*reg_it;
			for(int i = 0; i < size; ++i)
			{
				Gmpz2share(gm_products + (product_idx++), *(insert_point + i));
			}
		}
	}
	else
	{
		int sz=reg.size();
		for(int i = 0; i < sz; ++i)
		{
			Gmpz2share(gm_products + i, get_S_ref<gf2n>(reg[i]));
		}
	}
}

void Processor::GAddm_Ext_64(Share<gf2n>& a, gf2n& b, Share<gf2n>& c)
{
	mpz_t share_value, arg;
	mpz_init(share_value);
	mpz_init(arg);
	mpz_set_ui(share_value, a.get_share().get_word());
	mpz_set_ui(arg, b.get_word());
	if(0 == (*the_ext_lib.ext_mix_add)(spdz_gf2n_ext_handle, &share_value, &arg))
	{
		Gmpz2share(&share_value, c);
	}
	else
	{
		cerr << "Processor::GAddm_Ext_64 extension library mix_add failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(share_value);
	mpz_clear(arg);
}

void Processor::GSubml_Ext_64(Share<gf2n>& a, gf2n& b, Share<gf2n>& c)
{
	mpz_t share_value, arg;
	mpz_init(share_value);
	mpz_init(arg);
	mpz_set_ui(share_value, a.get_share().get_word());
	mpz_set_ui(arg, b.get_word());
	if(0 == (*the_ext_lib.ext_mix_sub_scalar)(spdz_gf2n_ext_handle, &share_value, &arg))
	{
		Gmpz2share(&share_value, c);
	}
	else
	{
		cerr << "Processor::GSubml_Ext_64 extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(share_value);
	mpz_clear(arg);
}

void Processor::GSubmr_Ext_64(gf2n& a, Share<gf2n>& b, Share<gf2n>& c)
{
	mpz_t share_value, arg;
	mpz_init(share_value);
	mpz_init(arg);
	mpz_set_ui(share_value, b.get_share().get_word());
	mpz_set_ui(arg, a.get_word());
	if(0 == (*the_ext_lib.ext_mix_sub_share)(spdz_gf2n_ext_handle, &arg, &share_value))
	{
		Gmpz2share(&share_value, c);
	}
	else
	{
		cerr << "Processor::GSubmr_Ext_64 extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(share_value);
	mpz_clear(arg);
}

void Processor::GLdsi_Ext_64(gf2n& value, Share<gf2n>& share)
{
	mpz_t mpz_value, mpz_share;
	mpz_init(mpz_value);
	mpz_init(mpz_share);
	mpz_set_ui(mpz_value, value.get_word());
	if(0 == (*the_ext_lib.ext_share_immediate)(spdz_gf2n_ext_handle, &mpz_value, &mpz_share))
	{
		Gmpz2share(&mpz_share, share);
	}
	else
	{
		cerr << "Processor::GLdsi_Ext_64 extension library share_immediate failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(mpz_value);
	mpz_clear(mpz_share);
}

void Processor::GBit_Ext_64(Share<gf2n>& share)
{
	mpz_t mpz_share;
	mpz_init(mpz_share);
	if(0 == (*the_ext_lib.ext_bit)(spdz_gf2n_ext_handle, &mpz_share))
	{
		Gmpz2share(&mpz_share, share);
	}
	else
	{
		cerr << "Processor::GBit_Ext_64 extension library bit failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(mpz_share);
}

void Processor::GInverse_Ext_64(Share<gf2n>& share_value, Share<gf2n>& share_inverse)
{
	mpz_t mpz_share_value, mpz_share_inverse;
	mpz_init(mpz_share_value);
	mpz_init(mpz_share_inverse);
	if(0 == (*the_ext_lib.ext_inverse)(spdz_gf2n_ext_handle, &mpz_share_value, &mpz_share_inverse))
	{
		Gmpz2share(&mpz_share_value, share_value);
		Gmpz2share(&mpz_share_inverse, share_inverse);
	}
	else
	{
		cerr << "Processor::GInverse_Ext_64 extension library inverse failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	mpz_clear(mpz_share_value);
	mpz_clear(mpz_share_inverse);
}

void Processor::GShares2mpz(const vector< Share<gf2n> >& shares, mpz_t * share_values)
{
	size_t count = shares.size();
	for(size_t i = 0; i < count; i++)
	{
		mpz_set_ui(share_values[i], shares[i].get_share().get_word());
	}
}

void Processor::Gmpz2gf2ns(const mpz_t * mpz_values, vector<gf2n>& gf2ns)
{
	size_t count = gf2ns.size();
	for(size_t i = 0; i < count; i++)
	{
		gf2ns[i].assign(mpz_get_ui(mpz_values[i]));
	}
}

void Processor::Gmpz2share(const mpz_t * mpzv, Share<gf2n> & shv)
{
	gf2n mac, value;
	value.assign(mpz_get_ui(*mpzv));
	mac.mul(MC2.get_alphai(), value);
	shv.set_share(value);
	shv.set_mac(mac);

}

#define LOAD_LIB_METHOD(Name,Proc)	\
if(0 != load_extension_method(Name, (void**)(&Proc), ext_lib_handle)) { dlclose(ext_lib_handle); abort(); }

spdz_ext_ifc::spdz_ext_ifc()
{
	ext_lib_handle = NULL;
	*(void**)(&ext_init) = NULL;
	*(void**)(&ext_term) = NULL;
	*(void**)(&ext_offline) = NULL;
	*(void**)(&ext_start_open) = NULL;
	*(void**)(&ext_stop_open) = NULL;
	*(void**)(&ext_triple) = NULL;
	*(void**)(&ext_input) = NULL;
	*(void**)(&ext_start_verify) = NULL;
	*(void**)(&ext_stop_verify) = NULL;
	*(void**)(&ext_start_input) = NULL;
	*(void**)(&ext_stop_input) = NULL;
	*(void**)(&ext_start_mult) = NULL;
	*(void**)(&ext_stop_mult) = NULL;
	*(void**)(&ext_mix_add) = NULL;
	*(void**)(&ext_mix_sub_scalar) = NULL;
	*(void**)(&ext_mix_sub_share) = NULL;
	*(void**)(&ext_bit) = NULL;
	*(void**)(&ext_inverse) = NULL;

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
	ext_lib_handle = dlopen(spdz_ext_lib, RTLD_NOW);
	if(NULL == ext_lib_handle)
	{
		const char * dlopen_err_msg = dlerror();
		cerr << "failed to load extension library [" << ((NULL != dlopen_err_msg)? dlopen_err_msg: "") << "]" << endl;
		abort();
	}

	//loading the SPDZ-2 extension library methods
	LOAD_LIB_METHOD("init",ext_init)
	LOAD_LIB_METHOD("term",ext_term)
	LOAD_LIB_METHOD("offline",ext_offline)
	LOAD_LIB_METHOD("start_open",ext_start_open)
	LOAD_LIB_METHOD("stop_open",ext_stop_open)
	LOAD_LIB_METHOD("triple",ext_triple)
	LOAD_LIB_METHOD("input",ext_input)
	LOAD_LIB_METHOD("start_verify",ext_start_verify)
	LOAD_LIB_METHOD("stop_verify",ext_stop_verify)
	LOAD_LIB_METHOD("start_input",ext_start_input)
	LOAD_LIB_METHOD("stop_input",ext_stop_input)
	LOAD_LIB_METHOD("start_mult",ext_start_mult)
	LOAD_LIB_METHOD("stop_mult",ext_stop_mult)
	LOAD_LIB_METHOD("mix_add",ext_mix_add)
	LOAD_LIB_METHOD("mix_sub_scalar",ext_mix_sub_scalar)
	LOAD_LIB_METHOD("mix_sub_share",ext_mix_sub_share)
	LOAD_LIB_METHOD("share_immediate",ext_share_immediate)
	LOAD_LIB_METHOD("bit", ext_bit)
	LOAD_LIB_METHOD("inverse", ext_inverse)
}

spdz_ext_ifc::~spdz_ext_ifc()
{
	dlclose(ext_lib_handle);
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
template void Processor::read_socket_private<gfp>(int client_id, const vector<int>& registers, bool send_macs);
template void Processor::read_socket_vector<gfp>(int client_id, const vector<int>& registers);
template void Processor::read_shares_from_file<gfp>(int start_file_pos, int end_file_pos_register, const vector<int>& data_registers);
template void Processor::write_shares_to_file<gfp>(const vector<int>& data_registers);
