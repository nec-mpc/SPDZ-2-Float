// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


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
  private_input_filename(get_filename(PREP_DIR "Private-Input-",true)),
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
  private_input.open(private_input_filename.c_str());
  public_output.open(get_filename(PREP_DIR "Public-Output-",true).c_str(), ios_base::out);
  private_output.open(get_filename(PREP_DIR "Private-Output-",true).c_str(), ios_base::out);

#if defined(EXTENDED_SPDZ)
    spdz_gfp_ext_handle = NULL;
	cout << "Processor " << thread_num << " SPDZ GFP extension library initializing." << endl;
	if(0 != (*the_ext_lib.x_init)(&spdz_gfp_ext_handle, P.my_num(), P.num_players(), thread_num, "gfp127", 700000, 700000, 619200))
	{
		cerr << "SPDZ GFP extension library initialization failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
	cout << "SPDZ GFP extension library initialized." << endl;

	alloc_po_mpz(5000);
	alloc_pm_mpz(5000);
	mpz_init(mpz_share_aux);
	mpz_init(mpz_arg_aux);
#endif
}

Processor::~Processor()
{
  cerr << "Sent " << sent << " elements in " << rounds << " rounds" << endl;
#if defined(EXTENDED_SPDZ)
	(*the_ext_lib.x_term)(spdz_gfp_ext_handle);
	dlclose(the_ext_lib.x_lib_handle);

	free_po_mpz();
	free_pm_mpz();
	mpz_clear(mpz_share_aux);
	mpz_clear(mpz_arg_aux);
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

	//the share values are saved as mpz
	if(Sh_PO.size() > po_size)
	{
		free_po_mpz();
		alloc_po_mpz(Sh_PO.size());
	}
	PShares2mpz(Sh_PO, po_shares);

	//the extension library is given the shares' values and returns opens' values
	if(0 != (*the_ext_lib.x_opens)(spdz_gfp_ext_handle, Sh_PO.size(), po_shares, po_opens, 1))
	{
		cerr << "Processor::POpen_Ext extension library start_open failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	Pmpz2gfps(po_opens, PO);
	POpen_Stop_prep_opens(dest, PO, C, size);

	sent += dest.size() * size;
	rounds++;
}

/*
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
*/

void Processor::PTriple_Ext(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c)
{
	mpz_t ma, mb, mc;

	mpz_init(ma);
	mpz_init(mb);
	mpz_init(mc);

	if(0 != (*the_ext_lib.x_triple)(spdz_gfp_ext_handle, ma, mb, mc))
	{
		cerr << "Processor::PTriple_Ext extension library triple failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
	Pmpz2share(&ma, a);
	Pmpz2share(&mb, b);
	Pmpz2share(&mc, c);

	mpz_clear(ma);
	mpz_clear(mb);
	mpz_clear(mc);
}

void Processor::PInput_Ext(Share<gfp>& input_value, const int input_party_id)
{
	if(0 != (*the_ext_lib.x_input)(spdz_gfp_ext_handle, input_party_id, 1, &mpz_share_aux))
	{
		cerr << "Processor::PInput_Ext extension library input failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
	Pmpz2share(&mpz_share_aux, input_value);
}

/*
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
*/

/*
void Processor::PInput_Stop_Ext_64(int player, vector<int> targets)
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
*/

void Processor::PMult_Ext(const vector<int>& reg, int size)
{
	vector<int> sources, dest;
	int n = reg.size() / 3;
	sources.reserve(2 * n);
	dest.reserve(n);
	for (int i = 0; i < n; i++)
	{
		dest.push_back(reg[3 * i]);
		sources.push_back(reg[3 * i + 1]);
		sources.push_back(reg[3 * i + 2]);
	}

	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sources.size()*size);

	prep_shares(sources, Sh_PO, size);

	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sources.size()*size);

	//the share values are saved as mpz
	if(Sh_PO.size() > pm_size)
	{
		free_pm_mpz();
		alloc_pm_mpz(Sh_PO.size());
	}
	PShares2mpz(Sh_PO, pm_shares);

	if(0 != (*the_ext_lib.x_mult)(spdz_gfp_ext_handle, Sh_PO.size(), pm_shares, pm_products, 1))
	{
		cerr << "Processor::PMult_Start_Ext_64 extension library start_mult failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}

	PMult_Stop_prep_products(dest, size);

	sent += dest.size() * size;
	rounds++;
}

/*
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
	if(Sh_PO.size() > pm_size)
	{
		free_pm_mpz();
		alloc_pm_mpz(Sh_PO.size());
	}
	PShares2mpz(Sh_PO, pm_shares);

	if(0 != (*the_ext_lib.ext_start_mult)(spdz_gfp_ext_handle, Sh_PO.size(), pm_shares, pm_products, 1))
	{
		cerr << "Processor::PMult_Start_Ext_64 extension library start_mult failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}
*/

/*
void Processor::PMult_Stop_Ext_64(const vector<int>& reg, int size)
{
	if(0 != (*the_ext_lib.ext_stop_mult)(spdz_gfp_ext_handle))
	{
		cerr << "Processor::PMult_Stop_Ext_64 extension library stop_mult failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	PMult_Stop_prep_products(reg, size);

	sent += reg.size() * size;
	rounds++;
}
*/

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

void Processor::PAddm_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	to_bigint(*((bigint*)(&mpz_share_aux)), a.get_share());
	to_bigint(*((bigint*)(&mpz_arg_aux)), b);
	if(0 == (*the_ext_lib.x_mix_add)(spdz_gfp_ext_handle, mpz_share_aux, mpz_arg_aux))
	{
		Pmpz2share(&mpz_share_aux, c);
	}
	else
	{
		cerr << "Processor::PAddm_Ext extension library mix_add failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PSubml_Ext(Share<gfp>& a, gfp& b, Share<gfp>& c)
{
	to_bigint(*((bigint*)(&mpz_share_aux)), a.get_share());
	to_bigint(*((bigint*)(&mpz_arg_aux)), b);
	if(0 == (*the_ext_lib.x_mix_sub_scalar)(spdz_gfp_ext_handle, mpz_share_aux, mpz_arg_aux))
	{
		Pmpz2share(&mpz_share_aux, c);
	}
	else
	{
		cerr << "Processor::PSubml_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PSubmr_Ext(gfp& a, Share<gfp>& b, Share<gfp>& c)
{
	to_bigint(*((bigint*)(&mpz_share_aux)), b.get_share());
	to_bigint(*((bigint*)(&mpz_arg_aux)), a);
	if(0 == (*the_ext_lib.x_mix_sub_share)(spdz_gfp_ext_handle, mpz_arg_aux, mpz_share_aux))
	{
		Pmpz2share(&mpz_share_aux, c);
	}
	else
	{
		cerr << "Processor::PSubmr_Ext extension library mix_sub_scalar failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PLdsi_Ext(gfp& value, Share<gfp>& share)
{
	to_bigint(*((bigint*)(&mpz_arg_aux)), value);
	if(0 == (*the_ext_lib.x_share_immediates)(spdz_gfp_ext_handle, 0, 1, &mpz_arg_aux, &mpz_share_aux))
	{
		Pmpz2share(&mpz_share_aux, share);
	}
	else
	{
		cerr << "Processor::PLdsi_Ext extension library share_immediates failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PBit_Ext(Share<gfp>& share)
{
	if(0 == (*the_ext_lib.x_bit)(spdz_gfp_ext_handle, mpz_share_aux))
	{
		Pmpz2share(&mpz_share_aux, share);
	}
	else
	{
		cerr << "Processor::PBit_Ext extension library bit failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PInverse_Ext(Share<gfp>& share_value, Share<gfp>& share_inverse)
{
	if(0 == (*the_ext_lib.x_inverse)(spdz_gfp_ext_handle, mpz_share_aux, mpz_arg_aux))
	{
		Pmpz2share(&mpz_share_aux, share_value);
		Pmpz2share(&mpz_arg_aux, share_inverse);
	}
	else
	{
		cerr << "Processor::PInverse_Ext extension library inverse failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PMulm_Ext(Share<gfp>& sec_product, const Share<gfp>& sec_factor, const gfp & clr_factor)
{
	to_bigint(*((bigint*)(&mpz_share_aux)), sec_factor.get_share());
	to_bigint(*((bigint*)(&mpz_arg_aux)), clr_factor);
	if(0 == (*the_ext_lib.x_mix_mul)(spdz_gfp_ext_handle, mpz_share_aux, mpz_arg_aux))
	{
		Pmpz2share(&mpz_share_aux, sec_product);
	}
	else
	{
		cerr << "Processor::PMulm_Ext extension library mix_mul failed." << endl;
		dlclose(the_ext_lib.x_lib_handle);
		abort();
	}
}

void Processor::PShares2mpz(const vector< Share<gfp> >& shares, mpz_t * share_values)
{
	size_t count = shares.size();
	for(size_t i = 0; i < count; i++)
	{
		to_bigint(*((bigint*)(share_values + i)), shares[i].get_share());
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
	*(void**)(&x_share_immediates) = NULL;
	*(void**)(&x_bit) = NULL;
	*(void**)(&x_inverse) = NULL;


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
	LOAD_LIB_METHOD("share_immediates", x_share_immediates)
	LOAD_LIB_METHOD("bit", x_bit)
	LOAD_LIB_METHOD("inverse", x_inverse)

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
