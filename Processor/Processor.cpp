// (C) 2017 University of Bristol. See License.txt


#include "Processor/Processor.h"
#include "Networking/STS.h"
#include "Auth/MAC_Check.h"

#include "Auth/fake-stuff.h"
#include <sodium.h>
#include <string>

#if defined(EXTENDED_SPDZ_32) || defined(EXTENDED_SPDZ_64)
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
{
  reset(program,0);

  public_input.open(get_filename("Programs/Public-Input/",false).c_str());
  private_input.open(get_filename("Player-Data/Private-Input-",true).c_str());
  public_output.open(get_filename("Player-Data/Public-Output-",true).c_str(), ios_base::out);
  private_output.open(get_filename("Player-Data/Private-Output-",true).c_str(), ios_base::out);

#if defined(EXTENDED_SPDZ_32) || defined(EXTENDED_SPDZ_64)
    spdz_ext_handle = NULL;
	cout << "SPDZ extension library initializing." << endl;
	stringstream ss;
	ss << gfp::pr();
	if(0 != (*the_ext_lib.ext_init)(&spdz_ext_handle, P.my_num(), P.num_players(), ss.str().c_str(), 10))
	{
		cerr << "SPDZ extension library initialization failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	cout << "SPDZ extension library initialized." << endl;
#endif
}

Processor::~Processor()
{
  cerr << "Sent " << sent << " elements in " << rounds << " rounds" << endl;
#if defined(EXTENDED_SPDZ_32) || defined(EXTENDED_SPDZ_64)
	(*the_ext_lib.ext_term)(spdz_ext_handle);
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
void Processor::POpen_Start_prep_shares(const vector<int>& reg, vector< Share<T> >& shares, int size)
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

	POpen_Start_prep_shares(reg, Sh_PO, size);

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

#if defined(EXTENDED_SPDZ_32) || defined(EXTENDED_SPDZ_64)

template <class T>
void uint2gfps(vector<gfp> & values, const T * uint_values, const size_t uint_value_count)
{
	values.resize(uint_value_count);
	for(size_t i = 0; i < uint_value_count; i++)
	{
		values[i].assign((u_int64_t)uint_values[i]);
	}
}

template <class T>
void gfp2uint(const gfp & gfp_value, T & t)
{
	bigint bi_value;
	to_bigint(bi_value, gfp_value);
	t = mpz_get_ui(bi_value.get_mpz_t());
}

template <class T>
void Processor::uint2share(const T in_value, Share<gfp> & out_value)
{
	gfp mac, value;
	value.assign((u_int64_t)in_value);
	mac.mul(MCp.get_alphai(), value);
	out_value.set_share(value);
	out_value.set_mac(mac);
}

#endif

#if defined(EXTENDED_SPDZ_32)

void Processor::POpen_Start_Ext_32(const vector<int>& reg, int size)
{
	int sz=reg.size();

	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	POpen_Start_prep_shares(reg, Sh_PO, size);

	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sz*size);

	//the share values are saved as unsigned long
	std::vector<u_int32_t> ui_share_values;
	shares2ui(Sh_PO, ui_share_values);
	if(Sh_PO.size() == ui_share_values.size())
	{
		//the extension library is given the shares' values and returns opens' values
		if(0 != (*the_ext_lib.ext_start_open)(spdz_ext_handle, ui_share_values.size(), &ui_share_values[0], 1))
		{
			cerr << "Processor::POpen_Start_Ext_32 extension library start_open failed." << endl;
			dlclose(the_ext_lib.ext_lib_handle);
			abort();
		}
		else
		{
			cout << "Processor::POpen_Start_Ext_32 extension start open launched." << endl;
		}
	}
	else
	{
		cout << "Processor::POpen_Start_Ext_32 ui_share_values size mismatch with PO_shares." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}

void Processor::POpen_Stop_Ext_32(const vector<int>& reg, int size)
{
	//vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	vector<gfp>& PO = get_PO<gfp>();
	vector<gfp>& C = get_C<gfp>();
	int sz=reg.size();
	PO.resize(sz*size);

	size_t open_count = 0;
	u_int32_t * opens = NULL;
	if(0 != (*the_ext_lib.ext_stop_open)(spdz_ext_handle, &open_count, &opens))
	{
		cerr << "Processor::POpen_Stop_Ext_32 extension library stop_open failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(NULL != opens)
	{
		if(PO.size() != open_count)
		{
			cerr << "Processor::POpen_Stop_Ext_32 size mismatch between share and open values array." << endl;
			dlclose(the_ext_lib.ext_lib_handle);
			abort();
		}
		uint2gfps(PO, opens, open_count);
		delete []opens;
		opens = NULL;
		open_count = 0;
	}
	else
	{
		cerr << "Processor::POpen_Stop_Ext_32 null open values array returned." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	POpen_Stop_prep_opens(reg, PO, C, size);

	sent += reg.size() * size;
	rounds++;
}

void Processor::Triple_Ext_32(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c)
{
	u_int32_t ui_a, ui_b, ui_c;
	if(0 != (*the_ext_lib.ext_triple)(spdz_ext_handle, &ui_a, &ui_b, &ui_c))
	{
		cerr << "Processor::Triple_Ext_32 extension library triple failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	uint2share(ui_a, a);
	uint2share(ui_b, b);
	uint2share(ui_c, c);
}

void Processor::Input_Ext_32(Share<gfp>& input_value, const int input_party_id)
{
	u_int32_t ui_input_value;
	if(0 != (*the_ext_lib.ext_input)(spdz_ext_handle, input_party_id, &ui_input_value))
	{
		cerr << "Processor::Input_Ext_32 extension library input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	uint2share(ui_input_value, input_value);
}

void Processor::Input_Start_Ext_32(int player, int n_inputs)
{
	if(0 != (*the_ext_lib.ext_start_input)(spdz_ext_handle, player, n_inputs))
	{
		cerr << "Processor::Input_Start_Ext_32 extension library start input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}

void Processor::Input_Stop_Ext_32(int /*player*/, vector<int> targets)
{
	size_t input_count = 0;
	u_int32_t * inputs = NULL;
	if(0 != (*the_ext_lib.ext_stop_input)(spdz_ext_handle, &input_count, &inputs))
	{
		cerr << "Processor::Input_Stop_Ext_32 extension library stop input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(NULL == inputs)
	{
		cerr << "Processor::Input_Stop_Ext_32 extension library stop input returned null ptr." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(targets.size() == input_count)
	{
		for(size_t i = 0; i < input_count; ++i)
		{
			Share<gfp>& share = get_S_ref<gfp>(targets[i]);
			uint2share(inputs[i], share);
		}
	}
	else
	{
		cerr << "Processor::Input_Stop_Ext_32 extension library stop input returned mismatched number of inputs " << targets.size() << "/" << input_count << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	delete []inputs;
}

/*
u_int32_t Processor::gfp2ui(const gfp & gfp_value)
{
	bigint bi_value;
	to_bigint(bi_value, gfp_value);
	return mpz_get_ui(bi_value.get_mpz_t());
}*/

void Processor::shares2ui(const vector< Share<gfp> > & shares, std::vector< u_int32_t > & ui_values)
{
	ui_values.clear();
	for(vector< Share<gfp> >::const_iterator i = shares.begin(); i != shares.end(); ++i)
	{
		u_int32_t v;
		gfp2uint(i->get_share(), v);
		ui_values.push_back(v);
	}
}

void Processor::test_extension_conversion(const gfp & original_gfp_value)
{
	u_int32_t outward_ui_value;
	gfp2uint(original_gfp_value, outward_ui_value);

	u_int32_t inward_ui_value = (*the_ext_lib.ext_test_conversion)(outward_ui_value);

	if(inward_ui_value != outward_ui_value)
	{
		cerr << "Processor::test_extension_conversion failed at unsigned long level " << inward_ui_value << " != " << outward_ui_value << endl;
		abort();
	}

	gfp restored_gfp_value;
	restored_gfp_value.assign((long)inward_ui_value);
	if(!original_gfp_value.equal(restored_gfp_value))
	{
		cerr << "Processor::test_extension_conversion failed at gfp level " << restored_gfp_value << " != " << original_gfp_value << endl;
		abort();
	}
}

#elif defined(EXTENDED_SPDZ_64)

void Processor::POpen_Start_Ext_64(const vector<int>& reg, int size)
{
	int sz=reg.size();

	vector< Share<gfp> >& Sh_PO = get_Sh_PO<gfp>();
	Sh_PO.clear();
	Sh_PO.reserve(sz*size);

	POpen_Start_prep_shares(reg, Sh_PO, size);

	vector<gfp>& PO = get_PO<gfp>();
	PO.resize(sz*size);

	//the share values are saved as unsigned long
	std::vector<u_int64_t> ul_share_values;
	shares2ul(Sh_PO, ul_share_values);
	if(Sh_PO.size() == ul_share_values.size())
	{
		//the extension library is given the shares' values and returns opens' values
		if(0 != (*the_ext_lib.ext_start_open)(spdz_ext_handle, ul_share_values.size(), &ul_share_values[0], 1))
		{
			cerr << "Processor::POpen_Start_Ext_64 extension library start_open failed." << endl;
			dlclose(the_ext_lib.ext_lib_handle);
			abort();
		}
		else
		{
			cout << "Processor::POpen_Start_Ext_64 extension start open launched." << endl;
		}
	}
	else
	{
		cout << "Processor::POpen_Start_Ext_64 ui_share_values size mismatch with PO_shares." << endl;
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

	size_t open_count = 0;
	u_int64_t * opens = NULL;
	if(0 != (*the_ext_lib.ext_stop_open)(spdz_ext_handle, &open_count, &opens))
	{
		cerr << "Processor::POpen_Stop_Ext_64 extension library stop_open failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(NULL != opens)
	{
		if(PO.size() != open_count)
		{
			cerr << "Processor::POpen_Stop_Ext_64 size mismatch between share and open values array." << endl;
			dlclose(the_ext_lib.ext_lib_handle);
			abort();
		}
		uint2gfps(PO, opens, open_count);
		delete []opens;
		opens = NULL;
		open_count = 0;
	}
	else
	{
		cerr << "Processor::POpen_Stop_Ext_64 null open values array returned." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	POpen_Stop_prep_opens(reg, PO, C, size);

	sent += reg.size() * size;
	rounds++;
}

void Processor::Triple_Ext_64(Share<gfp>& a, Share<gfp>& b, Share<gfp>& c)
{
	u_int64_t ul_a, ul_b, ul_c;
	if(0 != (*the_ext_lib.ext_triple)(spdz_ext_handle, &ul_a, &ul_b, &ul_c))
	{
		cerr << "Processor::Triple_Ext_64 extension library triple failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	uint2share(ul_a, a);
	uint2share(ul_b, b);
	uint2share(ul_c, c);
}

void Processor::Input_Ext_64(Share<gfp>& input_value, const int input_party_id)
{
	u_int64_t ul_input_value;
	if(0 != (*the_ext_lib.ext_input)(spdz_ext_handle, input_party_id, &ul_input_value))
	{
		cerr << "Processor::Input_Ext_64 extension library input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
	uint2share(ul_input_value, input_value);
}

void Processor::Input_Start_Ext_64(int player, int n_inputs)
{
	if(0 != (*the_ext_lib.ext_start_input)(spdz_ext_handle, player, n_inputs))
	{
		cerr << "Processor::Input_Start_Ext_64 extension library start input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}
}

void Processor::Input_Stop_Ext_64(int /*player*/, vector<int> targets)
{
	size_t input_count = 0;
	u_int64_t * inputs = NULL;
	if(0 != (*the_ext_lib.ext_stop_input)(spdz_ext_handle, &input_count, &inputs))
	{
		cerr << "Processor::Input_Stop_Ext_64 extension library stop input failed." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(NULL == inputs)
	{
		cerr << "Processor::Input_Stop_Ext_64 extension library stop input returned null ptr." << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	if(targets.size() == input_count)
	{
		for(size_t i = 0; i < input_count; ++i)
		{
			Share<gfp>& share = get_S_ref<gfp>(targets[i]);
			uint2share(inputs[i], share);
		}
	}
	else
	{
		cerr << "Processor::Input_Stop_Ext_64 extension library stop input returned mismatched number of inputs " << targets.size() << "/" << input_count << endl;
		dlclose(the_ext_lib.ext_lib_handle);
		abort();
	}

	delete []inputs;
}

void Processor::shares2ul(const vector< Share<gfp> > & shares, std::vector< u_int64_t > & ul_values)
{
	ul_values.clear();
	for(vector< Share<gfp> >::const_iterator i = shares.begin(); i != shares.end(); ++i)
	{
		u_int64_t v;
		gfp2uint(i->get_share(), v);
		ul_values.push_back(v);
	}
}

void Processor::test_extension_conversion(const gfp & original_gfp_value)
{
	u_int64_t outward_ul_value;
	gfp2uint(original_gfp_value, outward_ul_value);

	u_int64_t inward_ul_value = (*the_ext_lib.ext_test_conversion)(outward_ul_value);

	if(inward_ul_value != outward_ul_value)
	{
		cerr << "Processor::test_extension_conversion failed at unsigned long level " << inward_ul_value << " != " << outward_ul_value << endl;
		abort();
	}

	gfp restored_gfp_value;
	restored_gfp_value.assign(inward_ul_value);
	if(!original_gfp_value.equal(restored_gfp_value))
	{
		cerr << "Processor::test_extension_conversion failed at gfp level " << restored_gfp_value << " != " << original_gfp_value << endl;
		abort();
	}
}

#endif

#if defined(EXTENDED_SPDZ_32) || defined(EXTENDED_SPDZ_64)


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
	LOAD_LIB_METHOD("test_conversion",ext_test_conversion)
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

#endif //EXTENDED_SPDZ

template void Processor::POpen_Start(const vector<int>& reg,const Player& P,MAC_Check<gf2n>& MC,int size);
template void Processor::POpen_Start(const vector<int>& reg,const Player& P,MAC_Check<gfp>& MC,int size);
template void Processor::POpen_Stop(const vector<int>& reg,const Player& P,MAC_Check<gf2n>& MC,int size);
template void Processor::POpen_Stop(const vector<int>& reg,const Player& P,MAC_Check<gfp>& MC,int size);
template void Processor::read_socket_private<gfp>(int client_id, const vector<int>& registers, bool send_macs);
template void Processor::read_socket_vector<gfp>(int client_id, const vector<int>& registers);
template void Processor::read_shares_from_file<gfp>(int start_file_pos, int end_file_pos_register, const vector<int>& data_registers);
template void Processor::write_shares_to_file<gfp>(const vector<int>& data_registers);
