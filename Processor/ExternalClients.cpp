// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

#include "Processor/ExternalClients.h"
#include <netinet/in.h>
#include <arpa/inet.h>
#include <thread>

ExternalClients::ExternalClients(int party_num, const string& prep_data_dir):
   party_num(party_num), prep_data_dir(prep_data_dir), server_connection_count(-1)
{
}

ExternalClients::~ExternalClients() 
{
  // close client sockets
  for (map<int,int>::iterator it = external_client_sockets.begin();
    it != external_client_sockets.end(); it++)
  {
    if (close(it->second))
    {
       error("failed to close external client connection socket)");
    }
  }
  for (map<int,AnonymousServerSocket*>::iterator it = client_connection_servers.begin();
    it != client_connection_servers.end(); it++)
  {
    delete it->second;
  }
  for (map<int,octet*>::iterator it = symmetric_client_keys.begin();
    it != symmetric_client_keys.end(); it++)
  {
    delete[] it->second;
  }
  for (map<int, pair<vector<octet>,uint64_t> >::iterator it_cs = symmetric_client_commsec_send_keys.begin();
    it_cs != symmetric_client_commsec_send_keys.end(); it_cs++)
  {
    memset(&(it_cs->second.first[0]), 0, it_cs->second.first.size());
  }
  for (map<int, pair<vector<octet>,uint64_t> >::iterator it_cs = symmetric_client_commsec_recv_keys.begin();
    it_cs != symmetric_client_commsec_recv_keys.end(); it_cs++)
  {
    memset(&(it_cs->second.first[0]), 0, it_cs->second.first.size());
  }
}

void ExternalClients::start_listening(int portnum_base)
{
  client_connection_servers[portnum_base] = new AnonymousServerSocket(portnum_base + get_party_num());
  client_connection_servers[portnum_base]->init();
  cerr << "Start listening on thread " << this_thread::get_id() << endl;
  cerr << "Party " << get_party_num() << " is listening on port " << (portnum_base + get_party_num()) 
        << " for external client connections." << endl;
}

int ExternalClients::get_client_connection(int portnum_base)
{
  map<int,AnonymousServerSocket*>::iterator it = client_connection_servers.find(portnum_base);
  if (it == client_connection_servers.end())
  {
    cerr << "Thread " << this_thread::get_id() << " didn't find server." << endl; 
    return -1;
  }
  cerr << "Thread " << this_thread::get_id() << " found server." << endl; 
  int client_id, socket;
  socket = client_connection_servers[portnum_base]->get_connection_socket(client_id);
  external_client_sockets[client_id] = socket;
  cerr << "Party " << get_party_num() << " received external client connection from client id: " << dec << client_id << endl;
  return client_id;
}

int ExternalClients::connect_to_server(int portnum_base, int ipv4_address)
{
  struct in_addr addr = { (unsigned int)ipv4_address };
  int csocket;
  const char* address_str = inet_ntoa(addr);
  cerr << "Party " << get_party_num() << " connecting to server at " << address_str << " on port " << portnum_base + get_party_num() << endl;
  set_up_client_socket(csocket, address_str, portnum_base + get_party_num());
  cerr << "Party " << get_party_num() << " connected to server at " << address_str << " on port " << portnum_base + get_party_num() << endl;
  int server_id = server_connection_count;
  // server identifiers are -1, -2, ... to avoid conflict with client identifiers
  server_connection_count--;
  external_client_sockets[server_id] = csocket;
  return server_id;
}


void ExternalClients::require_ed25519_keys()
{
    if (!ed25519_keys_loaded)
        throw "Ed25519 keys required but not found in player key files";
}

int ExternalClients::get_party_num() 
{
  return party_num;
}

