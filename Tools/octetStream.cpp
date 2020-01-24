// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#include <fcntl.h>

#include "octetStream.h"

#include <string.h>
#include "Networking/sockets.h"
#include "Tools/sha1.h"
#include "Exceptions/Exceptions.h"
#include "Networking/data.h"
#include "Tools/time-func.h"
#include "Tools/int.h"

#include <openssl/sha.h>
#include <openssl/rand.h>


void octetStream::clear()
{
    if (data)
        delete[] data;
    data = 0;
    len = mxlen = ptr = 0;
}

void octetStream::assign(const octetStream& os)
{
  if (os.len>=mxlen)
    {
      if (data)
        delete[] data;
      mxlen=os.mxlen;  
      data=new octet[mxlen];
    }
  len=os.len;
  memcpy(data,os.data,len*sizeof(octet));
  ptr=os.ptr;
}


void octetStream::swap(octetStream& os)
{
  const size_t size = sizeof(octetStream);
  char tmp[size];
  memcpy(tmp, this, size);
  memcpy(this, &os, size);
  memcpy(&os, tmp, size);
}


octetStream::octetStream(size_t maxlen)
{
  mxlen=maxlen; len=0; ptr=0;
  data=new octet[mxlen];
}


octetStream::octetStream(const octetStream& os)
{
  mxlen=os.mxlen;
  len=os.len;
  data=new octet[mxlen];
  memcpy(data,os.data,len*sizeof(octet));
  ptr=os.ptr;
}


void octetStream::hash(octetStream& output) const
{
  SHA256(data, len, output.data);
  output.len=16;
}


octetStream octetStream::hash() const
{
  octetStream h(16);
  hash(h);
  return h;
}


bool octetStream::equals(const octetStream& a) const
{
  if (len!=a.len) { return false; }
  for (size_t i=0; i<len; i++)
    { if (data[i]!=a.data[i]) { return false; } }
  return true;
}


void octetStream::append_random(size_t num)
{
  resize(len+num);
  RAND_bytes(data+len, num);
  len+=num;
}


void octetStream::concat(const octetStream& os)
{
  resize(len+os.len);
  memcpy(data+len,os.data,os.len*sizeof(octet));
  len+=os.len;
}


void octetStream::store_bytes(octet* x, const size_t l)
{
  resize(len+4+l); 
  encode_length(data+len,l,4); len+=4;
  memcpy(data+len,x,l*sizeof(octet));
  len+=l;
}

void octetStream::get_bytes(octet* ans, size_t& length)
{
  length=decode_length(data+ptr,4); ptr+=4;
  memcpy(ans,data+ptr,length*sizeof(octet));
  ptr+=length;
}

void octetStream::store_int(size_t l, int n_bytes)
{
  resize(len+n_bytes);
  encode_length(data+len,l,n_bytes);
  len+=n_bytes;
}

void octetStream::store(int l)
{
  resize(len+4);
  encode_length(data+len,l,4);
  len+=4;
}


size_t octetStream::get_int(int n_bytes)
{
  size_t res=decode_length(data+ptr,n_bytes);
  ptr+=n_bytes;
  return res;
}

void octetStream::get(int& l)
{
  l=decode_length(data+ptr,4);
  ptr+=4;
}


void octetStream::exchange(int send_socket, int receive_socket, octetStream& receive_stream)
{
  send(send_socket, len, LENGTH_SIZE);
  const size_t buffer_size = 100000;
  size_t sent = 0, received = 0;
  bool length_received = false;
  size_t new_len = 0;
#ifdef TIME_ROUNDS
  Timer recv_timer;
#endif
  while (received < new_len or sent < len or not length_received)
    {
      if (sent < len)
        {
          size_t to_send = min(buffer_size, len - sent);
          send(send_socket, data + sent, to_send);
          sent += to_send;
        }

      // avoid extra branching, false before length received
      if (received < new_len)
        {
          // same buffer for sending and receiving
          // only receive up to already sent data
          // or when all is sent
          size_t to_receive = 0;
          if (sent == len)
            to_receive = new_len - received;
          else if (sent > received)
            to_receive = sent - received;
          if (to_receive > 0)
            {
#ifdef TIME_ROUNDS
              TimeScope ts(recv_timer);
#endif
              received += receive_non_blocking(receive_socket,
                  receive_stream.data + received, to_receive);
            }
        }
      else if (not length_received)
        {
#ifdef TIME_ROUNDS
          TimeScope ts(recv_timer);
#endif
          octet blen[LENGTH_SIZE];
          if (receive_all_or_nothing(receive_socket,blen,LENGTH_SIZE) == LENGTH_SIZE)
            {
              new_len=decode_length(blen,sizeof(blen));
              receive_stream.resize(max(new_len, len));
              length_received = true;
            }
        }
    }

#ifdef TIME_ROUNDS
  cout << "Exchange time: " << recv_timer.elapsed() << " seconds to receive "
      << 1e-3 * new_len << " KB" << endl;
#endif
  receive_stream.len = new_len;
  receive_stream.reset_read_head();
}


void octetStream::store(const vector<int>& v)
{
  store(v.size());
  for (int x : v)
    store(x);
}


void octetStream::get(vector<int>& v)
{
  size_t size;
  get(size);
  v.resize(size);
  for (int& x : v)
    get(x);
}

void octetStream::input(istream& s)
{
  size_t size;
  s.read((char*)&size, sizeof(size));
  resize_precise(size);
  s.read((char*)data, size);
}

void octetStream::output(ostream& s)
{
  s.write((char*)&len, sizeof(len));
  s.write((char*)data, len);
}

ostream& operator<<(ostream& s,const octetStream& o)
{
  for (size_t i=0; i<o.len; i++)
    { int t0=o.data[i]&15;
      int t1=o.data[i]>>4;
      s << hex << t1 << t0 << dec;
    }
  return s;
}




