// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

// Client key file format:
//      X25519  Public Key
//      X25519  Secret Key
//      Ed25519 Public Key
//      Ed25519 Secret Key
//      Server 1 X25519 Public Key
//      Server 1 Ed25519 Public Key
//      ... 
//      Server N Public Key
//      Server N Ed25519 Public Key
//
// Player key file format:
//      X25519  Public Key
//      X25519  Secret Key
//      Ed25519 Public Key
//      Ed25519 Secret Key
//      Number of clients [64 bit little endian]
//      Client 1 X25519 Public Key
//      Client 1 Ed25519 Public Key
//      ...
//      Client N X25519 Public Key
//      Client N Ed25519 Public Key
//      Number of servers [64 bit little endian]
//      Server 1 X25519 Public Key
//      Server 1 Ed25519 Public Key
//      ...
//      Server N X25519 Public Key
//      Server N Ed25519 Public Key
#include "Tools/octetStream.h"
#include "Networking/Player.h"
#include "Math/gf2n.h"
#include "Config.h"
#include <vector>
#include <iomanip>

namespace Config {
    class ConfigError : public std::exception
    {
        std::string s;

        public:
        ConfigError(std::string ss) : s(ss) {}
        ~ConfigError() throw () {}
        const char* what() const throw() { return s.c_str(); }
    };


    void print_vector(const vector<octet> &vec)
    {
        cerr << hex;
        for(size_t i = 0; i < vec.size(); i ++ ) {
            cerr << setfill('0') << setw(2) << (int)vec[i];
        }
        cerr << dec << endl;
    }

    uint64_t getW64le(ifstream &infile)
    {
        uint8_t buf[8];
        uint64_t res=0;
        infile.read((char*)buf,sizeof buf);

        if (!infile.good())
            throw ConfigError("getW64le: could not read from config file");

        for(size_t i = 0; i < sizeof buf ; i ++ ) {
            res |= ((uint64_t)buf[i]) << i*8;
        }

        return res;
    }

    void putW64le(ofstream &outf, uint64_t nr)
    {
        char buf[8];
        for(int i=0;i<8;i++) {
            char byte = (uint8_t)(nr >> (i*8));
            buf[i] = (char)byte;
        }
        outf.write(buf,sizeof buf);
    }

    const string default_player_config_file_prefix = "Player-SPDZ-Keys-P";
    string player_config_file(int player_number)
    {
        stringstream filename;
        filename << default_player_config_file_prefix << player_number;
        return filename.str();
    }

}
