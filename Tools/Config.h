// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

#include "Tools/octetStream.h"
#include "Networking/Player.h"
namespace Config {
    typedef vector<octet> public_key;
    typedef vector<octet> public_signing_key;
    typedef vector<octet> secret_key;
    typedef vector<octet> secret_signing_key;

    uint64_t getW64le(ifstream &infile);
    void putW64le(ofstream &outf, uint64_t nr);
    extern const string default_player_config_file_prefix;
    string player_config_file(int player_number);
    void print_vector(const vector<octet> &vec);
}
