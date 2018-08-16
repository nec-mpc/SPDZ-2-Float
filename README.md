A fork of the University Of Bristol SPDZ-2 Repository, with changes to support extending the SPDZ-2 Framework to run additional protocols. Changes performed by Bar Ilan Cryptography Research Group, NEC  Security Research Labs, and Data61 (Marcel Keller). This code is used in the publication "Generalizing the SPDZ Compiler For Other Protocols" accpeted for ACM-CCS 2018. A link to the eprint will follow once available.

We would like to thank to the team behind the SPDZ-2 framework, which is an extensive effort and an excellent contribution to the MPC community. Special thanks to Marcel Keller (now at data61) for his insights and explanantions making this work possible.   

(C) 2017 University of Bristol. See License.txt
Software for the SPDZ and MASCOT secure multi-party computation protocols.
See `Programs/Source/` for some example MPC programs, and `tutorial.md` for
a basic tutorial.
See also https://www.cs.bris.ac.uk/Research/CryptographySecurity/SPDZ

## SPDZ-2 With Extensions
The SPDZ-2 extensions is a mechanism that enables substitution of the original implementation of various operations with an alternate external implementation. This is done by dynamically loading a configured library and prescribed API function pointers. 
In runtime, the SPDZ-2 processor will call the loaded API functions instead of the original implementation and provide it with the required parameters.

### Detecting code changes
Code changes in the fork can be detected by searching for ```EXTENDED_SPDZ``` compiler directive. 
If this directive is swiched off, the original SPDZ-2 code will be compiled.
Here is an example of a code change, in [Instruction.cpp](https://github.com/cryptobiu/SPDZ-2/edit/master/Processor/Instruction.cpp)
```
      case STARTOPEN:
#if defined(EXTENDED_SPDZ)
    	  Proc.POpen_Start_Ext_64(start, size);
#else
    	  Proc.POpen_Start(start,Proc.P,Proc.MCp,size);
#endif
```

if the directive is defined, the code will call the SPDZ extension code instrad of the standard processing for a open.
### Extension API

The API for an extension is defined [in the following include file](https://github.com/cryptobiu/SPDZ-2-Extension-MpcHonestMajority/blob/master/spdzext.h)

```
int init(void ** handle, const int pid, const int num_of_parties, const int thread_id,
			const char * field, const int open_count, const int mult_count, const int bits_count)

int term(void * handle);

int offline(void * handle, const int offline_size);

int start_open(void * handle, const size_t share_count, const mpz_t * shares, mpz_t * opens, int verify);

int stop_open(void * handle);

int triple(void * handle, mpz_t * a, mpz_t * b, mpz_t * c);

int input(void * handle, const int input_of_pid, mpz_t * input_value);

int start_verify(void * handle, int * error);

int stop_verify(void * handle);

int start_input(void * handle, const int input_of_pid, const size_t num_of_inputs, mpz_t * inputs);

int stop_input(void * handle);

int start_mult(void * handle, const size_t share_count, const mpz_t * shares, mpz_t * products, int verify);

int stop_mult(void * handle);

int mix_add(void * handle, mpz_t * share, const mpz_t * scalar);

int mix_sub_scalar(void * handle, mpz_t * share, const mpz_t * scalar);

int mix_sub_share(void * handle, const mpz_t * scalar, mpz_t * share);

int start_share_immediates(void * handle, const size_t value_count, const mpz_t * values, mpz_t * shares);

int stop_share_immediates(void * handle);
 
int share_immediate(void * handle, const mpz_t * value, mpz_t * share);

int bit(void * handle, mpz_t * share);

int inverse(void * handle, mpz_t * share_value, mpz_t * share_inverse);

```
### Example of SPDZ-2 extension
See https://github.com/cryptobiu/SPDZ-2-Extension-MpcHonestMajority for an example of such implemented extension library.

## SPDZ-2 

#### Requirements:
 - GCC
 - MPIR library, compiled with C++ support (use flag --enable-cxx when running configure)
 - libsodium library, tested against 1.0.11
 - CPU supporting AES-NI and PCLMUL
 - Python 2.x, ideally with `gmpy` package (for testing)

#### OS X:
 - `g++` might actually refer to clang, in which case you need to change `CONFIG` to use GCC instead.
 - It has been reported that MPIR has to be compiled with GCC for the linking to work:
   ```./configure CC=<path to GCC gcc> CXX=<path to GCC g++> --enable-cxx```

#### To compile SPDZ:

1) Optionally, edit CONFIG and CONFIG.mine so that the following variables point to the right locations:
 - PREP_DIR: this should be a local, unversioned directory to store preprocessing data (defaults to Player-Data in the working directory)

2) Run make (use the flag -j for faster compilation with multiple threads)


#### To setup for the online phase

Run:

`Scripts/setup-online.sh`

This sets up parameters for the online phase for 2 parties with a 128-bit prime field and 40-bit binary field, and creates fake offline data (multiplication triples etc.) for these parameters.

Parameters can be customised by running

`Scripts/setup-online.sh <nparties> <nbitsp> <nbits2>`


#### To compile a program

To compile the program in `./Programs/Source/tutorial.mpc`, run:

`./compile.py tutorial`

This creates the bytecode and schedule files in Programs/Bytecode/ and Programs/Schedules/

#### To run a program

To run the above program (on one machine), first run:

`./Server.x 2 5000 &`

(or replace `5000` with your desired port number)

Then run both parties' online phase:

`./Player-Online.x -pn 5000 0 tutorial`

`./Player-Online.x -pn 5000 1 tutorial` (in a separate terminal)

Or, you can use a script to do the above automatically:

`Scripts/run-online.sh tutorial`

To run a program on two different machines, firstly the preprocessing data must be
copied across to the second machine (or shared using sshfs), and secondly, Player-Online.x
needs to be passed the machine where Server.x is running.
e.g. if this machine is name `diffie` on the local network:

`./Player-Online.x -pn 5000 -h diffie 0 tutorial`

`./Player-Online.x -pn 5000 -h diffie 1 tutorial`

#### Compiling and running programs from external directories

Programs can also be edited, compiled and run from any directory with the above basic structure. So for a source file in `./Programs/Source/`, all SPDZ scripts must be run from `./`. The `setup-online.sh` script must also be run from `./` to create the relevant data. For example:

```
spdz$ cd ../
$ mkdir myprogs
$ cd myprogs
$ mkdir -p Programs/Source
$ vi Programs/Source/test.mpc
$ ../spdz/compile.py test.mpc
$ ls Programs/
Bytecode  Public-Input  Schedules  Source
$ ../spdz/Scripts/setup-online.sh
$ ls
Player-Data Programs
$ ../spdz/Scripts/run-online.sh test
```

#### Offline phase (MASCOT)

In order to compile the MASCOT code, the following must be set in CONFIG or CONFIG.mine:

`USE_GF2N_LONG = 1`

It also requires SimpleOT:
```
git submodule update --init SimpleOT
cd SimpleOT
make
```

If SPDZ has been built before, any compiled code needs to be removed:

`make clean`

HOSTS must contain the hostnames or IPs of the players, see HOSTS.example for an example.

Then, MASCOT can be run as follows:

`host1:$ ./ot-offline.x -p 0 -c`

`host2:$ ./ot-offline.x -p 1 -c`
