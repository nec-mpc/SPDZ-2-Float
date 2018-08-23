# SPDZ-2 With Extensions

A fork of the University Of Bristol SPDZ-2 Repository, with changes to support extending the SPDZ-2 Framework to run additional protocols. Changes performed by Bar Ilan Cryptography Research Group and NEC Security Research Labs. This code is used in the publication "Generalizing the SPDZ Compiler For Other Protocols" accpeted for ACM-CCS 2018. A link to the eprint will follow once available.

We would like to thank to the team behind the SPDZ-2 framework, which is an extensive effort and an excellent contribution to the MPC community. Special thanks to Marcel Keller for his numerous insights and explanantions making this work possible.   

(C) 2017 University of Bristol. See License.txt
Software for the SPDZ and MASCOT secure multi-party computation protocols.
See `Programs/Source/` for some example MPC programs, and `tutorial.md` for
a basic tutorial.
See also https://www.cs.bris.ac.uk/Research/CryptographySecurity/SPDZ

## SPDZ-2 With Extensions - rationale
The SPDZ-2 extensions is a mechanism that enables substitution of the original implementation of various operations with an alternate external implementation. This is done by dynamically loading a configured library and prescribed API function pointers. 
In runtime, the SPDZ-2 processor will call the loaded API functions instead of the original implementation and provide it with the required parameters.

### MPC programs source code
The [Programs/Source](https://github.com/cryptobiu/SPDZ-2/tree/master/Programs/Source) folder of this fork contains MPC programs added as paro of our work to evaluate different protocols under the framework. For example, the following program evaluates a decision tree.  
```
import util
#------------------------------------------------------------------------------
#definitions
c_FeaturesSetSize = 17
c_TreeDepth = 30
c_NodeSetSize = 1255

#user 0 the evaluator
#user 1 is the evaluee
#------------------------------------------------------------------------------
# Code for oblivious selection of an array member by a secure index
def oblivious_selection(sec_array, array_size, sec_index):
    bitcnt = util.log2(array_size)
    sec_index_bits = sec_index.bit_decompose(bitcnt)
    return obliviously_select(sec_array, array_size, 0, sec_index_bits, len(sec_index_bits) - 1)

def obliviously_select(array, size, offset, bits, bits_index):
    #print('size={}; offset={}; bi={};'.format(size, offset, bits_index))
    if offset >= size:
        return 0
    elif bits_index < 0:
        return array[offset]
    else:
        half_size = 2**(bits_index)
        msb = bits[bits_index]
        return msb.if_else(
            obliviously_select(array, size, offset + half_size, bits, bits_index-1) , 
            obliviously_select(array, size, offset, bits, bits_index-1) )

#------------------------------------------------------------------------------
# Reading feature set from user 1 (the evaluee)
#print_ln('user 1: please enter input offset:')
User1InputOffset = sint.get_input_from(1)
#print_ln('user 1: please enter feature set (%s feature values):', c_FeaturesSetSize)
FeaturesSet = Array(c_FeaturesSetSize, sint)
@for_range(c_FeaturesSetSize)
def init_features_set(i):
    FeaturesSet[i] = sint.get_input_from(1) - User1InputOffset
    #debug-print
    #print_ln('FeaturesSet[%s] = %s', i, FeaturesSet[i].reveal())

#------------------------------------------------------------------------------
def test(FeatureIdx, Operator, Threshold):
    feature_value = oblivious_selection(FeaturesSet, c_FeaturesSetSize, FeatureIdx)
    return Operator.if_else(feature_value > Threshold, feature_value == Threshold)
    
#------------------------------------------------------------------------------
#print_ln('user 0: please enter input offset:')
User0InputOffset = sint.get_input_from(0)
def read_node(i):
    #print_ln('user 0: please enter node %s feature index:', i)
    FeatureIdx = sint.get_input_from(0) - User0InputOffset
    #debug-print
    #print_ln('FeatureIdx[%s] = %s', i, FeatureIdx.reveal())

    #print_ln('user 0: please enter node %s operator:', i)
    Operator = sint.get_input_from(0) - User0InputOffset
    #debug-print
    #print_ln('Operator[%s] = %s', i, Operator.reveal())

    #print_ln('user 0: please enter node %s Threshold:', i)
    Threshold = sint.get_input_from(0) - User0InputOffset
    #debug-print
    #print_ln('Threshold[%s] = %s', i, Threshold.reveal())
    
    #print_ln('user 0: please enter node %s GT/EQ:', i)
    GT_or_EQ = sint.get_input_from(0) - User0InputOffset
    #debug-print
    #print_ln('GT_or_EQ[%s] = %s', i, GT_or_EQ.reveal())

    #print_ln('user 0: please enter node %s LTE/NEQ:', i)
    LTE_or_NEQ = sint.get_input_from(0) - User0InputOffset
    #debug-print
    #print_ln('LTE_or_NEQ[%s] = %s', i, LTE_or_NEQ.reveal())

    NodePass = test(FeatureIdx, Operator, Threshold)

    #debug-print
    #print_ln('Node[%s] passage = %s', i, NodePass.reveal())
    
    return NodePass*GT_or_EQ + (1 - NodePass)*LTE_or_NEQ
#------------------------------------------------------------------------------
# Reading node set from user 0 (the evaluator)
NodeSet = Array(c_NodeSetSize, sint)
@for_range(c_NodeSetSize)
def read_node_loop(i):
    NodeSet[i] = read_node(i)
#------------------------------------------------------------------------------
#evaluation
NodePtr = MemValue(sint(0))
@for_range(c_TreeDepth)
def evaluation_loop(c_CurrLyr):
    NextNodePtr = oblivious_selection(NodeSet, c_NodeSetSize, NodePtr)
    CycleBack = (NextNodePtr < 0) * (c_CurrLyr < (c_TreeDepth-1))
    NodePtr.write(CycleBack.if_else(NodePtr, NextNodePtr))
    #debug-print
    #print_ln('CurrentLayer = %s; NodePtr = %s; NextNodePtr = %s', c_CurrLyr, NodePtr.reveal(), NextNodePtr.reveal())

NodePtr = (NodePtr + 1) * (-1)
print_ln('evaluation result = %s', NodePtr.reveal())
```
### marking of code changes
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
	 const char * field, const int open_count, const int mult_count, const int bits_count);

int term(void * handle);

int offline(void * handle, const int offline_size);

int opens(void * handle, const size_t share_count, const mpz_t * shares, mpz_t * opens, int verify);

int triple(void * handle, mpz_t a, mpz_t b, mpz_t c);

int verify(void * handle, int * error);

int input(void * handle, const int input_of_pid, const size_t num_of_inputs, mpz_t * inputs);

int mult(void * handle, const size_t share_count, const mpz_t * shares, mpz_t * products, int verify);

int mix_add(void * handle, mpz_t share, const mpz_t scalar);

int mix_sub_scalar(void * handle, mpz_t share, const mpz_t scalar);

int mix_sub_share(void * handle, const mpz_t scalar, mpz_t share);

int share_immediates(void * handle, const int party_id, const size_t value_count, const mpz_t * values, mpz_t * shares);

int bit(void * handle, mpz_t share);

int inverse(void * handle, mpz_t share_value, mpz_t share_inverse);

```
### Example of SPDZ-2 extension
See https://github.com/cryptobiu/SPDZ-2-Extension-MpcHonestMajority for an example of such implemented extension library.


## SPDZ-2 

#### Requirements:
 - GCC (tested with 7.2) or LLVM (tested with 3.8)
 - MPIR library, compiled with C++ support (use flag --enable-cxx when running configure)
 - libsodium library, tested against 1.0.11
 - CPU supporting AES-NI and PCLMUL
 - Python 2.x
 - If using macOS, Sierra or later

#### To compile:

1) Edit `CONFIG` or `CONFIG.mine`:

 - Add the following line at the top: `MY_CFLAGS = -DINSECURE`
 - For processors without AVX (e.g., Intel Atom) or for optimization, set `ARCH = -march=<architecture>`.

2) Run `make bmr` (use the flag -j for faster compilation multiple threads). Remember to run `make clean` first after changing `CONFIG` or `CONFIG.mine`.

#### Configure the parameters:

1) Edit `Program/Source/gc_oram.mpc` to change size and to choose Circuit ORAM or linear scan without ORAM.
2) Run `./compile.py -D gc_oram`.

#### Run the protocol:

- Run everything locally: `Scripts/bmr-program-run.sh gc_oram`.
- Run on different hosts: `Scripts/bmr-program-run-remote.sh gc_oram <host1> <host2> [...]`

To run with more than two parties, change `CFLAGS = -DN_PARTIES=<n>` in `CONFIG`, and compile again after `make clean`.
