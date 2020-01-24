# Confidential:
# (C) 2017 University of Bristol. See License.txt

import itertools
from random import randint
import time
import inspect
import functools
from Compiler.exceptions import *
from Compiler.config import *
from Compiler import util
import copy


###
### Opcode constants
###
### Whenever these are changed the corresponding enums in Processor/instruction.h
### MUST also be changed. (+ the documentation)
###
opcodes = dict(
    # Load/store
    LDI = 0x1,
    LDSI = 0x2,
    LDMC = 0x3,
    LDMS = 0x4,
    STMC = 0x5,
    STMS = 0x6,
    LDMCI = 0x7,
    LDMSI = 0x8,
    STMCI = 0x9,
    STMSI = 0xA,
    MOVC = 0xB,
    MOVS = 0xC,
    PROTECTMEMS = 0xD,
    PROTECTMEMC = 0xE,
    PROTECTMEMINT = 0xF,
    LDMINT = 0xCA,
    STMINT = 0xCB,
    LDMINTI = 0xCC,
    STMINTI = 0xCD,
    PUSHINT = 0xCE,
    POPINT = 0xCF,
    MOVINT = 0xD0,
    # Machine
    LDTN = 0x10,
    LDARG = 0x11,
    REQBL = 0x12,
    STARG = 0x13,
    TIME = 0x14,
    START = 0x15,
    STOP = 0x16,
    USE = 0x17,
    USE_INP = 0x18,
    RUN_TAPE = 0x19,
    JOIN_TAPE = 0x1A,
    CRASH = 0x1B,
    USE_PREP = 0x1C,
    # Addition
    ADDC = 0x20,
    ADDS = 0x21,
    ADDM = 0x22,
    ADDCI = 0x23,
    ADDSI = 0x24,
    SUBC = 0x25,
    SUBS = 0x26,
    SUBML = 0x27,
    SUBMR = 0x28,
    SUBCI = 0x29,
    SUBSI = 0x2A,
    SUBCFI = 0x2B,
    SUBSFI = 0x2C,
    # Multiplication/division
    MULC = 0x30,
    MULM = 0x31,
    MULCI = 0x32,
    MULSI = 0x33,
    DIVC = 0x34,
    DIVCI = 0x35,
    MODC = 0x36,
    MODCI = 0x37,
    LEGENDREC = 0x38,
    DIGESTC = 0x39,
    E_STARTMULT = 0x40,
    E_STOPMULT = 0x41,
    E_MULTI_STARTMULT = 0x42,
    E_MULTI_STOPMULT = 0x43,
    GMULBITC = 0x136,
    GMULBITM = 0x137,
    # Open
    STARTOPEN = 0xA0,
    STOPOPEN = 0xA1,
    E_STARTOPEN = 0xA2,
    E_STOPOPEN = 0xA3,
    E_MULT = 0xA4,
    OPEN = 0xA5,
    # Data access
    TRIPLE = 0x50,
    BIT = 0x51,
    SQUARE = 0x52,
    INV = 0x53,
    GBITTRIPLE = 0x154,
    GBITGF2NTRIPLE = 0x155,
    INPUTMASK = 0x56,
    PREP = 0x57,
    # Input
    INPUT = 0x60,
    STARTINPUT = 0x61,
    STOPINPUT = 0x62,  
    READSOCKETC = 0x63,
    READSOCKETS = 0x64,
    WRITESOCKETC = 0x65,
    WRITESOCKETS = 0x66,
    READSOCKETINT = 0x69,
    WRITESOCKETINT = 0x6a,
    WRITESOCKETSHARE = 0x6b,
    LISTEN = 0x6c,
    ACCEPTCLIENTCONNECTION = 0x6d,
    CONNECTIPV4 = 0x6e,
    READCLIENTPUBLICKEY = 0x6f,
    # Bitwise logic
    ANDC = 0x70,
    XORC = 0x71,
    ORC = 0x72,
    ANDCI = 0x73,
    XORCI = 0x74,
    ORCI = 0x75,
    NOTC = 0x76,
    # Bitwise shifts
    SHLC = 0x80,
    SHRC = 0x81,
    SHLCI = 0x82,
    SHRCI = 0x83,
    # Branching and comparison
    JMP = 0x90,
    JMPNZ = 0x91,
    JMPEQZ = 0x92,
    EQZC = 0x93,
    LTZC = 0x94,
    LTC = 0x95,
    GTC = 0x96,
    EQC = 0x97,
    JMPI = 0x98,
    # Integers
    LDINT = 0x9A,
    ADDINT = 0x9B,
    SUBINT = 0x9C,
    MULINT = 0x9D,
    DIVINT = 0x9E,
    PRINTINT = 0x9F,
    # Conversion
    CONVINT = 0xC0,
    CONVMODP = 0xC1,
    GCONVGF2N = 0x1C1,
    # IO
    PRINTMEM = 0xB0,
    PRINTREG = 0XB1,
    RAND = 0xB2,
    PRINTREGPLAIN = 0xB3,
    PRINTCHR = 0xB4,
    PRINTSTR = 0xB5,
    PUBINPUT = 0xB6,
    RAWOUTPUT = 0xB7,
    STARTPRIVATEOUTPUT = 0xB8,
    STOPPRIVATEOUTPUT = 0xB9,
    PRINTCHRINT = 0xBA,
    PRINTSTRINT = 0xBB,
    PRINTFLOATPLAIN = 0xBC,
    E_PRINTFIXEDPLAIN = 0x1BC,
    WRITEFILESHARE = 0xBD,    
    READFILESHARE = 0xBE,
    E_READ_FROM_FILE = 0xBF,
    GE_READ_FROM_FILE = 0x1BF,
    GBITDEC = 0x184,
    GBITCOM = 0x185,
    # bit de-compostion (BIU-NEC)
    E_BITDEC = 0x186,
    E_SKEW_DEC = 0x187,
    # bit injection (BIU-NEC)
    E_BITINJ = 0x188,
    E_SKEW_INJ = 0x189,
    # bit re-composition (BIU-NEC)
    E_BITREC = 0x190,
    E_SKEW_REC = 0x191,
    E_POST_REC = 0x192,
    # Secure socket
    INITSECURESOCKET = 0x1BA,
    RESPSECURESOCKET = 0x1BB,
    # added-extensions
    E_SKEW_BIT_DEC = 0x1D0,
    E_SKEW_BIT_REC = 0x1D1,
    E_SKEW_RING_REC = 0x1D2,
    E_SKEW_BIT_INJ = 0x1D3,
    E_INPUT_SHARE_INT = 0x203,
    GE_INPUT_SHARE_INT = 0x303,
    E_INPUT_SHARE_FIX = 0x204,
    E_INPUT_CLEAR_INT = 0x205,
    E_INPUT_CLEAR_FIX = 0x206,
    E_VERIFY_OPTIONAL_SUGGEST = 0x207,
    E_VERIFY_FINAL = 0x208,
    #E_START_MULT = 0x209,
    #E_STOP_MULT = 0x20A,
    E_START_OPEN = 0x20B,
    E_STOP_OPEN = 0x20C,
    E_MP_LDSI = 0x20D,
    E_MP_LDI = 0x20E,
    E_MP_ADDS = 0x300,
    E_MP_ADDM = 0x301,
    E_MP_ADDSI = 0x302,
    E_MP_SUBS=0x303,
    E_MP_SUBML = 0x304,
    E_MP_SUBMR = 0x305,
    E_MP_SUBSIL=0x306,
    E_MP_SUBSIR=0x307,
    E_MP_MULM=0x308,
    E_MP_MULSI=0x309,
    E_MP_OPEN=0x30A,
    E_MP_MULT=0x30B,
    E_MP_SKEW_BIT_DEC = 0x30C,
    E_MP_SKEW_RING_REC = 0x30D,
    E_MP_SKEW_BIT_INJ = 0x30E
)


def int_to_bytes(x):
    """ 32 bit int to big-endian 4 byte conversion. """
    return [(x >> 8*i) % 256 for i in (3,2,1,0)]


global_vector_size = 1
global_vector_size_depth = 0
global_instruction_type_stack = ['modp']

def set_global_vector_size(size):
    global global_vector_size, global_vector_size_depth
    if size == 1:
        return
    if global_vector_size == 1 or global_vector_size == size:
        global_vector_size = size
        global_vector_size_depth += 1
    else:
        raise CompilerError('Cannot set global vector size when already set')

def set_global_instruction_type(t):
    if t == 'modp' or t == 'gf2n':
        global_instruction_type_stack.append(t)
    else:
        raise CompilerError('Invalid type %s for setting global instruction type')

def reset_global_vector_size():
    global global_vector_size, global_vector_size_depth
    if global_vector_size_depth > 0:
        global_vector_size_depth -= 1
        if global_vector_size_depth == 0:
            global_vector_size = 1

def reset_global_instruction_type():
    global_instruction_type_stack.pop()

def get_global_vector_size():
    return global_vector_size

def get_global_instruction_type():
    return global_instruction_type_stack[-1]


def vectorize(instruction, global_dict=None):
    """ Decorator to vectorize instructions. """

    if global_dict is None:
        global_dict = inspect.getmodule(instruction).__dict__

    class Vectorized_Instruction(instruction):
        __slots__ = ['size']
        def __init__(self, size, *args, **kwargs):
            self.size = size
            super(Vectorized_Instruction, self).__init__(*args, **kwargs)
            for arg,f in zip(self.args, self.arg_format):
                if issubclass(ArgFormats[f], RegisterArgFormat):
                    arg.set_size(size)
        def get_code(self):
            return (self.size << 9) + self.code
        def get_pre_arg(self):
            return "%d, " % self.size
        def is_vec(self):
            return self.size > 1
        def get_size(self):
            return self.size
        def expand(self):
            set_global_vector_size(self.size)
            super(Vectorized_Instruction, self).expand()
            reset_global_vector_size()

    @functools.wraps(instruction)
    def maybe_vectorized_instruction(*args, **kwargs):
        if global_vector_size == 1:
            return instruction(*args, **kwargs)
        else:
            return Vectorized_Instruction(global_vector_size, *args, **kwargs)
    maybe_vectorized_instruction.vec_ins = Vectorized_Instruction
    maybe_vectorized_instruction.std_ins = instruction
    
    vectorized_name = 'v' + instruction.__name__
    Vectorized_Instruction.__name__ = vectorized_name
    global_dict[vectorized_name] = Vectorized_Instruction
    global_dict[instruction.__name__ + '_class'] = instruction
    return maybe_vectorized_instruction


def gf2n(instruction):
    """ Decorator to create GF_2^n instruction corresponding to a given
        modp instruction.

        Adds the new GF_2^n instruction to the globals dictionary. Also adds a
        vectorized GF_2^n instruction if a modp version exists. """
    global_dict = inspect.getmodule(instruction).__dict__

    if global_dict.has_key('v' + instruction.__name__):
        vectorized = True
    else:
        vectorized = False

    if isinstance(instruction, type) and issubclass(instruction, Instruction):
        instruction_cls = instruction
    else:
        try:
            instruction_cls = global_dict[instruction.__name__ + '_class']
        except KeyError:
            raise CompilerError('Cannot decorate instruction %s' % instruction)

    ### added to merge NEC and BIU (start)
    def reformat(arg_format):
        if isinstance(arg_format, list):
            __format = []
            for __f in arg_format:
                if __f in ('int', 'p', 'ci', 'str'):
                    __format.append(__f)
                else:
                    __format.append(__f[0] + 'g' + __f[1:])
            arg_format[:] = __format
        else:
            for __f in arg_format.args:
                reformat(__f)
    ### added to merge NEC and BIU (end)

    class GF2N_Instruction(instruction_cls):
        __doc__ = instruction_cls.__doc__.replace('c_', 'c^g_').replace('s_', 's^g_')
        __slots__ = []
        field_type = 'gf2n'
        if isinstance(instruction_cls.code, int):
            code = (1 << 8) + instruction_cls.code

        # set modp registers in arg_format to GF2N registers
        if 'gf2n_arg_format' in instruction_cls.__dict__:
            arg_format = instruction_cls.gf2n_arg_format
        elif isinstance(instruction_cls.arg_format, itertools.repeat):
            __f = instruction_cls.arg_format.next()
            if __f != 'int' and __f != 'p':
                arg_format = itertools.repeat(__f[0] + 'g' + __f[1:])
        else:
            arg_format = copy.deepcopy(instruction_cls.arg_format)
            reformat(arg_format)
            # __format = []
            # for __f in instruction_cls.arg_format:
            #     if __f in ('int', 'p', 'ci', 'str'):
            #         __format.append(__f)
            #     else:
            #         __format.append(__f[0] + 'g' + __f[1:])
            # arg_format = __format

        def is_gf2n(self):
            return True

        def expand(self):
            set_global_instruction_type('gf2n')
            super(GF2N_Instruction, self).expand()
            reset_global_instruction_type()

    GF2N_Instruction.__name__ = 'g' + instruction_cls.__name__
    if vectorized:
        vec_GF2N = vectorize(GF2N_Instruction, global_dict)

    @functools.wraps(instruction)
    def maybe_gf2n_instruction(*args, **kwargs):
        if get_global_instruction_type() == 'gf2n':
            if vectorized:
                return vec_GF2N(*args, **kwargs)
            else:
                return GF2N_Instruction(*args, **kwargs)
        else:
            return instruction(*args, **kwargs)
    
    # If instruction is vectorized, new GF2N instruction must also be
    if vectorized:
        global_dict[GF2N_Instruction.__name__] = vec_GF2N
    else:
        global_dict[GF2N_Instruction.__name__] = GF2N_Instruction

    global_dict[instruction.__name__ + '_class'] = instruction_cls
    return maybe_gf2n_instruction
    #return instruction


class RegType(object):
    """ enum-like static class for Register types """
    ClearModp = 'c'
    SecretModp = 's'
    ClearGF2N = 'cg'
    SecretGF2N = 'sg'
    ClearInt = 'ci'

    Types = [ClearModp, SecretModp, ClearGF2N, SecretGF2N, ClearInt]

    @staticmethod
    def create_dict(init_value_fn):
        """ Create a dictionary with all the RegTypes as keys """
        res = defaultdict(init_value_fn)
        # initialization for legacy
        for t in RegType.Types:
            res[t]
        return res

class ArgFormat(object):
    @classmethod
    def check(cls, arg):
        return NotImplemented

    @classmethod
    def encode(cls, arg):
        return NotImplemented

class RegisterArgFormat(ArgFormat):
    @classmethod
    def check(cls, arg):
        if not isinstance(arg, program.curr_tape.Register):
            raise ArgumentError(arg, 'Invalid register argument')
        if arg.i > REG_MAX:
            raise ArgumentError(arg, 'Register index too large')
        if arg.program != program.curr_tape:
            raise ArgumentError(arg, 'Register from other tape, trace: %s' % \
                                    util.format_trace(arg.caller))
        if arg.reg_type != cls.reg_type:
            raise ArgumentError(arg, "Wrong register type '%s', expected '%s'" % \
                                    (arg.reg_type, cls.reg_type))

    @classmethod
    def encode(cls, arg):
        return int_to_bytes(arg.i)

class ClearModpAF(RegisterArgFormat):
    reg_type = RegType.ClearModp

class SecretModpAF(RegisterArgFormat):
    reg_type = RegType.SecretModp

class ClearGF2NAF(RegisterArgFormat):
    reg_type = RegType.ClearGF2N

class SecretGF2NAF(RegisterArgFormat):
    reg_type = RegType.SecretGF2N

class ClearIntAF(RegisterArgFormat):
    reg_type = RegType.ClearInt

class IntArgFormat(ArgFormat):
    @classmethod
    def check(cls, arg):
        if not isinstance(arg, (int, long)):
            raise ArgumentError(arg, 'Expected an integer-valued argument')

    @classmethod
    def encode(cls, arg):
        return int_to_bytes(arg)

class ImmediateModpAF(IntArgFormat):
    @classmethod
    def check(cls, arg):
        super(ImmediateModpAF, cls).check(arg)
        if arg >= 2**31 or arg < -2**31:
            raise ArgumentError(arg, 'Immediate value outside of 32-bit range')

class ImmediateGF2NAF(IntArgFormat):
    @classmethod
    def check(cls, arg):
        # bounds checking for GF(2^n)???
        super(ImmediateGF2NAF, cls).check(arg)

class PlayerNoAF(IntArgFormat):
    @classmethod
    def check(cls, arg):
        super(PlayerNoAF, cls).check(arg)
        if arg > 256:
            raise ArgumentError(arg, 'Player number > 256')

class String(ArgFormat):
    length = 16

    @classmethod
    def check(cls, arg):
        if not isinstance(arg, str):
            raise ArgumentError(arg, 'Argument is not string')
        if len(arg) > cls.length:
            raise ArgumentError(arg, 'String longer than ' + cls.length)
        if '\0' in arg:
            raise ArgumentError(arg, 'String contains zero-byte')

    @classmethod
    def encode(cls, arg):
        return arg + '\0' * (cls.length - len(arg))

ArgFormats = {
    'c': ClearModpAF,
    's': SecretModpAF,
    'cw': ClearModpAF,
    'sw': SecretModpAF,
    'cg': ClearGF2NAF,
    'sg': SecretGF2NAF,
    'cgw': ClearGF2NAF,
    'sgw': SecretGF2NAF,
    'ci': ClearIntAF,
    'ciw': ClearIntAF,
    'i': ImmediateModpAF,
    'ig': ImmediateGF2NAF,
    'int': IntArgFormat,
    'p': PlayerNoAF,
    'str': String,
}

def format_str_is_reg(format_str):
    return issubclass(ArgFormats[format_str], RegisterArgFormat)

def format_str_is_writeable(format_str):
    return format_str_is_reg(format_str) and format_str[-1] == 'w'


class Instruction(object):
    """
    Base class for a RISC-type instruction. Has methods for checking arguments,
    getting byte encoding, emulating the instruction, etc.
    """
    __slots__ = ['args', 'arg_format', 'code', 'caller']
    count = 0

    def __init__(self, *args, **kwargs):
        """ Create an instruction and append it to the program list. """
        self.args = list(args)
        self.check_args()
        if not program.FIRST_PASS:
            if kwargs.get('add_to_prog', True):
                program.curr_block.instructions.append(self)
            if program.DEBUG:
                self.caller = [frame[1:] for frame in inspect.stack()[1:]]
            else:
                self.caller = None
            if program.EMULATE:
                self.execute()
        
        Instruction.count += 1
        if Instruction.count % 100000 == 0:
            print "Compiled %d lines at" % self.__class__.count, time.asctime()

    def get_code(self):
        return self.code
    
    def get_encoding(self):
        enc = int_to_bytes(self.get_code())
        # add the number of registers if instruction flagged as has var args
        if self.has_var_args():
            enc += int_to_bytes(len(self.args))
        for arg,format in zip(self.args, self.arg_format):
            enc += ArgFormats[format].encode(arg)
        return enc
    
    def get_bytes(self):
        return bytearray(self.get_encoding())
    
    def execute(self):
        """ Emulate execution of this instruction """
        raise NotImplementedError('execute method must be implemented')
    
    def check_args(self):
        """ Check the args match up with that specified in arg_format """
        for n,(arg,f) in enumerate(itertools.izip_longest(self.args, self.arg_format)):
            if arg is None:
                if not isinstance(self.arg_format, (list, tuple)):
                    break # end of optional arguments
                else:
                    raise CompilerError('Incorrect number of arguments for instruction %s' % (self))
            try:
                ArgFormats[f].check(arg)
            except ArgumentError as e:
                raise CompilerError('Invalid argument "%s" to instruction: %s'
                    % (e.arg, self) + '\n' + e.msg)
            except KeyError as e:
                raise CompilerError('Incorrect number of arguments for instruction %s' % (self))
    
    def get_used(self):
        """ Return the set of registers that are read in this instruction. """
        return set(arg for arg,w in zip(self.args, self.arg_format) if \
            format_str_is_reg(w) and not format_str_is_writeable(w))
    
    def get_def(self):
        """ Return the set of registers that are written to in this instruction. """
        return set(arg for arg,w in zip(self.args, self.arg_format) if \
            format_str_is_writeable(w))
    
    def get_pre_arg(self):
        return ""

    def has_var_args(self):
        return False

    def is_vec(self):
        return False

    def is_gf2n(self):
        return False

    def get_size(self):
        return 1

    def add_usage(self, req_node):
        pass

    # String version of instruction attempting to replicate encoded version
    def __str__(self):
        
        if self.has_var_args():
            varargCount = str(len(self.args)) + ', '
        else:
            varargCount = ''

        return self.__class__.__name__ + ' ' + self.get_pre_arg() + varargCount + ', '.join(str(a) for a in self.args)

    def __repr__(self):
        return self.__class__.__name__ + '(' + self.get_pre_arg() + ','.join(str(a) for a in self.args) + ')'

###
### Basic arithmetic
###

class AddBase(Instruction):
    __slots__ = []

    def execute(self):
        self.args[0].value = (self.args[1].value + self.args[2].value) % program.P

class SubBase(Instruction):
    __slots__ = []

    def execute(self):
        self.args[0].value = (self.args[1].value - self.args[2].value) % program.P

class MulBase(Instruction):
    __slots__ = []

    def execute(self):
        self.args[0].value = (self.args[1].value * self.args[2].value) % program.P

###
### Basic arithmetic with immediate values
###

class ImmediateBase(Instruction):
    __slots__ = ['op']

    def execute(self):
        exec('self.args[0].value = self.args[1].value.%s(self.args[2]) %% program.P' % self.op)

class SharedImmediate(ImmediateBase):
    __slots__ = []
    arg_format = ['sw', 's', 'i']

class ClearImmediate(ImmediateBase):
    __slots__ = []
    arg_format = ['cw', 'c', 'i']


###
### Memory access instructions
###

class DirectMemoryInstruction(Instruction):
    __slots__ = []
    def __init__(self, *args, **kwargs):
        super(DirectMemoryInstruction, self).__init__(*args, **kwargs)

class ReadMemoryInstruction(Instruction):
    __slots__ = []

class WriteMemoryInstruction(Instruction):
    __slots__ = []

class DirectMemoryWriteInstruction(DirectMemoryInstruction, \
                                       WriteMemoryInstruction):
    __slots__ = []
    def __init__(self, *args, **kwargs):
        if program.curr_tape.prevent_direct_memory_write:
            raise CompilerError('Direct memory writing prevented')
        super(DirectMemoryWriteInstruction, self).__init__(*args, **kwargs)

###
### I/O instructions
###

class DoNotEliminateInstruction(Instruction):
    """ What do you think? """
    __slots__ = []

class IOInstruction(DoNotEliminateInstruction):
    """ Instruction that uses stdin/stdout during runtime. These are linked
    to prevent instruction reordering during optimization. """
    __slots__ = []

    @classmethod
    def str_to_int(cls, s):
        """ Convert a 4 character string to an integer. """
        if len(s) > 4:
            raise CompilerError('String longer than 4 characters')
        n = 0
        for c in reversed(s.ljust(4)):
            n <<= 8
            n += ord(c)
        return n

class AsymmetricCommunicationInstruction(DoNotEliminateInstruction):
    """ Instructions involving sending from or to only one party. """
    __slots__ = []

class RawInputInstruction(AsymmetricCommunicationInstruction):
    """ Raw input instructions. """
    __slots__ = []

class PublicFileIOInstruction(DoNotEliminateInstruction):
    """ Instruction to reads/writes public information from/to files. """
    __slots__ = []

###
### Data access instructions
###

class DataInstruction(Instruction):
    __slots__ = []
    field_type = 'modp'

    def add_usage(self, req_node):
        req_node.increment((self.field_type, self.data_type), self.get_size())

###
### Integer operations
### 

class IntegerInstruction(Instruction):
    """ Base class for integer operations. """
    __slots__ = []
    arg_format = ['ciw', 'ci', 'ci']

class StackInstruction(Instruction):
    """ Base class for thread-local stack instructions. """
    __slots__ = []

###
### Clear comparison instructions
###

class UnaryComparisonInstruction(Instruction):
    """ Base class for unary comparisons. """
    __slots__ = []
    arg_format = ['ciw', 'ci']

### 
### Clear shift instructions
### 

class ClearShiftInstruction(ClearImmediate):
    __slots__ = []

    def check_args(self):
        super(ClearShiftInstruction, self).check_args()
        if program.galois_length > 64:
            bits = 127
        else:
            # assume 64-bit machine
            bits = 63
        if self.args[2] > bits:
            raise CompilerError('Shifting by more than %d bits '
                                'not implemented' % bits)

###
### Jumps etc
###

class dummywrite(Instruction):
    """ Dummy instruction to create source node in the dependency graph,
        preventing read-before-write warnings. """
    __slots__ = []
    
    def __init__(self, *args, **kwargs):
        self.arg_format = [arg.reg_type + 'w' for arg in args]
        super(dummywrite, self).__init__(*args, **kwargs)
    
    def execute(self):
        pass
    
    def get_encoding(self):
        return []

class JumpInstruction(Instruction):
    __slots__ = ['jump_arg']

    def set_relative_jump(self, value):
        if value == -1:
            raise CompilerException('Jump by -1 would cause infinite loop')
        self.args[self.jump_arg] = value

    def get_relative_jump(self):
        return self.args[self.jump_arg]


class VarArgsInstruction(Instruction):
    def has_var_args(self):
        return True


class CISC(Instruction):
    """
    Base class for a CISC instruction.
    
    Children must implement expand(self) to process the instruction.
    """
    __slots__ = []
    code = None

    def __init__(self, *args):
        self.args = args
        self.check_args()
        #if EMULATE:
        #    self.expand()
        if not program.FIRST_PASS:
            self.expand()
    
    def expand(self):
        """ Expand this into a sequence of RISC instructions. """
        raise NotImplementedError('expand method must be implemented')
