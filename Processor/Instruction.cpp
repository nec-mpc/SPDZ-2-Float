// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


#include "Processor/Instruction.h"
#include "Processor/Machine.h"
#include "Processor/Processor.h"
#include "Exceptions/Exceptions.h"
#include "Tools/time-func.h"
#include "Tools/parse.h"

#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <map>
#include <math.h>
#include <valgrind/callgrind.h>
#include <bitset>


// broken
#undef DEBUG

void Instruction::parse(istream& s)
{
  n=0; start.resize(0);
  r[0]=0; r[1]=0; r[2]=0; r[3]=0;

  int pos=s.tellg();
  opcode=get_int(s);
  size=opcode>>9;
//  opcode&=0x1FF;
  opcode&=0x3FF;
  
  if (size==0)
    size=1;

//  cout << "opcode = " << hex << opcode << dec << endl;

  parse_operands(s, pos);

}


void BaseInstruction::parse_operands(istream& s, int pos)
{
  int num_var_args = 0;
  mp::uint256_t tmp = 0;
  switch (opcode)
  {
      // instructions with 3 register operands
      case ADDC:
      case ADDS:
      case ADDM:
      case SUBC:
      case SUBS:
      case SUBML:
      case SUBMR:
      case MULC:
      case MULM:
      case DIVC:
      case MODC:
      case TRIPLE:
      case ANDC:
      case XORC:
      case ORC:
      case SHLC:
      case SHRC:
      case GADDC:
      case GADDS:
      case GADDM:
      case GSUBC:
      case GSUBS:
      case GSUBML:
      case GSUBMR:
      case GMULC:
      case GMULM:
      case GDIVC:
      case GTRIPLE:
      case GBITTRIPLE:
      case GBITGF2NTRIPLE:
      case GANDC:
      case GXORC:
      case GORC:
      case GMULBITC:
      case GMULBITM:
      case LTC:
      case GTC:
      case EQC:
      case ADDINT:
      case SUBINT:
      case MULINT:
      case DIVINT:
      case E_MP_ADDS:
      case E_MP_ADDM:
      case E_MP_SUBS:
      case E_MP_SUBMR:
      case E_MP_SUBML:
      case E_MP_MULM:
        r[0]=get_int(s);
        r[1]=get_int(s);
        r[2]=get_int(s);
        break;
      // instructions with 2 register operands
      case LDMCI:
      case LDMSI:
      case STMCI:
      case STMSI:
      case MOVC:
      case MOVS:
      case MOVINT:
      case LDMINTI:
      case STMINTI:
      case LEGENDREC:
      case SQUARE:
      case INV:
      case GINV:
      case CONVINT:
      case GLDMCI:
      case GLDMSI:
      case GSTMCI:
      case GSTMSI:
      case GMOVC:
      case GMOVS:
      case GSQUARE:
      case GNOTC:
      case GCONVINT:
      case GCONVGF2N:
      case LTZC:
      case EQZC:
      case RAND:
      case PROTECTMEMS:
      case PROTECTMEMC:
      case GPROTECTMEMS:
      case GPROTECTMEMC:
      case PROTECTMEMINT:
        r[0]=get_int(s);
        r[1]=get_int(s);
        break;
      // instructions with 1 register operand
      case BIT:
      case PRINTMEM:
      case PRINTREGPLAIN:
      case LDTN:
      case LDARG:
      case STARG:
      case JMPI:
      case GBIT:
      case GPRINTMEM:
      case GPRINTREGPLAIN:
      case JOIN_TAPE:
      case PUSHINT:
      case POPINT:
      case PUBINPUT:
      case RAWOUTPUT:
      case GRAWOUTPUT:
      case PRINTCHRINT:
      case PRINTSTRINT:
      case PRINTINT:
        r[0]=get_int(s);
        break;
      // instructions with 3 registers + 1 integer operand
        r[0]=get_int(s);
        r[1]=get_int(s);
        r[2]=get_int(s);
        n = get_int(s);
        break;
      // instructions with 2 registers + 1 integer operand
      case ADDCI:
      case ADDSI:
      case SUBCI:
      case SUBSI:
      case SUBCFI:
      case SUBSFI:
      case MULCI:
      case MULSI:
      case DIVCI:
      case MODCI:
      case ANDCI:
      case XORCI:
      case ORCI:
      case SHLCI:
      case SHRCI:
      case NOTC:
      case CONVMODP:
      case GADDCI:
      case GADDSI:
      case GSUBCI:
      case GSUBSI:
      case GSUBCFI:
      case GSUBSFI:
      case GMULCI:
      case GMULSI:
      case GDIVCI:
      case GANDCI:
      case GXORCI:
      case GORCI:
      case GSHLCI:
      case GSHRCI:
      case USE:
      case USE_INP:
      case RUN_TAPE:
      case STARTPRIVATEOUTPUT:
      case GSTARTPRIVATEOUTPUT:
      case DIGESTC:
      case CONNECTIPV4: // write socket handle, read IPv4 address, portnum
        r[0]=get_int(s);
        r[1]=get_int(s);
        n = get_int(s);
        break;
      // instructions with 1 register + 1 integer operand
      case LDI:
      case LDSI:
      case LDMC:
      case LDMS:
      case STMC:
      case STMS:
      case LDMINT:
      case STMINT:
      case INPUT:
      case JMPNZ:
      case JMPEQZ:
      case GLDI:
      case GLDSI:
      case GLDMC:
      case GLDMS:
      case GSTMC:
      case GSTMS:
      case GINPUT:
      case PRINTREG:
      case GPRINTREG:
      case LDINT:
      case STARTINPUT:
      case GSTARTINPUT:
      case STOPPRIVATEOUTPUT:
      case GSTOPPRIVATEOUTPUT:
      case INPUTMASK:
      case GINPUTMASK:
      case ACCEPTCLIENTCONNECTION:
        r[0]=get_int(s);
        n = get_int(s);
        break;
      // instructions with 1 integer operand
      case PRINTSTR:
      case PRINTCHR:
      case JMP:
      case START:
      case STOP:
      case LISTEN:
        n = get_int(s);
        break;
      // instructions with no operand
      case TIME:
      case CRASH:
      case STARTGRIND:
      case STOPGRIND:
      	break;
      // instructions with 4 register operands
      case PRINTFLOATPLAIN:
        get_vector(4, start, s);
        break;
      // open instructions + read/write instructions with variable length args
      case STARTOPEN:
      case STOPOPEN:
      case GSTARTOPEN:
      case GSTOPOPEN:
      case WRITEFILESHARE:
      case OPEN:
      case GOPEN:
#if defined(EXTENDED_SPDZ_GFP)
      case E_STARTMULT:
      case E_STOPMULT:
      case E_MULT:
      case E_MP_MULT:
      case E_MP_OPEN:
#endif
#if defined(EXTENDED_SPDZ_GF2N)
      case GE_MULT:
#endif
        num_var_args = get_int(s);
        get_vector(num_var_args, start, s);
        break;

      // read from file, input is opcode num_args, 
      //   start_file_posn (read), end_file_posn(write) var1, var2, ...
      case READFILESHARE:
        num_var_args = get_int(s) - 2;
        r[0] = get_int(s);
        r[1] = get_int(s);
        get_vector(num_var_args, start, s);
        break;

      // read from external client, input is : opcode num_args, client_id, var1, var2 ...
      case READSOCKETC:
      case READSOCKETS:
      case READSOCKETINT:
      case READCLIENTPUBLICKEY:   
        num_var_args = get_int(s) - 1;
        r[0] = get_int(s);
        get_vector(num_var_args, start, s);
        break;

      // write to external client, input is : opcode num_args, client_id, message_type, var1, var2 ...
      case WRITESOCKETC:
      case WRITESOCKETS:
      case WRITESOCKETSHARE:
      case WRITESOCKETINT:
        num_var_args = get_int(s) - 2;
        r[0] = get_int(s);
        r[1] = get_int(s);
        get_vector(num_var_args, start, s);
        break;
      case INITSECURESOCKET:
      case RESPSECURESOCKET:
        num_var_args = get_int(s) - 1;
        r[0] = get_int(s);
        get_vector(num_var_args, start, s);
        break;
      // raw input
      case STOPINPUT:
      case GSTOPINPUT:
        // subtract player number argument
        num_var_args = get_int(s) - 1;
        n = get_int(s);
        get_vector(num_var_args, start, s);
        break;
      case GBITDEC:
      case GBITCOM:
        num_var_args = get_int(s) - 2;
        r[0] = get_int(s);
        n = get_int(s);
        get_vector(num_var_args, start, s);
        break;
      case BITDECINT:
          num_var_args = get_int(s) - 1;
          r[0] = get_int(s);
          get_vector(num_var_args, start, s);
          break;
      case PREP:
      case GPREP:
        // subtract extra argument
        num_var_args = get_int(s) - 1;
        s.read((char*)r, sizeof(r));
        start.resize(num_var_args);
        for (int i = 0; i < num_var_args; i++)
        { start[i] = get_int(s); }
        break;
      case USE_PREP:
      case GUSE_PREP:
        s.read((char*)r, sizeof(r));
        n = get_int(s);
        break;
      case REQBL:
        n = get_int(s);
        break;
      case GREQBL:
        n = get_int(s);
        /*
        if (n > 0 && gf2n::degree() < int(n))
          {
            stringstream ss;
            ss << "Tape requires prime of bit length " << n << endl;
            throw Processor_Error(ss.str());
          }
          */
        break;
#if defined(EXTENDED_SPDZ_Z2N)
      case E_PRINTFIXEDPLAIN:
    	  r[0] = get_int(s);
    	  n = get_int(s);
    	  break;
      case E_SKEW_BIT_DEC:
      case E_MP_SKEW_BIT_DEC:
    	  r[0] = get_int(s);
    	  n = get_int(s);
    	  num_var_args = n * 3;
    	  get_vector(num_var_args, start, s);
    	  break;
      case E_SKEW_BIT_INJ:
      case E_SKEW_BIT_REC:
      case E_MP_SKEW_BIT_INJ:
         r[0] = get_int(s);
         get_vector(3, start, s);
         break;
      case E_SKEW_RING_REC:
      case E_MP_SKEW_RING_REC:
    	  r[0] = get_int(s);
    	  num_var_args = get_int(s);
    	  get_vector(num_var_args, start, s);
    	  break;
      case E_MP_LDSI:
      case E_MP_LDI:
         r[0]=get_int(s);
         mp_n=0;
         for(size_t i=0; i<8; ++i)
          {
        	tmp = (mp::uint256_t) get_int(s);
        	mp_n +=  tmp << (31 * i);
          }
         tmp = (mp::uint256_t) get_int(s);
         mp_n +=  tmp << 248;
         break;
      case E_MP_ADDSI:
      case E_MP_MULSI:
    	  r[0]=get_int(s);
    	  r[1]=get_int(s);
    	  mp_n=0;
    	  for(size_t i=0; i<8; ++i) {
    		  tmp = (mp::uint256_t) get_int(s);
    		  mp_n +=  tmp << (31 * i);
    	  }
    	  tmp = (mp::uint256_t) get_int(s);
    	  mp_n +=  tmp << 248;
    	  break;
      case E_MP_SUBSIL:
               r[0]=get_int(s);
               r[1]=get_int(s);
               mp_n=0;
               for(size_t i=0; i<8; ++i) {
                	tmp = (mp::uint256_t) get_int(s);
                	mp_n +=  tmp << (31 * i);
                }
               tmp = (mp::uint256_t) get_int(s);
               mp_n +=  tmp << 248;
               break;
      case E_MP_SUBSIR:
         r[0]=get_int(s);
         mp_n=0;
         for(size_t i=0; i<8; ++i) {
          	tmp = (mp::uint256_t) get_int(s);
          	mp_n +=  tmp << (31 * i);
          }
         tmp = (mp::uint256_t) get_int(s);
         mp_n +=  tmp << 248;
         r[1]=get_int(s);
         break;

#endif
      default:
        ostringstream os;
        os << "Invalid instruction " << hex << showbase << opcode << " at " << dec << pos;
        throw Invalid_Instruction(os.str());
  }
}


bool Instruction::get_offline_data_usage(DataPositions& usage)
{
  switch (opcode)
  {
    case USE:
      if (r[0] >= N_DATA_FIELD_TYPE)
        throw invalid_program();
      if (r[1] >= N_DTYPE)
        throw invalid_program();
      usage.files[r[0]][r[1]] = n;
      return int(n) >= 0;
    case USE_INP:
      if (r[0] >= N_DATA_FIELD_TYPE)
        throw invalid_program();
      if ((unsigned)r[1] >= usage.inputs.size())
        throw Processor_Error("Player number too high");
      usage.inputs[r[1]][r[0]] = n;
      return int(n) >= 0;
    case USE_PREP:
      usage.extended[gfp::field_type()][r] = n;
      return int(n) >= 0;
    case GUSE_PREP:
      usage.extended[gf2n::field_type()][r] = n;
      return int(n) >= 0;
    default:
      return true;
  }
}

int BaseInstruction::get_reg_type() const
{
  switch (opcode) {
    case LDMINT:
    case STMINT:
    case LDMINTI:
    case STMINTI:
    case PUSHINT:
    case POPINT:
    case MOVINT:
    case READSOCKETINT:
    case WRITESOCKETINT:
    case READCLIENTPUBLICKEY:
    case INITSECURESOCKET:
    case RESPSECURESOCKET:
    case LDARG:
    case LDINT:
    case CONVMODP:
    case GCONVGF2N:
    case RAND:
      return INT;
    case PREP:
    case USE_PREP:
    case GUSE_PREP:
      // those use r[] for a string
      return NONE;
    default:
      if (is_gf2n_instruction())
        return GF2N;
      else if (opcode >> 4 == 0x9)
        return INT;
      else
        return MODP;
  }
}

int BaseInstruction::get_max_reg(int reg_type) const
{
  if (get_reg_type() != reg_type) { return 0; }

  if (start.size())
    return *max_element(start.begin(), start.end()) + size;
  else
    return *max_element(r, r + 3) +  size;
}

int Instruction::get_mem(RegType reg_type, SecrecyType sec_type) const
{
  if (get_reg_type() == reg_type and is_direct_memory_access(sec_type))
    return n + size;
  else
    return 0;
}

bool BaseInstruction::is_direct_memory_access(SecrecyType sec_type) const
{
  if (sec_type == SECRET)
  {
    switch (opcode)
    {
      case LDMS:
      case STMS:
      case GLDMS:
      case GSTMS:
        return true;
      default:
        return false;
    }
  }
  else
  {
    switch (opcode)
    {
      case LDMC:
      case STMC:
      case GLDMC:
      case GSTMC:
      case LDMINT:
      case STMINT:
        return true;
      default:
        return false;
    }
  }
}



ostream& operator<<(ostream& s,const Instruction& instr)
{
  s << instr.opcode << " : ";
  for (int i=0; i<3; i++)
    { s << instr.r[i] << " "; }
  s << " : " << instr.n;
  if (instr.start.size()!=0)
    { s << " : " << instr.start.size() << " : ";
      for (unsigned int i=0; i<instr.start.size(); i++)
	{ s << instr.start[i] << " "; }
    }
  return s;
} 


void Instruction::execute(Processor& Proc) const
{
  Proc.PC+=1;

#ifndef EXTENDED_SPDZ_GF2N
#ifndef DEBUG
  // optimize some instructions
  switch (opcode)
  {
    case GADDC:
      for (int i = 0; i < size; i++)
        Proc.get_C2_ref(r[0] + i).add(Proc.read_C2(r[1] + i),Proc.read_C2(r[2] + i));
      return;
    case GADDS:
      for (int i = 0; i < size; i++)
        Proc.get_S2_ref(r[0] + i).add(Proc.read_S2(r[1] + i),Proc.read_S2(r[2] + i));
      return;
    case GMOVC:
      for (int i = 0; i < size; i++)
        Proc.write_C2(r[0] + i, Proc.read_C2(r[1] + i));
      return;
    case GANDC:
      for (int i = 0; i < size; i++)
         Proc.get_C2_ref(r[0] + i).AND(Proc.read_C2(r[1] + i),Proc.read_C2(r[2] + i));
      return;
    case GSHLCI:
      for (int i = 0; i < size; i++)
         Proc.get_C2_ref(r[0] + i).SHL(Proc.read_C2(r[1] + i),n);
      return;
    case GSHRCI:
      for (int i = 0; i < size; i++)
        Proc.get_C2_ref(r[0] + i).SHR(Proc.read_C2(r[1] + i),n);
      return;
    case GMULM:
      for (int i = 0; i < size; i++)
         Proc.get_S2_ref(r[0] + i).mul(Proc.read_S2(r[1] + i),Proc.read_C2(r[2] + i));
      return;
  }
#endif
#endif

  int r[3] = {this->r[0], this->r[1], this->r[2]};
  int n = this->n;
  for (int i = 0; i < size; i++) 
  { switch (opcode)
    { case LDI: 
        Proc.temp.ansp.assign(n);
        Proc.write_Cp(r[0],Proc.temp.ansp);
        break;
      case E_MP_LDI:
            Proc.temp.ansp.assign(mp_n);
            Proc.write_Cp(r[0],Proc.temp.ansp);
            break;
      case GLDI: 
        Proc.temp.ans2.assign(n);
        Proc.write_C2(r[0],Proc.temp.ans2);
        break;
      case LDSI:
    	  Proc.temp.ansp.assign(n);
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.PLdsi_Ext(Proc.temp.ansp, Proc.get_Sp_ref(r[0]));
#else
        { if (Proc.P.my_num()==0)
            Proc.get_Sp_ref(r[0]).set_share(Proc.temp.ansp);
          else
            Proc.get_Sp_ref(r[0]).assign_zero();
          gfp& tmp=Proc.temp.tmpp;
          tmp.mul(Proc.MCp.get_alphai(),Proc.temp.ansp);
          Proc.get_Sp_ref(r[0]).set_mac(tmp);
        }
#endif
        break;
      case GLDSI:
    	Proc.temp.ans2.assign(n);
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GLdsi_Ext(Proc.temp.ans2, Proc.get_S2_ref(r[0]));
#else
        {
          if (Proc.P.my_num()==0)
            Proc.get_S2_ref(r[0]).set_share(Proc.temp.ans2);
          else
            Proc.get_S2_ref(r[0]).assign_zero();
          gf2n& tmp=Proc.temp.tmp2;
          tmp.mul(Proc.MC2.get_alphai(),Proc.temp.ans2);
          Proc.get_S2_ref(r[0]).set_mac(tmp);
        }
#endif
        break;
      case E_MP_LDSI:
          	  Proc.temp.ansp.assign(mp_n);
          	  Proc.MPLdsi_Ext(Proc.temp.ansp, Proc.get_Sp_ref(r[0]));
      break;
      case LDMC:
        Proc.write_Cp(r[0],Proc.machine.Mp.read_C(n));
        n++;
        break;
      case GLDMC:
        Proc.write_C2(r[0],Proc.machine.M2.read_C(n));
        n++;
        break;
      case LDMS:
        Proc.write_Sp(r[0],Proc.machine.Mp.read_S(n));
        n++;
        break;
      case GLDMS:
        Proc.write_S2(r[0],Proc.machine.M2.read_S(n));
        n++;
        break;
      case LDMINT:
        Proc.write_Ci(r[0],Proc.machine.Mi.read_C(n).get());
        n++;
        break;
      case LDMCI:
        Proc.write_Cp(r[0], Proc.machine.Mp.read_C(Proc.read_Ci(r[1])));
        break;
      case GLDMCI:
        Proc.write_C2(r[0], Proc.machine.M2.read_C(Proc.read_Ci(r[1])));
        break;
      case LDMSI:
        Proc.write_Sp(r[0], Proc.machine.Mp.read_S(Proc.read_Ci(r[1])));
        break;
      case GLDMSI:
        Proc.write_S2(r[0], Proc.machine.M2.read_S(Proc.read_Ci(r[1])));
        break;
      case LDMINTI:
        Proc.write_Ci(r[0], Proc.machine.Mi.read_C(Proc.read_Ci(r[1])).get());
        break;
      case STMC:
        Proc.machine.Mp.write_C(n,Proc.read_Cp(r[0]),Proc.PC);
        n++;
        break;
      case GSTMC:
        Proc.machine.M2.write_C(n,Proc.read_C2(r[0]),Proc.PC);
        n++;
        break;
      case STMS:
        Proc.machine.Mp.write_S(n,Proc.read_Sp(r[0]),Proc.PC);
        n++;
        break;
      case GSTMS:
        Proc.machine.M2.write_S(n,Proc.read_S2(r[0]),Proc.PC);
        n++;
        break;
      case STMINT:
        Proc.machine.Mi.write_C(n,Integer(Proc.read_Ci(r[0])),Proc.PC);
        n++;
        break;
      case STMCI:
        Proc.machine.Mp.write_C(Proc.read_Ci(r[1]), Proc.read_Cp(r[0]),Proc.PC);
        break;
      case GSTMCI:
        Proc.machine.M2.write_C(Proc.read_Ci(r[1]), Proc.read_C2(r[0]),Proc.PC);
        break;
      case STMSI:
        Proc.machine.Mp.write_S(Proc.read_Ci(r[1]), Proc.read_Sp(r[0]),Proc.PC);
        break;
      case GSTMSI:
        Proc.machine.M2.write_S(Proc.read_Ci(r[1]), Proc.read_S2(r[0]),Proc.PC);
        break;
      case STMINTI:
        Proc.machine.Mi.write_C(Proc.read_Ci(r[1]), Integer(Proc.read_Ci(r[0])),Proc.PC);
        break;
      case MOVC:
        Proc.write_Cp(r[0],Proc.read_Cp(r[1]));
        break;
      case GMOVC:
        Proc.write_C2(r[0],Proc.read_C2(r[1]));
        break;
      case MOVS:
#if defined(EXTENDED_SPDZ_GFP)
    	Proc.PMovs_Ext(Proc.get_Sp_ref(r[0]), Proc.get_Sp_ref(r[1]));
#else
        Proc.write_Sp(r[0],Proc.read_Sp(r[1]));
#endif
        break;
      case GMOVS:
#if defined(EXTENDED_SPDZ_GF2N)
    	Proc.GMovs_Ext(Proc.get_S2_ref(r[0]), Proc.get_S2_ref(r[1]));
#else
        Proc.write_S2(r[0],Proc.read_S2(r[1]));
#endif
        break;
      case MOVINT:
        Proc.write_Ci(r[0],Proc.read_Ci(r[1]));
        break;
      case PROTECTMEMS:
        Proc.machine.Mp.protect_s(Proc.read_Ci(r[0]), Proc.read_Ci(r[1]));
        break;
      case PROTECTMEMC:
        Proc.machine.Mp.protect_c(Proc.read_Ci(r[0]), Proc.read_Ci(r[1]));
        break;
      case GPROTECTMEMS:
        Proc.machine.M2.protect_s(Proc.read_Ci(r[0]), Proc.read_Ci(r[1]));
        break;
      case GPROTECTMEMC:
        Proc.machine.M2.protect_c(Proc.read_Ci(r[0]), Proc.read_Ci(r[1]));
        break;
      case PROTECTMEMINT:
        Proc.machine.Mi.protect_c(Proc.read_Ci(r[0]), Proc.read_Ci(r[1]));
        break;
      case PUSHINT:
        Proc.pushi(Proc.read_Ci(r[0]));
        break;
      case POPINT:
        Proc.popi(Proc.get_Ci_ref(r[0]));
        break;
      case LDTN:
        Proc.write_Ci(r[0],Proc.get_thread_num());
        break;
      case LDARG:
        Proc.write_Ci(r[0],Proc.get_arg());
        break;
      case STARG:
        Proc.set_arg(Proc.read_Ci(r[0]));
        break;
      case ADDC:
	#ifdef DEBUG
           Proc.temp.ansp.add(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
           Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).add(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	#endif
        break;
      case GADDC:
	#ifdef DEBUG
           ans2.add(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).add(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	#endif
        break;
      case ADDS:
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.PAdds_Ext(Proc.get_Sp_ref(r[0]), Proc.get_Sp_ref(r[1]), Proc.get_Sp_ref(r[2]));
#else
		#ifdef DEBUG
           Sansp.add(Proc.read_Sp(r[1]),Proc.read_Sp(r[2]));
           Proc.write_Sp(r[0],Sansp);
        #else
           Proc.get_Sp_ref(r[0]).add(Proc.read_Sp(r[1]),Proc.read_Sp(r[2]));
        #endif
#endif
        break;
      case E_MP_ADDS:
         Proc.MPAdds_Ext(Proc.get_Sp_ref(r[0]), Proc.get_Sp_ref(r[1]), Proc.get_Sp_ref(r[2]));
        break;
      case GADDS:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GAdds_Ext(Proc.get_S2_ref(r[0]), Proc.get_S2_ref(r[1]), Proc.get_S2_ref(r[2]));
#else
	#ifdef DEBUG
           Sans2.add(Proc.read_S2(r[1]),Proc.read_S2(r[2]));
           Proc.write_S2(r[0],Sans2);
        #else
           Proc.get_S2_ref(r[0]).add(Proc.read_S2(r[1]),Proc.read_S2(r[2]));
        #endif
#endif
        break;
      case ADDM:
#if defined(EXTENDED_SPDZ_GFP)
	   Proc.PAddm_Ext(Proc.get_Sp_ref(r[1]), Proc.get_Cp_ref(r[2]), Proc.get_Sp_ref(r[0]));
#else
       #ifdef DEBUG
           Sansp.add(Proc.read_Sp(r[1]),Proc.read_Cp(r[2]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
           Proc.write_Sp(r[0],Sansp);
       #else
           Proc.get_Sp_ref(r[0]).add(Proc.read_Sp(r[1]),Proc.read_Cp(r[2]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
       #endif
#endif
        break;
      case E_MP_ADDM:
      	 Proc.MPAddm_Ext(Proc.get_Sp_ref(r[1]), Proc.get_Cp_ref(r[2]), Proc.get_Sp_ref(r[0]));
        break;
      case GADDM:
#if defined(EXTENDED_SPDZ_GF2N)
	   Proc.GAddm_Ext(Proc.get_S2_ref(r[1]), Proc.get_C2_ref(r[2]), Proc.get_S2_ref(r[0]));
#else
        #ifdef DEBUG
           Sans2.add(Proc.read_S2(r[1]),Proc.read_C2(r[2]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
	   Proc.write_S2(r[0],Sans2);
        #else
           Proc.get_S2_ref(r[0]).add(Proc.read_S2(r[1]),Proc.read_C2(r[2]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
        #endif
#endif
        break;
      case SUBC:
	#ifdef DEBUG
          ansp.sub(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
          Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
          Proc.get_Cp_ref(r[0]).sub(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	#endif
        break;
      case GSUBC:
	#ifdef DEBUG
          Proc.temp.ans2.sub(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
          Proc.write_C2(r[0],Proc.temp.ans2);
	#else
          Proc.get_C2_ref(r[0]).sub(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	#endif
        break;
      case SUBS:
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.PSubs_Ext(Proc.get_Sp_ref(r[0]), Proc.get_Sp_ref(r[1]), Proc.get_Sp_ref(r[2]));
#else
	#ifdef DEBUG
    	  Sansp.sub(Proc.read_Sp(r[1]),Proc.read_Sp(r[2]));
          Proc.write_Sp(r[0],Sansp);
    #else
           Proc.get_Sp_ref(r[0]).sub(Proc.read_Sp(r[1]),Proc.read_Sp(r[2]));
	#endif
#endif
        break;
      case E_MP_SUBS:
          	Proc.MPSubs_Ext(Proc.get_Sp_ref(r[0]), Proc.get_Sp_ref(r[1]), Proc.get_Sp_ref(r[2]));
        break;
      case GSUBS:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GSubs_Ext(Proc.get_S2_ref(r[0]), Proc.get_S2_ref(r[1]), Proc.get_S2_ref(r[2]));
#else
	#ifdef DEBUG
           Sans2.sub(Proc.read_S2(r[1]),Proc.read_S2(r[2]));
	   Proc.write_S2(r[0],Sans2);
        #else
           Proc.get_S2_ref(r[0]).sub(Proc.read_S2(r[1]),Proc.read_S2(r[2]));
	#endif
#endif
        break;
      case SUBML:
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.PSubml_Ext(Proc.get_Sp_ref(r[1]), Proc.get_Cp_ref(r[2]), Proc.get_Sp_ref(r[0]));
#else
   	   #ifdef DEBUG
           Sansp.sub(Proc.read_Sp(r[1]),Proc.read_Cp(r[2]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
           Proc.write_Sp(r[0],Sansp);
        #else
           Proc.get_Sp_ref(r[0]).sub(Proc.read_Sp(r[1]),Proc.read_Cp(r[2]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
        #endif
#endif
        break;
      case E_MP_SUBML:
          	Proc.MPSubml_Ext(Proc.get_Sp_ref(r[1]), Proc.get_Cp_ref(r[2]), Proc.get_Sp_ref(r[0]));
        break;
      case GSUBML:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GSubml_Ext(Proc.get_S2_ref(r[1]), Proc.get_C2_ref(r[2]), Proc.get_S2_ref(r[0]));
#else
	#ifdef DEBUG
    	  Sans2.sub(Proc.read_S2(r[1]),Proc.read_C2(r[2]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
    	  Proc.write_S2(r[0],Sans2);
	#else
    	  Proc.get_S2_ref(r[0]).sub(Proc.read_S2(r[1]),Proc.read_C2(r[2]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
	#endif
#endif
        break;
      case SUBMR:
#if defined(EXTENDED_SPDZ_GFP)
	   	   Proc.PSubmr_Ext(Proc.get_Cp_ref(r[1]), Proc.get_Sp_ref(r[2]), Proc.get_Sp_ref(r[0]));
#else
        #ifdef DEBUG
           Sansp.sub(Proc.read_Cp(r[1]),Proc.read_Sp(r[2]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
           Proc.write_Sp(r[0],Sansp);
        #else
           Proc.get_Sp_ref(r[0]).sub(Proc.read_Cp(r[1]),Proc.read_Sp(r[2]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
	#endif
#endif
        break;
      case E_MP_SUBMR:
      	   	Proc.MPSubmr_Ext(Proc.get_Cp_ref(r[1]), Proc.get_Sp_ref(r[2]), Proc.get_Sp_ref(r[0]));
        break;
      case GSUBMR:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GSubmr_Ext(Proc.get_C2_ref(r[1]), Proc.get_S2_ref(r[2]), Proc.get_S2_ref(r[0]));
#else
        #ifdef DEBUG
           Sans2.sub(Proc.read_C2(r[1]),Proc.read_S2(r[2]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
	   Proc.write_S2(r[0],Sans2);
        #else
           Proc.get_S2_ref(r[0]).sub(Proc.read_C2(r[1]),Proc.read_S2(r[2]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
		#endif
#endif
        break;
      case MULC:
	#ifdef DEBUG
          ansp.mul(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
          Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
          Proc.get_Cp_ref(r[0]).mul(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	#endif
        break;
      case GMULC:
	#ifdef DEBUG
          Proc.temp.ans2.mul(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
          Proc.write_C2(r[0],Proc.temp.ans2);
	#else
          Proc.get_C2_ref(r[0]).mul(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	#endif
        break;
      case MULM:
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.PMulm_Ext(Proc.get_Sp_ref(r[0]), Proc.read_Sp(r[1]), Proc.read_Cp(r[2]));
#else
	#ifdef DEBUG
       Sansp.mul(Proc.read_Sp(r[1]),Proc.read_Cp(r[2]));
	   Proc.write_Sp(r[0],Sansp);
	#else
           Proc.get_Sp_ref(r[0]).mul(Proc.read_Sp(r[1]),Proc.read_Cp(r[2]));
	#endif
#endif
        break;
      case E_MP_MULM:
          	Proc.MPMulm_Ext(Proc.get_Sp_ref(r[0]), Proc.read_Sp(r[1]), Proc.read_Cp(r[2]));
        break;
      case GMULM:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GMulm_Ext(Proc.get_S2_ref(r[0]), Proc.read_S2(r[1]), Proc.read_C2(r[2]));
#else
	#ifdef DEBUG
    	  Sans2.mul(Proc.read_S2(r[1]),Proc.read_C2(r[2]));
    	  Proc.write_S2(r[0],Sans2);
	#else
    	  Proc.get_S2_ref(r[0]).mul(Proc.read_S2(r[1]),Proc.read_C2(r[2]));
	#endif
#endif
        break;
      case DIVC:
        if (Proc.read_Cp(r[2]).is_zero())
          throw Processor_Error("Division by zero from register");
        Proc.temp.ansp.invert(Proc.read_Cp(r[2]));
        Proc.temp.ansp.mul(Proc.read_Cp(r[1]));
        Proc.write_Cp(r[0],Proc.temp.ansp);
        break;
      case GDIVC:
        if (Proc.read_C2(r[2]).is_zero())
          throw Processor_Error("Division by zero from register");
        Proc.temp.ans2.invert(Proc.read_C2(r[2]));
        Proc.temp.ans2.mul(Proc.read_C2(r[1]));
        Proc.write_C2(r[0],Proc.temp.ans2);
        break;
      case MODC:
        exit(EBADR);
        break;
      case LEGENDREC:
    	 exit(EBADR);
        break;
      case DIGESTC:
    	 exit(EBADR);
        break;
      case DIVCI:
        if (n == 0)
          throw Processor_Error("Division by immediate zero");
        Proc.temp.ansp.assign(n);
        Proc.temp.ansp.invert();
        Proc.temp.ansp.mul(Proc.read_Cp(r[1]));
        Proc.write_Cp(r[0],Proc.temp.ansp);
        break;
      case GDIVCI:
        if (n == 0)
          throw Processor_Error("Division by immediate zero");
        Proc.temp.ans2.assign(n);
        Proc.temp.ans2.invert();
        Proc.temp.ans2.mul(Proc.read_C2(r[1]));
        Proc.write_C2(r[0],Proc.temp.ans2);
        break;
      case MODCI:
    	 exit(EBADR);
        break;
      case GMULBITC:
  #ifdef DEBUG
          Proc.temp.ans2.mul_by_bit(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
          Proc.write_C2(r[0],Proc.temp.ans2);
  #else
          Proc.get_C2_ref(r[0]).mul_by_bit(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
  #endif
        break;
      case GMULBITM:
  #ifdef DEBUG
          Sans2.mul_by_bit(Proc.read_S2(r[1]),Proc.read_C2(r[2]));
          Proc.write_S2(r[0],Sans2);
  #else
          Proc.get_S2_ref(r[0]).mul_by_bit(Proc.read_S2(r[1]),Proc.read_C2(r[2]));
  #endif
        break;
      case ADDCI:
        Proc.temp.ansp.assign(n);
	#ifdef DEBUG
           Proc.temp.ansp.add(Proc.temp.ansp,Proc.read_Cp(r[1]));
           Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).add(Proc.temp.ansp,Proc.read_Cp(r[1]));
	#endif
        break;
      case GADDCI:
        Proc.temp.ans2.assign(n);
	#ifdef DEBUG
           Proc.temp.ans2.add(Proc.temp.ans2,Proc.read_C2(r[1]));
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).add(Proc.temp.ans2,Proc.read_C2(r[1]));
	#endif
        break;
      case ADDSI:
        Proc.temp.ansp.assign(n);
#if defined(EXTENDED_SPDZ_GFP)
        Proc.PAddm_Ext(Proc.get_Sp_ref(r[1]), Proc.temp.ansp, Proc.get_Sp_ref(r[0]));
#else
	#ifdef DEBUG
        Sansp.add(Proc.read_Sp(r[1]),Proc.temp.ansp,Proc.P.my_num()==0,Proc.MCp.get_alphai());
        Proc.write_Sp(r[0],Sansp);
	#else
        Proc.get_Sp_ref(r[0]).add(Proc.read_Sp(r[1]),Proc.temp.ansp,Proc.P.my_num()==0,Proc.MCp.get_alphai());
	#endif
#endif
        break;
      case E_MP_ADDSI:
        Proc.temp.ansp.assign(mp_n);
        Proc.MPAddm_Ext(Proc.get_Sp_ref(r[1]), Proc.temp.ansp, Proc.get_Sp_ref(r[0]));
        break;
      case GADDSI:
        Proc.temp.ans2.assign(n);
#if defined(EXTENDED_SPDZ_GF2N)
        Proc.GAddm_Ext(Proc.get_S2_ref(r[1]), Proc.temp.ans2, Proc.get_S2_ref(r[0]));
#else
	#ifdef DEBUG
           Sans2.add(Proc.read_S2(r[1]),Proc.temp.ans2,Proc.P.my_num()==0,Proc.MC2.get_alphai());
	   Proc.write_S2(r[0],Sans2);
        #else
           Proc.get_S2_ref(r[0]).add(Proc.read_S2(r[1]),Proc.temp.ans2,Proc.P.my_num()==0,Proc.MC2.get_alphai());
	#endif
#endif
        break;
      case SUBCI:
        Proc.temp.ansp.assign(n);
	#ifdef DEBUG
           Proc.temp.ansp.sub(Proc.read_Cp(r[1]),Proc.temp.ansp);
           Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).sub(Proc.read_Cp(r[1]),Proc.temp.ansp);
	#endif
        break;
      case GSUBCI:
        Proc.temp.ans2.assign(n);
	#ifdef DEBUG
           Proc.temp.ans2.sub(Proc.read_C2(r[1]),Proc.temp.ans2);
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).sub(Proc.read_C2(r[1]),Proc.temp.ans2);
	#endif
        break;
      case SUBSI:
    	  Proc.temp.ansp.assign(n);
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.PSubml_Ext(Proc.get_Sp_ref(r[1]), Proc.temp.ansp, Proc.get_Sp_ref(r[0]));
#else
	#ifdef DEBUG
    	  Sansp.sub(Proc.read_Sp(r[1]),Proc.temp.ansp,Proc.P.my_num()==0,Proc.MCp.get_alphai());
    	  Proc.write_Sp(r[0],Sansp);
	#else
    	  Proc.get_Sp_ref(r[0]).sub(Proc.read_Sp(r[1]),Proc.temp.ansp,Proc.P.my_num()==0,Proc.MCp.get_alphai());
	#endif
#endif
        break;
      case E_MP_SUBSIL:
         Proc.temp.ansp.assign(mp_n);
         Proc.MPSubml_Ext(Proc.get_Sp_ref(r[1]), Proc.temp.ansp, Proc.get_Sp_ref(r[0]));
         break;
      case E_MP_SUBSIR:
    	  Proc.temp.ansp.assign(mp_n);
         Proc.MPSubmr_Ext(Proc.temp.ansp, Proc.get_Sp_ref(r[1]), Proc.get_Sp_ref(r[0]));
        break;
      case GSUBSI:
        Proc.temp.ans2.assign(n);
#if defined(EXTENDED_SPDZ_GFP)
    	  Proc.GSubml_Ext(Proc.get_S2_ref(r[1]), Proc.temp.ans2, Proc.get_S2_ref(r[0]));
#else
  	#ifdef DEBUG
           Sans2.sub(Proc.read_S2(r[1]),Proc.temp.ans2,Proc.P.my_num()==0,Proc.MC2.get_alphai());
	   Proc.write_S2(r[0],Sans2);
        #else
           Proc.get_S2_ref(r[0]).sub(Proc.read_S2(r[1]),Proc.temp.ans2,Proc.P.my_num()==0,Proc.MC2.get_alphai());
        #endif
#endif
        break;
      case SUBCFI:
        Proc.temp.ansp.assign(n);
	#ifdef DEBUG
           Proc.temp.ansp.sub(Proc.temp.ansp,Proc.read_Cp(r[1]));
           Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).sub(Proc.temp.ansp,Proc.read_Cp(r[1]));
	#endif
        break;
      case GSUBCFI:
        Proc.temp.ans2.assign(n);
	#ifdef DEBUG
           Proc.temp.ans2.sub(Proc.temp.ans2,Proc.read_C2(r[1]));
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).sub(Proc.temp.ans2,Proc.read_C2(r[1]));
	#endif
        break;
      case SUBSFI:
        Proc.temp.ansp.assign(n);
#if defined(EXTENDED_SPDZ_GFP)
        Proc.PSubmr_Ext(Proc.temp.ansp, Proc.get_Sp_ref(r[1]), Proc.get_Sp_ref(r[0]));
#else
	#ifdef DEBUG
        Sansp.sub(Proc.temp.ansp,Proc.read_Sp(r[1]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
        Proc.write_Sp(r[0],Sansp);
	#else
        Proc.get_Sp_ref(r[0]).sub(Proc.temp.ansp,Proc.read_Sp(r[1]),Proc.P.my_num()==0,Proc.MCp.get_alphai());
	#endif
#endif
        break;
      case GSUBSFI:
    	  Proc.temp.ans2.assign(n);
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GSubmr_Ext(Proc.temp.ans2, Proc.get_S2_ref(r[1]), Proc.get_S2_ref(r[0]));
#else
	#ifdef DEBUG
    	  Sans2.sub(Proc.temp.ans2,Proc.read_S2(r[1]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
    	  Proc.write_S2(r[0],Sans2);
	#else
    	  Proc.get_S2_ref(r[0]).sub(Proc.temp.ans2,Proc.read_S2(r[1]),Proc.P.my_num()==0,Proc.MC2.get_alphai());
	#endif
#endif
        break;
      case MULCI:
        Proc.temp.ansp.assign(n);
	#ifdef DEBUG
           Proc.temp.ansp.mul(Proc.temp.ansp,Proc.read_Cp(r[1]));
           Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).mul(Proc.temp.ansp,Proc.read_Cp(r[1]));
	#endif
        break;
      case GMULCI:
        Proc.temp.ans2.assign(n);
	#ifdef DEBUG
           Proc.temp.ans2.mul(Proc.temp.ans2,Proc.read_C2(r[1]));
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).mul(Proc.temp.ans2,Proc.read_C2(r[1]));
	#endif
        break;
      case MULSI:
        Proc.temp.ansp.assign(n);
#if defined(EXTENDED_SPDZ_GFP)
        Proc.PMulm_Ext(Proc.get_Sp_ref(r[0]), Proc.read_Sp(r[1]), Proc.temp.ansp);
#else
	#ifdef DEBUG
        Sansp.mul(Proc.read_Sp(r[1]),Proc.temp.ansp);
        Proc.write_Sp(r[0],Sansp);
	#else
        Proc.get_Sp_ref(r[0]).mul(Proc.read_Sp(r[1]),Proc.temp.ansp);
	#endif
#endif
        break;
      case E_MP_MULSI:
         Proc.temp.ansp.assign(mp_n);
         Proc.MPMulm_Ext(Proc.get_Sp_ref(r[0]), Proc.read_Sp(r[1]), Proc.temp.ansp);
        break;
      case GMULSI:
    	  Proc.temp.ans2.assign(n);
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GMulm_Ext(Proc.get_S2_ref(r[0]), Proc.read_S2(r[1]), Proc.temp.ans2);
#else
	#ifdef DEBUG
    	  Sans2.mul(Proc.read_S2(r[1]),Proc.temp.ans2);
    	  Proc.write_S2(r[0],Sans2);
	#else
    	  Proc.get_S2_ref(r[0]).mul(Proc.read_S2(r[1]),Proc.temp.ans2);
	#endif
#endif
        break;
      case TRIPLE:
#if defined(EXTENDED_SPDZ_GFP)
    	Proc.PTriple_Ext(Proc.get_Sp_ref(r[0]),Proc.get_Sp_ref(r[1]),Proc.get_Sp_ref(r[2]));
#else
    	Proc.DataF.get_three(DATA_MODP, DATA_TRIPLE, Proc.get_Sp_ref(r[0]),Proc.get_Sp_ref(r[1]),Proc.get_Sp_ref(r[2]));
#endif
        break;
      case GTRIPLE:
#if defined(EXTENDED_SPDZ_GF2N)
    	Proc.GTriple_Ext(Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]),Proc.get_S2_ref(r[2]));
#else
        Proc.DataF.get_three(DATA_GF2N, DATA_TRIPLE, Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]),Proc.get_S2_ref(r[2]));
#endif
        break;
      case GBITTRIPLE:
        Proc.DataF.get_three(DATA_GF2N, DATA_BITTRIPLE, Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]),Proc.get_S2_ref(r[2]));
        break;
      case GBITGF2NTRIPLE:
        Proc.DataF.get_three(DATA_GF2N, DATA_BITGF2NTRIPLE, Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]),Proc.get_S2_ref(r[2]));
        break;
      case SQUARE:
        Proc.DataF.get_two(DATA_MODP, DATA_SQUARE, Proc.get_Sp_ref(r[0]),Proc.get_Sp_ref(r[1]));
        break;
      case GSQUARE:
        Proc.DataF.get_two(DATA_GF2N, DATA_SQUARE, Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]));
        break;
      case BIT:
#if defined(EXTENDED_SPDZ_GFP)
    	Proc.PBit_Ext(Proc.get_Sp_ref(r[0]));
#else
        Proc.DataF.get_one(DATA_MODP, DATA_BIT, Proc.get_Sp_ref(r[0]));
#endif
        break;
      case GBIT:
#if defined(EXTENDED_SPDZ_GF2N)
    	Proc.GBit_Ext(Proc.get_S2_ref(r[0]));
#else
        Proc.DataF.get_one(DATA_GF2N, DATA_BIT, Proc.get_S2_ref(r[0]));
#endif
        break;
      case INV:
#if defined(EXTENDED_SPDZ_GFP)
    	Proc.PInverse_Ext(Proc.get_Sp_ref(r[0]),Proc.get_Sp_ref(r[1]));
#else
        Proc.DataF.get_two(DATA_MODP, DATA_INVERSE, Proc.get_Sp_ref(r[0]),Proc.get_Sp_ref(r[1]));
#endif
        break;
      case GINV:
#if defined(EXTENDED_SPDZ_GF2N)
    	Proc.GInverse_Ext(Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]));
#else
        Proc.DataF.get_two(DATA_GF2N, DATA_INVERSE, Proc.get_S2_ref(r[0]),Proc.get_S2_ref(r[1]));
#endif
        break;
      case INPUTMASK:
        Proc.DataF.get_input(Proc.get_Sp_ref(r[0]), Proc.temp.ansp, n);
        if (n == Proc.P.my_num())
          Proc.temp.ansp.output(Proc.private_output, false);
        break;
      case GINPUTMASK:
        Proc.DataF.get_input(Proc.get_S2_ref(r[0]), Proc.temp.ans2, n);
        if (n == Proc.P.my_num())
          Proc.temp.ans2.output(Proc.private_output, false);
        break;
      case INPUT:
#if defined(EXTENDED_SPDZ_GFP)
        Proc.PInput_Ext(Proc.get_Sp_ref(r[0]), n);
#else
        { gfp& rr=Proc.temp.rrp; gfp& t=Proc.temp.tp; gfp& tmp=Proc.temp.tmpp;
          Proc.DataF.get_input(Proc.get_Sp_ref(r[0]),rr,n);
          octetStream o;
          if (n==Proc.P.my_num())
            { gfp& xi=Proc.temp.xip;
	      #ifdef DEBUG
	         printf("Enter your input : \n");
	      #endif
              word x;
              cin >> x;
              t.assign(x);
              t.sub(t,rr);
              t.pack(o);
              Proc.P.send_all(o);
              xi.add(t,Proc.get_Sp_ref(r[0]).get_share());
              Proc.get_Sp_ref(r[0]).set_share(xi);
            }
          else
            { Proc.P.receive_player(n,o);
              t.unpack(o);
            }
          tmp.mul(Proc.MCp.get_alphai(),t);
          tmp.add(Proc.get_Sp_ref(r[0]).get_mac(),tmp);
          Proc.get_Sp_ref(r[0]).set_mac(tmp);
        }
#endif
        break;
      case GINPUT:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GInput_Ext(Proc.get_S2_ref(r[0]), n);
#else
        { gf2n& rr=Proc.temp.rr2; gf2n& t=Proc.temp.t2; gf2n& tmp=Proc.temp.tmp2;
          Proc.DataF.get_input(Proc.get_S2_ref(r[0]),rr,n);
          octetStream o;
          if (n==Proc.P.my_num())
            { gf2n& xi=Proc.temp.xi2;
	      #ifdef DEBUG
	         printf("Enter your input : \n");
              #endif
              word x;
              cin >> x;
              t.assign(x);
              t.sub(t,rr);
              t.pack(o);
              Proc.P.send_all(o);
              xi.add(t,Proc.get_S2_ref(r[0]).get_share());
              Proc.get_S2_ref(r[0]).set_share(xi);
            }
          else
            { Proc.P.receive_player(n,o);
              t.unpack(o);
            }
          tmp.mul(Proc.MC2.get_alphai(),t);
          tmp.add(Proc.get_S2_ref(r[0]).get_mac(),tmp);
          Proc.get_S2_ref(r[0]).set_mac(tmp);
        }
#endif
        break;
      case STARTINPUT:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "STARTINPUT instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
    	  Proc.inputp.start(r[0],n);
#endif
        break;
      case GSTARTINPUT:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "GSTARTINPUT instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
        Proc.input2.start(r[0],n);
#endif
        break;
      case STOPINPUT:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "STOPINPUT instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
    	  Proc.inputp.stop(n,start);
#endif
        break;
      case GSTOPINPUT:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "GSTOPINPUT instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
        Proc.input2.stop(n,start);
#endif
        break;
      case ANDC:
	#ifdef DEBUG
           ansp.AND(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	   Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).AND(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	#endif
        break;
      case GANDC:
	#ifdef DEBUG
           Proc.temp.ans2.AND(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	   Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).AND(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	#endif
        break;
      case XORC:
	#ifdef DEBUG
           Proc.temp.ansp.XOR(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	   Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).XOR(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	#endif
        break;
      case GXORC:
	#ifdef DEBUG
           Proc.temp.ans2.XOR(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	   Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).XOR(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	#endif
        break;
      case ORC:
	#ifdef DEBUG
           Proc.temp.ansp.OR(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	   Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).OR(Proc.read_Cp(r[1]),Proc.read_Cp(r[2]));
	#endif
        break;
      case GORC:
	#ifdef DEBUG
           Proc.temp.ans2.OR(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	   Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).OR(Proc.read_C2(r[1]),Proc.read_C2(r[2]));
	#endif
        break;
      case ANDCI:
    	 Proc.temp.ansp.assign(n);
    	 Proc.get_Cp_ref(r[0]).AND(Proc.read_Cp(r[1]),Proc.temp.ansp);
        break;
      case GANDCI:
        Proc.temp.ans2.assign(n);
        Proc.get_C2_ref(r[0]).AND(Proc.temp.ans2,Proc.read_C2(r[1]));
        break;
      case XORCI:
    	 Proc.temp.ansp.assign(n);
    	 Proc.get_Cp_ref(r[0]).XOR(Proc.read_Cp(r[1]),Proc.temp.ansp);
        break;
      case GXORCI:
        Proc.temp.ans2.assign(n);
	#ifdef DEBUG
           Proc.temp.ans2.XOR(Proc.temp.ans2,Proc.read_C2(r[1]));
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).XOR(Proc.temp.ans2,Proc.read_C2(r[1]));
	#endif
        break;
      case ORCI:
    	  Proc.temp.ansp.assign(n);
    	  Proc.get_Cp_ref(r[0]).OR(Proc.read_Cp(r[1]),Proc.temp.ansp);
        break;
      case GORCI:
        Proc.temp.ans2.assign(n);
	#ifdef DEBUG
           Proc.temp.ans2.OR(Proc.temp.ans2,Proc.read_C2(r[1]));
           Proc.write_C2(r[0],Proc.temp.ans2);
	#else
	   Proc.get_C2_ref(r[0]).OR(Proc.temp.ans2,Proc.read_C2(r[1]));
	#endif
        break;
      // Note: Fp version has different semantics for NOTC than GNOTC
      case NOTC:
    	 exit(EBADR);
        break;
      case GNOTC:
	#ifdef DEBUG
           Proc.temp.ans2.NOT(Proc.read_C2(r[1]));
	   Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).NOT(Proc.read_C2(r[1]));
	#endif
        break;
      case SHLC:
    	 exit(EBADR);
        break;
      case SHRC:
     	 exit(EBADR);
        break;
      case SHLCI:
	#ifdef DEBUG
           Proc.temp.ansp.SHL(Proc.read_Cp(r[1]),Proc.temp.aa);
	   Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).SHL(Proc.read_Cp(r[1]),n);
	#endif
        break;
      case GSHLCI:
	#ifdef DEBUG
           Proc.temp.ans2.SHL(Proc.read_C2(r[1]),n);
	   Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).SHL(Proc.read_C2(r[1]),n);
	#endif
        break;
      case SHRCI:
	#ifdef DEBUG
           Proc.temp.ansp.SHR(Proc.read_Cp(r[1]),Proc.temp.aa);
	   Proc.write_Cp(r[0],Proc.temp.ansp);
	#else
           Proc.get_Cp_ref(r[0]).SHR(Proc.read_Cp(r[1]),n);
	#endif
        break;
      case GSHRCI:
	#ifdef DEBUG
           Proc.temp.ans2.SHR(Proc.read_C2(r[1]),n);
	   Proc.write_C2(r[0],Proc.temp.ans2);
	#else
           Proc.get_C2_ref(r[0]).SHR(Proc.read_C2(r[1]),n);
	#endif
        break;
      case GBITDEC:
        for (int j = 0; j < size; j++)
          {
            gf2n::internal_type a = Proc.read_C2(r[0] + j).get();
            for (unsigned int i = 0; i < start.size(); i++)
              {
                Proc.get_C2_ref(start[i] + j) = a & 1;
                a >>= n;
              }
          }
        return;
      case GBITCOM:
        for (int j = 0; j < size; j++)
          {
            gf2n::internal_type a = 0;
            for (unsigned int i = 0; i < start.size(); i++)
              {
                a ^= Proc.read_C2(start[i] + j).get() << (i * n);
              }
            Proc.get_C2_ref(r[0] + j) = a;
          }
        return;
      case BITDECINT:
        {
          long a = Proc.read_Ci(r[0]);
          for (unsigned int i = 0; i < start.size(); i++)
            {
              Proc.get_Ci_ref(start[i]) = (a >> i) & 1;
            }
          break;
        }
      case STARTOPEN:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "STARTOPEN instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
    	  Proc.POpen_Start(start,Proc.P,Proc.MCp,size);
#endif
        return;
      case GSTARTOPEN:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "GSTARTOPEN instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
        Proc.POpen_Start(start,Proc.P,Proc.MC2,size);
#endif
        return;
      case STOPOPEN:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "STOPOPEN instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
    	  Proc.POpen_Stop(start,Proc.P,Proc.MCp,size);
#endif
        return;
      case GSTOPOPEN:
#if defined(EXTENDED_SPDZ_GFP)
    	  std::cerr << "GSTOPOPEN instruction is not implemented in extension" << std::endl;
    	  exit(-1);
#else
        Proc.POpen_Stop(start,Proc.P,Proc.MC2,size);
#endif
        return;
      case OPEN:
#if defined(EXTENDED_SPDZ_GFP)
          Proc.POpen_Ext(start, size);
#else
          Proc.POpen(start, Proc.P, Proc.MCp, size);
#endif
        break;
      case GOPEN:
#if defined(EXTENDED_SPDZ_GF2N)
    	  Proc.GOpen_Ext(start, size);
#else
          Proc.POpen(start, Proc.P, Proc.MC2, size);
#endif
        break;
      case E_MP_OPEN:
    	  Proc.MPOpen_Ext(start, size);
    	  break;
#if defined(EXTENDED_SPDZ_GFP)
      case E_STARTMULT:
    	  std::cerr << "E_STARTMULT instruction is not implemented in extension" << std::endl;
    	  exit(-1);
        return;
      case E_STOPMULT:
    	  std::cerr << "E_STOPMULT instruction is not implemented in extension" << std::endl;
    	  exit(-1);
        return;
      case E_MULT:
        Proc.PMult_Ext(start, size);
        break;
      case E_MP_MULT:
        Proc.MPMult_Ext(start, size);
        break;
#endif
#if defined(EXTENDED_SPDZ_GF2N)
      case GE_MULT:
        Proc.GMult_Ext(start, size);
        break;
#endif
#if defined(EXTENDED_SPDZ_Z2N)
      case E_SKEW_BIT_DEC:
    	  Proc.Skew_Bit_Decomp_Ext(start, Proc.get_Sp_ref(r[0]), size);
    	  break;
      case E_MP_SKEW_BIT_DEC:
    	  Proc.MP_Skew_Bit_Decomp_Ext(start, Proc.get_Sp_ref(r[0]), size);
    	  break;
      case E_SKEW_BIT_REC:
		  Proc.Skew_Bit_Recomp_Ext(start, Proc.get_S2_ref(r[0]), size);
    	  break;
      case E_SKEW_RING_REC:
		  Proc.Skew_Ring_Recomp_Ext(Proc.get_Sp_ref(r[0]), start, size);
    	  break;
      case E_MP_SKEW_RING_REC:
		  Proc.MP_Skew_Ring_Recomp_Ext(Proc.get_Sp_ref(r[0]), start, size);
    	  break;
      case E_SKEW_BIT_INJ:
		  Proc.Skew_Bit_Inject_Ext(start, Proc.get_S2_ref(r[0]), size);
    	  break;
      case E_MP_SKEW_BIT_INJ:
      	  Proc.MP_Skew_Bit_Inject_Ext(start, Proc.get_S2_ref(r[0]), size);
         break;
#endif
      case JMP:
        Proc.PC += (signed int) n;
        break;
      case JMPI:
        Proc.PC += (signed int) Proc.read_Ci(r[0]);
        break;
      case JMPNZ:
        if (Proc.read_Ci(r[0]) != 0)
          { Proc.PC += (signed int) n; }
        break;
      case JMPEQZ:
        if (Proc.read_Ci(r[0]) == 0)
          { Proc.PC += (signed int) n; }
        break;
      case EQZC:
        if (Proc.read_Ci(r[1]) == 0)
          Proc.write_Ci(r[0], 1);
        else
          Proc.write_Ci(r[0], 0);
        break;
      case LTZC:
        if (Proc.read_Ci(r[1]) < 0)
          Proc.write_Ci(r[0], 1);
        else
          Proc.write_Ci(r[0], 0);
       break;
      case LTC:
        if (Proc.read_Ci(r[1]) < Proc.read_Ci(r[2]))
          Proc.write_Ci(r[0], 1);
        else
          Proc.write_Ci(r[0], 0);
        break;
      case GTC:
        if (Proc.read_Ci(r[1]) > Proc.read_Ci(r[2]))
          Proc.write_Ci(r[0], 1);
        else
          Proc.write_Ci(r[0], 0);
        break;
      case EQC:
        if (Proc.read_Ci(r[1]) == Proc.read_Ci(r[2]))
          Proc.write_Ci(r[0], 1);
        else
          Proc.write_Ci(r[0], 0);
        break;
      case LDINT:
        Proc.write_Ci(r[0], n);
        break;
      case ADDINT:
        Proc.get_Ci_ref(r[0]) = Proc.read_Ci(r[1]) + Proc.read_Ci(r[2]);
        break;
      case SUBINT:
        Proc.get_Ci_ref(r[0]) = Proc.read_Ci(r[1]) - Proc.read_Ci(r[2]);
        break;
      case MULINT:
        Proc.get_Ci_ref(r[0]) = Proc.read_Ci(r[1]) * Proc.read_Ci(r[2]);
        break;
      case DIVINT:
        Proc.get_Ci_ref(r[0]) = Proc.read_Ci(r[1]) / Proc.read_Ci(r[2]);
        break;
      case CONVINT:
        Proc.get_Cp_ref(r[0]).assign(Proc.read_Ci(r[1]));
        break;
      case GCONVINT:
        Proc.get_C2_ref(r[0]).assign((word)Proc.read_Ci(r[1]));
        break;
      case CONVMODP:
        Proc.write_Ci(r[0], (long)(Proc.read_Cp(r[1]).get()));
        break;
      case GCONVGF2N:
        Proc.write_Ci(r[0], Proc.read_C2(r[1]).get_word());
        break;
      case PRINTMEM:
        if (Proc.P.my_num() == 0)
	  { cout << "Mem[" <<  r[0] << "] = " << Proc.machine.Mp.read_C(r[0]) << endl; }
        break;
      case GPRINTMEM:
        if (Proc.P.my_num() == 0)
	  { cout << "Mem[" <<  r[0] << "] = " << Proc.machine.M2.read_C(r[0]) << endl; }
        break;
      case PRINTREG:
        if (Proc.P.my_num() == 0)
           {
             cout << "Reg[" << r[0] << "] = " << Proc.read_Cp(r[0])
              << " # " << string((char*)&n,sizeof(n)) << endl;
           }
        break;
      case GPRINTREG:
        if (Proc.P.my_num() == 0)
           {
             cout << "Reg[" << r[0] << "] = " << Proc.read_C2(r[0])
              << " # " << string((char*)&n,sizeof(n)) << endl;
           }
        break;
      case PRINTREGPLAIN:
        if (Proc.P.my_num() == 0)
           {
             cout << Proc.read_Cp(r[0]) << flush;
           }
        break;
      case GPRINTREGPLAIN:
        if (Proc.P.my_num() == 0)
           {
             cout << Proc.read_C2(r[0]) << flush;
           }
        break;
      case E_PRINTFIXEDPLAIN:
    	  if (Proc.P.my_num() == 0)
    	  {
    		  uint64_t tmp_val = Proc.read_Cp(r[0]).get();
    		  uint64_t val;
    		  double res =0;
    		  double int_part =0;
    		  double dec_part =0;
    		  int sign = (tmp_val >> 63) & 1;
    		  if (sign == 1)
    		  {
    			  val = 0 - tmp_val;
    		  }
    		  else
    		  {
    			  val = tmp_val;
    		  }
    		  for(int i = 0; i < 64; i++)
    		  {
    			  if(i < n)
    			  {
    				  dec_part += ((val >> i) & 1) * pow(2.0, (double) (-(n-i)));
    			  }
    			  else
    			  {
    				  int_part += ((val >> i) & 1) * pow(2.0, (double) (i-n));
    			  }
    		  }
    		  res = int_part + dec_part;
    		  if (sign == 1)
    		  {
    			  cout << "-" << fixed << res << flush;;
    		  }
    		  else
    		  {
    			  cout << fixed << res << flush;;
    		  }
    	  }
    	  break;
      case PRINTFLOATPLAIN:
        if (Proc.P.my_num() == 0)

        {
#ifdef DEBUG
        	int flag = 0;
        	//flag = 1: binary, flag = 0 : decimal, flag = 2: hexadecimal
        	uint64_t v = Proc.read_Cp(start[0]).get();
        	uint64_t p = Proc.read_Cp(start[1]).get();
        	uint64_t z = Proc.read_Cp(start[2]).get();
        	uint64_t s = Proc.read_Cp(start[3]).get();
        	int p2 = p;
        	std::string res;
        	//std::string str_v2;
        	//cout << p << flush;
        	//std::string str_v = std::to_string(v2);
        	//std::string res2;
        	std::string bin_v = std::bitset<24>(v).to_string();
        	std::stringstream hex_v;
        	hex_v << std::hex << v;
        	if(flag == 0){
        		if (z == 0){
        			//if zero flag is 1
        	       if(s == 1){
        	    	   res = "-" + std::to_string(v);
        	    	   res = res + " * 2^(" + std::to_string(p2) + ")";
        	       }
        	       else{
        	    	   //for Debug
        	       	//p = 1 - p;
        	       	res = std::to_string(v) + " * 2^(" + std::to_string(p2) + ")";
        	       	}
        	       cout << res << flush;
        		}
        		else{
        			cout << 1-z <<flush;
        		}
        	}
        	if(flag == 1){
        		if (z == 0){
        			//if zero flag is 1
        			if(s == 1){
        				res = "-" + bin_v;
        				res = res + " * 2^(" + std::to_string(p2) + ")";
        			}
        			else{
        				res = bin_v + " * 2^(" + std::to_string(p2) + ")";
        			}
        			cout << res << flush;
        		}
        		else{
        			cout << 1-z <<flush;
        		}
        	}
        	if(flag == 2){
        		if (z == 0){
        			//if zero flag is 1
        			if(s == 1){
        				res = "- 0x" + hex_v.str();
        				res = res + " * 2^(" + std::to_string(p2) + ")";
        			}
        			else{
        				res = "0x" + hex_v.str() + " * 2^(" + std::to_string(p2) + ")";
        			}
        			cout << res << flush;
        		}
        		else{
        			cout << 1-z <<flush;
        		}
        	}
#else
        	uint64_t v = Proc.read_Cp(start[0]).get();
        	int64_t p = (int64_t) Proc.read_Cp(start[1]).get();
        	uint64_t z = Proc.read_Cp(start[2]).get();
        	uint64_t s = Proc.read_Cp(start[3]).get();
        	double p2 = 1.0 * p;
        	double v_val = 1.0 * v;
        	double p_val = pow(2.0, p2);
        	if (z == 0){
        		//if zero flag is 1
        		if(s == 1){
        			v_val = (-1.0) * v_val * p_val;
        		}
        		else{
        			//for Debug
        			//p = 1 - p;
        			v_val = v_val * p_val;
        		}
        		cout << std::scientific << v_val << flush;
        	}
        	else{
        		cout << 1-z <<flush;
        	}
#endif
        }
      break;
      case PRINTINT:
        if (Proc.P.my_num() == 0)
           {
             cout << Proc.read_Ci(r[0]) << flush;
           }
        break;
      case PRINTSTR:
        if (Proc.P.my_num() == 0)
           {
             cout << string((char*)&n,sizeof(n)) << flush;
           }
        break;
      case PRINTCHR:
        if (Proc.P.my_num() == 0)
           {
             cout << string((char*)&n,1) << flush;
           }
        break;
      case PRINTCHRINT:
        if (Proc.P.my_num() == 0)
           {
             cout << string((char*)&(Proc.read_Ci(r[0])),1) << flush;
           }
        break;
      case PRINTSTRINT:
        if (Proc.P.my_num() == 0)
           {
             cout << string((char*)&(Proc.read_Ci(r[0])),sizeof(int)) << flush;
           }
        break;
      case RAND:
        Proc.write_Ci(r[0], Proc.prng.get_uint() % (1 << Proc.read_Ci(r[1])));
        break;
      case REQBL:
      case GREQBL:
      case USE:
      case USE_INP:
      case USE_PREP:
      case GUSE_PREP:
        break;
      case TIME:
        Proc.machine.time();
	break;
      case START:
        Proc.machine.start(n);
        break;
      case STOP:
        Proc.machine.stop(n);
        break;
      case RUN_TAPE:
        Proc.DataF.skip(Proc.machine.run_tape(r[0], n, r[1], -1));
        break;
      case JOIN_TAPE:
        Proc.machine.join_tape(r[0]);
        break;
      case CRASH:
        throw crash_requested();
        break;
      case STARTGRIND:
        CALLGRIND_START_INSTRUMENTATION;
        break;
      case STOPGRIND:
        CALLGRIND_STOP_INSTRUMENTATION;
        break;
      // ***
      // TODO: read/write shared GF(2^n) data instructions
      // ***
      case LISTEN:
        // listen for connections at port number n
        Proc.external_clients.start_listening(n);
        break;
      case ACCEPTCLIENTCONNECTION:
      {
        // get client connection at port number n + my_num())
        int client_handle = Proc.external_clients.get_client_connection(n);
        if (client_handle == -1)
        {
          stringstream ss;
          ss << "No connection on port " << r[0] << endl;
          throw Processor_Error(ss.str());
        }
        Proc.write_Ci(r[0], client_handle);
        break;
      }
      case CONNECTIPV4:
      {
        // connect to server at port n + my_num()
        int ipv4 = Proc.read_Ci(r[1]);
        int server_handle = Proc.external_clients.connect_to_server(n, ipv4);
        Proc.write_Ci(r[0], server_handle);
        break;
      }
      case READCLIENTPUBLICKEY:
        Proc.read_client_public_key(Proc.read_Ci(r[0]), start);
        break;
      case INITSECURESOCKET:
    	  abort();
    	  break;
      case RESPSECURESOCKET:
    	  abort();
    	  break;
      case READSOCKETINT:
        Proc.read_socket_ints(Proc.read_Ci(r[0]), start);
        break;
      case READSOCKETC:
        Proc.read_socket_vector<gfp>(Proc.read_Ci(r[0]), start);
        break;
      case READSOCKETS:
        // read shares and MAC shares
        Proc.read_socket_private<gfp>(Proc.read_Ci(r[0]), start, true);
        break;
      case GREADSOCKETS:
        //Proc.get_S2_ref(r[0]).get_share().pack(socket_octetstream);
        //Proc.get_S2_ref(r[0]).get_mac().pack(socket_octetstream);
        break;
      case WRITESOCKETINT:
        Proc.write_socket(INT, CLEAR, false, Proc.read_Ci(r[0]), r[1], start);
        break;
      case WRITESOCKETC:
        Proc.write_socket(MODP, CLEAR, false, Proc.read_Ci(r[0]), r[1], start);
        break;
      case WRITESOCKETS:
        // Send shares + MACs
        Proc.write_socket(MODP, SECRET, true, Proc.read_Ci(r[0]), r[1], start);
        break;
      case WRITESOCKETSHARE:
        // Send only shares, no MACs
        // N.B. doesn't make sense to have a corresponding read instruction for this
        Proc.write_socket(MODP, SECRET, false, Proc.read_Ci(r[0]), r[1], start);
        break;
      /*case GWRITESOCKETS:
        Proc.get_S2_ref(r[0]).get_share().pack(socket_octetstream);
        Proc.get_S2_ref(r[0]).get_mac().pack(socket_octetstream);
        break;*/
      case WRITEFILESHARE:
        // Write shares to file system
        Proc.write_shares_to_file<gfp>(start);
        break;
      case READFILESHARE:
        // Read shares from file system
        Proc.read_shares_from_file<gfp>(Proc.read_Ci(r[0]), r[1], start);
        break;        
      case PUBINPUT:
        Proc.public_input >> Proc.get_Ci_ref(r[0]);
        break;
      case RAWOUTPUT:
        Proc.read_Cp(r[0]).output(Proc.public_output, false);
        break;
      case GRAWOUTPUT:
        Proc.read_C2(r[0]).output(Proc.public_output, false);
        break;
      case STARTPRIVATEOUTPUT:
        Proc.privateOutputp.start(n,r[0],r[1]);
        break;
      case GSTARTPRIVATEOUTPUT:
        Proc.privateOutput2.start(n,r[0],r[1]);
        break;
      case STOPPRIVATEOUTPUT:
        Proc.privateOutputp.stop(n,r[0]);
        break;
      case GSTOPPRIVATEOUTPUT:
        Proc.privateOutput2.stop(n,r[0]);
        break;
      case PREP:
        Proc.DataF.get<gfp>(Proc, r, start, size);
        return;
      case GPREP:
        Proc.DataF.get<gf2n>(Proc, r, start, size);
        return;
      default:
        printf("Case of opcode=%d not implemented yet\n",opcode);
        throw not_implemented();
        break;
    }
  if (size > 1)
    {
      r[0]++; r[1]++; r[2]++;
    }
  }
}
