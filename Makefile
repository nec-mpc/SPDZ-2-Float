# (C) 2018 University of Bristol, Bar-Ilan University. See License.txt


include CONFIG

MATH = $(patsubst %.cpp,%.o,$(wildcard Math/*.cpp))

TOOLS = $(patsubst %.cpp,%.o,$(wildcard Tools/*.cpp))

NETWORK = $(patsubst %.cpp,%.o,$(wildcard Networking/*.cpp))

AUTH = $(patsubst %.cpp,%.o,$(wildcard Auth/*.cpp))

PROCESSOR = $(patsubst %.cpp,%.o,$(wildcard Processor/*.cpp))

GC = $(patsubst %.cpp,%.o,$(wildcard GC/*.cpp))

# OT stuff needs GF2N_LONG, so only compile if this is enabled
#ifeq ($(USE_GF2N_LONG),1)
#OT = $(patsubst %.cpp,%.o,$(filter-out OT/OText_main.cpp,$(wildcard OT/*.cpp)))
#OT_EXE = ot.x ot-offline.x
#endif

COMMON = $(MATH) $(TOOLS) $(NETWORK) $(AUTH)
COMPLETE = $(COMMON) $(PROCESSOR) $(TINYOTOFFLINE) $(OT)
BMR = $(patsubst %.cpp,%.o,$(wildcard BMR/*.cpp BMR/network/*.cpp)) $(COMMON) $(PROCESSOR) $(GC)


LIB = libSPDZ.a
LIBSIMPLEOT = SimpleOT/libsimpleot.a

# used for dependency generation
OBJS = $(BMR) $(TINYOTOFFLINE)
DEPS := $(OBJS:.o=.d)


all: gen_input online 

ifeq ($(USE_NTL),1)
all: overdrive 
endif

-include $(DEPS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -MMD -c -o $@ $<

online: Server.x Player-Online.x

#offline: $(OT_EXE)

gen_input: gen_input_f2n.x gen_input_fp.x

bmr: bmr-program-party.x bmr-program-tparty.x

overdrive: simple-offline.x pairwise-offline.x cnc-offline.x

Server.x: Server.cpp $(COMMON)
	$(CXX) $(CFLAGS) Server.cpp -o Server.x $(COMMON) $(LDLIBS)

Player-Online.x: Player-Online.cpp $(COMMON) $(PROCESSOR)
	$(CXX) $(CFLAGS) Player-Online.cpp -o Player-Online.x $(COMMON) $(PROCESSOR) $(LDLIBS)

ifeq ($(USE_GF2N_LONG),1)
ot.x: $(OT) $(COMMON) OT/OText_main.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(LIBSIMPLEOT)

ot-check.x: $(OT) $(COMMON)
	$(CXX) $(CFLAGS) -o ot-check.x OT/BitVector.o OT/OutputCheck.cpp $(COMMON) $(LDLIBS)

ot-bitmatrix.x: $(OT) $(COMMON) OT/BitMatrixTest.cpp
	$(CXX) $(CFLAGS) -o ot-bitmatrix.x OT/BitMatrixTest.cpp OT/BitMatrix.o OT/BitVector.o $(COMMON) $(LDLIBS)

ot-offline.x: $(OT) $(COMMON) ot-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(LIBSIMPLEOT)
endif

check-passive.x: $(COMMON) check-passive.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

gen_input_f2n.x: Scripts/gen_input_f2n.cpp $(COMMON)
	$(CXX) $(CFLAGS) Scripts/gen_input_f2n.cpp	-o gen_input_f2n.x $(COMMON) $(LDLIBS)

gen_input_fp.x: Scripts/gen_input_fp.cpp $(COMMON)
	$(CXX) $(CFLAGS) Scripts/gen_input_fp.cpp	-o gen_input_fp.x $(COMMON) $(LDLIBS)



gc-emulate.x: $(GC) $(COMMON) $(PROCESSOR) gc-emulate.cpp $(BMR)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(BOOST)

bmr-program-party.x: $(BMR) bmr-program-party.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(BOOST)

bmr-program-tparty.x: $(BMR) bmr-program-tparty.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS) $(BOOST)

bmr-clean:
	-rm BMR/*.o BMR/*/*.o GC/*.o

ifeq ($(USE_NTL),1)
simple-offline.x: $(COMMON) simple-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

pairwise-offline.x: $(COMMON) pairwise-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

cnc-offline.x: $(COMMON) cnc-offline.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLIBS)

endif

clean:
	-rm */*.o *.o */*.d *.d *.x core.* *.a gmon.out */*/*.o
