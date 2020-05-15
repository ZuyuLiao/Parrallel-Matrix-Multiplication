CLASS_DIR=/project/linuxlab/class/cse539
CC=$(CLASS_DIR)/tapir/build/bin/clang
CXX=$(CLASS_DIR)/tapir/build/bin/clang++

# vanilla cilk plus runtime location
CILK_LIBS=$(CLASS_DIR)/cilkplus-rts/libs

CFLAGS = -ggdb -O3 -fcilkplus
CXXFLAGS = -ggdb -O3 -fcilkplus
LIBS = -L$(CILK_LIBS) -Wl,-rpath -Wl,$(CILK_LIBS) -lcilkrts -lpthread -lrt -lm
PROGS = mm_dac mm_dac_inst

INST_RTS_LIBS=$(CLASS_DIR)/cilkplus-rts-instrumented/lib
INST_RTS_INCLUDE=$(CLASS_DIR)/cilkplus-rts-instrumented/include
INST_LIBS=-L$(INST_RTS_LIBS) -Wl,-rpath -Wl,$(INST_RTS_LIBS) -lcilkrts -lpthread -lrt -lm -ldl
INST_FLAGS=-DINSTRUMENT_RTS=1 -I$(INST_RTS_INCLUDE)

all:: $(PROGS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

%.inst.o: %.c
	$(CC) $(CFLAGS) $(INST_FLAGS) -o $@ -c $<

%.inst.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INST_FLAGS) -o $@ -c $<

mm_dac: ktiming.o getoptions.o mm_dac.o
	$(CXX) -o $@ $^ $(LIBS)

mm_dac_inst: ktiming.inst.o getoptions.inst.o mm_dac.inst.o
	$(CXX) -o $@ $^ $(INST_LIBS)


clean::
	-rm -f $(PROGS) *.o *.inst.o
