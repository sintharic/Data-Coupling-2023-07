COMPILER := /usr/local/Cellar/open-mpi/4.1.5/bin/mpic++
INCLUDE  := /usr/local/include
CPPFLAGS := -O2 -lHYPRE -ferror-limit=32 -std=c++11
HYPRELIB := /usr/local/Cellar/hypre/2.28.0/lib
LDFLAGS  := "-L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++"
LLVMLIB  := /usr/local/opt/llvm/lib/c++

default_target: parallel

parallel:
	$(COMPILER) $(LDFLAGS) -I $(INCLUDE) -L $(HYPRELIB) $(CPPFLAGS) demo_mpi.cpp -o currentCalc.exe