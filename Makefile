
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)

ionizer: ionizer.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o ionizer  ionizer.cc $(ROOTLIBS)
	@echo 'done: ionizer (with ROOT)'

ionizerpp: ionizer.cpp Makefile
	g++ -O2 -Wall -Wextra -std=c++17 -m64 -o ionizerpp  ionizer.cpp
	@echo 'done: ionizerpp (without ROOT)'
