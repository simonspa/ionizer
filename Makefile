
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)

ionizer: ionizer.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o ionizer  ionizer.cc $(ROOTLIBS)
	@echo 'done: ionizer'

mfp: mfp.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o mfp  mfp.cc $(ROOTLIBS)
	@echo 'done: mfp'
