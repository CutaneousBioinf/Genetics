CXX := g++
CXXFLAGS := -Isrc -std=c++17 -Wall -Werror -Wextra -Wpedantic
PRODFLAGS := -O3
DEBUGFLAGS := -g -fsanitize=address,leak,undefined
PROFFLAGS := -g3 -pg

SRC := $(wildcard src/*.cpp) src/diskhash/src/diskhash.o
TST := $(wildcard tst/*.cpp) $(SRC)

build: diskhash build_ldLookup tests

debug: diskhash debug_ldLookup tests

profile: diskhash profile_ldLookup tests

diskhash:
	cd src/diskhash/src && $(MAKE)

build_ldLookup: main.cpp $(SRC)
	$(CXX) $(CXXFLAGS) $(PRODFLAGS) main.cpp $(SRC) -o ldLookup

debug_ldLookup: main.cpp $(SRC)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) main.cpp $(SRC) -o ldLookup

profile_ldLookup: main.cpp $(SRC)
	$(CXX) $(CXXFLAGS) $(PROFFLAGS) main.cpp $(SRC) -o ldLookup

tests: $(TST)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(TST) -o tests

pull_diskhash:
	git clone https://github.com/luispedro/diskhash.git
	rm -rf diskhash/.travis diskhash/haskell diskhash/python
	rm diskhash/.travis.yml diskhash/diskhash.cabal diskhash/Makefile diskhash/MANIFEST.in diskhash/README.md diskhash/Setup.hs diskhash/setup.py diskhash/stack.yaml

clean :
	rm -rvf *.exe *~ *.so *.o *.out *.dSYM *.stackdump
	rm -f ldLookup tests

.PHONY: build debug profile build_ldLookup debug_ldLookup profile_ldLookup tests pull_diskhash clean