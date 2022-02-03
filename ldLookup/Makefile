DEFAULT_GOAL = release
CXXFLAGS_ALL = -std=c++17 -Wextra -Wall -Wpedantic -Isrc
CXXFLAGS_DEBUG = -g3 -fsanitize=address,leak,undefined
CXXFLAGS_PROFILE = -g3 -pg
CXXFLAGS_RELEASE = -O3

BINARY_NAME = ldLookup
BINARY_DEPS = src/diskhash/diskhash.o main.cpp $(wildcard src/*.cpp)
TESTS_NAME = tests
TESTS_DEPS = src/diskhash/diskhash.o $(wildcard tst/*.cpp)

release:
	@$(MAKE) -C src/diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_RELEASE) $(BINARY_DEPS) -o $(BINARY_NAME)
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_RELEASE) $(TESTS_DEPS) -o $(TESTS_NAME)

debug:
	@$(MAKE) -C src/diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_DEBUG) $(BINARY_DEPS) -o $(BINARY_NAME)
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_DEBUG) $(TESTS_DEPS) -o $(TESTS_NAME)

profile:
	@$(MAKE) -C src/diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_PROFILE) $(BINARY_DEPS) -o $(BINARY_NAME)
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_PROFILE) $(TESTS_DEPS) -o $(TESTS_NAME)

clean:
	find . -name '*.o *.out *.so' -delete
	rm -f src/diskhash/disktest $(BINARY_NAME) $(TESTS_NAME)

.PHONY: release debug profile clean