DEFAULT_GOAL = release
CXXFLAGS_ALL = -std=c++17 -Wextra -Wall -Wpedantic
CXXFLAGS_DEBUG = -g3 -fsanitize=address,leak,undefined
CXXFLAGS_PROFILE = -g3 -pg
CXXFLAGS_RELEASE = -o3

BINARY_NAME = ldLookup
BINARY_DEPS = diskhash/diskhash.o geneticData.hpp main.cpp tables.cpp utils.cpp vdh.cpp
TESTS_NAME = tests
TESTS_DEPS = diskhash/diskhash.o geneticData.hpp tables.cpp tests.cpp utils.cpp vdh.cpp

release:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_RELEASE) $(BINARY_DEPS) -o $(BINARY_NAME)
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_RELEASE) $(TESTS_DEPS) -o $(TESTS_NAME)

debug:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_DEBUG) $(BINARY_DEPS) -o $(BINARY_NAME)
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_DEBUG) $(TESTS_DEPS) -o $(TESTS_NAME)

profile:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_PROFILE) $(BINARY_DEPS) -o $(BINARY_NAME)
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_PROFILE) $(TESTS_DEPS) -o $(TESTS_NAME)

clean:
	find . -name '*.o *.out *.so' -delete
	rm -f diskhash/disktest $(BINARY_NAME) $(TESTS_NAME)
