DEFAULT_GOAL = release
CXXFLAGS_ALL = -std=c++17 -Wextra -Wall -Wpedantic
CXXFLAGS_DEBUG = -g3
CXXFLAGS_PROFILE = $(CXXFLAGS_DEBUG) -pg
CXXFLAGS_RELEASE = -o3

BINARY_NAME = ldLookupx
BINARY_DEPS = diskhash/diskhash.o geneticData.hpp main.cpp tables.cpp vdh.cpp
TESTS_NAME = testsx
TESTS_DEPS = diskhash/diskhash.o geneticData.hpp tables.cpp tests.cpp vdh.cpp

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
	find . -name '*.o' -delete
	find . -name '*.so' -delete
	rm -f diskhash/disktest $(BINARY_NAME) $(TESTS_NAME)
