DEFAULT_GOAL = release
CXXFLAGS_ALL = -std=c++17 -Wextra -Wall -Wpedantic
CXXFLAGS_DEBUG = -g3
CXXFLAGS_PROFILE = $(CXXFLAGS_DEBUG) -pg
CXXFLAGS_RELEASE = -o3
BINARY_NAME = ldLookupx

release:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_RELEASE) diskhash/diskhash.o ldLookup.cpp -o $(BINARY_NAME)

debug:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_DEBUG) diskhash/diskhash.o ldLookup.cpp -o $(BINARY_NAME)

profile:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_PROFILE) diskhash/diskhash.o ldLookup.cpp -o $(BINARY_NAME)


clean :
	find . -name '*.o' -delete
	find . -name '*.so' -delete
	rm -f diskhash/disktest $(BINARY_NAME)
