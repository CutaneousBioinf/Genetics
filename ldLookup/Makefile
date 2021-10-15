DEFAULT_GOAL = release
CXXFLAGS_ALL = -std=c++17 -Wextra -Wall -Wpedantic
CXXFLAGS_DEBUG = -g3
CXXFLAGS_PROFILE = $(CXXFLAGS_DEBUG) -pg
CXXFLAGS_RELEASE = -o3

release:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_RELEASE) diskhash/diskhash.o ldlookup.cpp -o ldlookupx

debug:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_DEBUG) diskhash/diskhash.o ldlookup.cpp -o ldlookupx

profile:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_ALL) $(CXXFLAGS_PROFILE) diskhash/diskhash.o ldlookup.cpp -o ldlookupx


clean :
	rm -vf ldlookupx
