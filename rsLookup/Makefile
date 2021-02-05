DEFAULT_GOAL = release
CXXFLAGS_DEBUG = -g
CXXFLAGS_RELEASE = -o3

# Compile the main executable
release: genome_hash.cpp
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_RELEASE) diskhash/diskhash.o genome_hash.cpp -o rsLookup

debug:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_DEBUG) diskhash/diskhash.o genome_hash.cpp -o rsLookup


# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump

