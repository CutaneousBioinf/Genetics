DEFAULT_GOAL = release
CXXFLAGS_DEBUG = -g3
CXXFLAGS_RELEASE = -o3
CXXFLAGS_LIBS = -lboost_system -lboost_iostreams -lboost_filesystem

# Compile the main executable
release: genome_hash.cpp
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_RELEASE) diskhash/diskhash.o genome_hash.cpp -o rsLookup $(CXXFLAGS_LIBS)

debug:
	@$(MAKE) -C diskhash
	g++ $(CXXFLAGS_DEBUG) diskhash/diskhash.o genome_hash.cpp -o rsLookup $(CXXFLAGS_LIBS)


# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump

