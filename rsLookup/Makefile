DEFAULT_GOAL = release
CXXFLAGS_DEBUG = -ldht -g
CXXFLAGS_RELEASE = -ldht -o3

# Compile the main executable
release: genome_hash.cpp
	g++ $(CXXFLAGS_RELEASE) genome_hash.cpp -o rsLookup

debug:
	g++ $(CXXFLAGS_DEBUG) genome_hash.cpp -o rsLookup


# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump

