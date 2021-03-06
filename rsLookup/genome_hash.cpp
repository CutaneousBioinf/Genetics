#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <cmath>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include "diskhash/diskhash.hpp"
#include "CLI11.hpp"
using namespace std;
using namespace dht;

// Genome rsID hash table program

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;

const int MAX_DATA_LENGTH = 60;
const int MAX_BIG_DATA_LENGTH = 60;

// Defines the struct containing genome data struct
struct SNPData {
	char data[MAX_DATA_LENGTH];
};

// snp150Common column indices (zero-indexed)
const int chromosome_col = 1;
const int startpos_col = 2;
const int rsid_col = 4;
const int ref_allele_col = 8;
const int alleles_col = 9;

// rsID hash table maximum key length
const int key_maxlen = 15;
const int key_maxlen_big = 20;


// Check if string in list
bool in_list(vector<string> s, string key) {
	for (size_t i = 0; i < s.size(); ++i) {
		if (s[i] == key) {
			return true;
		}
	}
	return false;
}

// Check if vector<string> a is a subset of b
bool is_subset(vector<string> a, vector<string> b) {
	for (size_t i = 0; i < a.size(); ++i) {
		bool found = false;
		for (size_t j = 0; j < b.size(); ++j) {
			if (a[i] == b[j]) {
				found = true;
				break;
			}
		}
		if (!found) return false;
	}
	return true;
}

// Join string with separator
string join_str(const vector<string> &pieces, const string &separator) {
	string str = "";
	for (int i = 0; i < pieces.size(); ++i) {
		if (str == "") str = pieces[i];
		else str += separator + pieces[i];
	}
	return str;
}

// Splits string by delimiter into vector of strings
vector<string> str_split(const std::string& str, char delim = ' ')
{
	vector<string> cont;
	size_t current, previous = 0;
	current = str.find(delim);
	while (current != string::npos) {
		cont.push_back(str.substr(previous, current - previous));
		previous = current + 1;
		current = str.find(delim, previous);
	}
	cont.push_back(str.substr(previous, current - previous));
	return cont;
}

// Creates the hash table from the tab-separated file using RSID number string as key
int create_rsid_table(const char *source_name, const char* rsid_table_name, ostream *log_file) {
	cout << "Creating table from rsID to position...\n";
	DiskHash<SNPData> ht(rsid_table_name, key_maxlen, dht::DHOpenRW);
	string line;
	ifstream filebase;
	filebase.open(string(source_name));
	if (!filebase.is_open()) {
		cout << "Source file could not be opened" << endl;
		exit(1);
	}
	bool compressed = boost::filesystem::path(string(source_name)).extension() == ".gz";
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	if (compressed) {
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(filebase);
	}
	istream compressed_stream(&inbuf);
	istream &file = compressed ? compressed_stream : filebase;
	int max_rsid_length = 0;
	int max_content_length = 0;
	while (std::getline(file, line)) {
		SNPData item;
		vector<string> tsv = str_split(line, '\t');
		string chromosome = tsv[chromosome_col];
		string startpos = tsv[startpos_col];
		string rsid_str = tsv[rsid_col];
		if (rsid_str.length() - 2 > max_rsid_length) {
			max_rsid_length = rsid_str.length() - 2;
		}
		string alleles = tsv[alleles_col];
		string stored_data = chromosome + "\t" + startpos + "\t" + alleles;
		if (stored_data.length() > max_content_length) {
			max_content_length = stored_data.length();
		}
		char stored_arr[MAX_DATA_LENGTH] = {};
		strcpy(stored_arr, stored_data.c_str());
		strcpy(item.data, stored_arr);
		const char *rsid_num_part = rsid_str.substr(2, rsid_str.length()).c_str();
		const bool inserted = ht.insert(rsid_num_part, item);
	}
	filebase.close();
	return 0;
}

// Convert chromosome string "chr[1...22,x,y}" into shorter string for storage purposes
string chromosome_to_key(string chromosome) {
	if (chromosome.length() < 4 || chromosome.length() > 5) return "-1";
	if (chromosome.substr(0, 3) != "chr") return "-1";
	string chr_part = chromosome.substr(3, 2);
	int num = 0;
	if (tolower(chr_part.at(0)) == 'x') {
		num = 23;
	}
	else if (tolower(chr_part.at(0)) == 'y') {
		num = 24;
	}
	else {
		num = atoi(chr_part.c_str());
	}
	return to_string(num);
}

// Returns whether string s is nitrogenous base (A, C, T, G)
bool is_allele_letter(const string &s) {
	return s == "A" || s == "C" || s == "G" || s == "T";
}

// Check for valid allele sequence
bool is_allele(const string &s) {
	if (s.length() == 0) return false;
	for (size_t i = 0; i < s.length(); ++i) {
		string c(1, s[i]);
		if (!is_allele_letter(c)) return false;
	}
	return true;
}

// Get complement of allele (A->T, C->G)
string get_complement_letter(string allele) {
	if (allele == "A") return "T";
	else if (allele == "C") return "G";
	else if (allele == "G") return "C";
	else if (allele == "T") return "A";
	else return "";
}

// Get complement of allele
string get_complement(const string &allele) {
	string complement = "";
	for (size_t i = 0; i < allele.length(); ++i) {
		string c(1, allele[i]);
		complement += get_complement_letter(c);
	}
	return complement;
}

// Get complement of vector of allele sequences
vector<string> get_complement(const vector<string> &alleles) {
	vector<string> complement;
	for (size_t i = 0; i < alleles.size(); ++i) {
		complement.push_back(get_complement(alleles[i]));
	}
	return complement;
}

bool all_alleles(const vector<string> &alleles) {
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (!is_allele(alleles[i])) return false;
	}
	return alleles.size() > 0;
}

// Check if allele string is SNP on the condition that all slash-separated values are either A, C, G, or T
bool is_snp(const string &alleles_str) {
	if (alleles_str == "") return false;
	vector<string> alleles = str_split(alleles_str, '/');
	return all_alleles(alleles);
}

// Check if allele sequence is basic indel 
bool is_indel(const string &alleles_str) {
	if (alleles_str == "" || alleles_str.find("/") == string::npos) return false;
	vector<string> alleles = str_split(alleles_str, '/');
	return alleles.size() >= 2 && alleles[0] == "-" && all_alleles(vector<string>(alleles.begin() + 1, alleles.end()));
}

// Check if user inDEL input matches short form marker in source file
int is_indel_shortform(const string &user_input, const vector<string> &marker_split) {
	if (marker_split.size() != 2) return 0;
	if (marker_split[0].length() >= marker_split[1].length()) return 0;
	if (marker_split[0] == marker_split[1].substr(0, marker_split[0].length())) {
		if (user_input == "-/" + marker_split[1].substr(marker_split[0].length(), marker_split[1].length() - marker_split[0].length()))
			return marker_split[0].length();
	}
	return 0;
}

int is_indel_shortform_loop(const string &user_input, const vector<string> &marker_split) {
	vector<string> allele_split = str_split(user_input, '/');
	if (marker_split.size() < 2 || allele_split.size() < 2 || allele_split[0] != "-") return false;
	for (size_t i = 1; i < allele_split.size(); ++i) {
		int result = is_indel_shortform(allele_split[0] + "/" + allele_split[i], marker_split);
		if (result) return result;
	}
	return 0;
}

int is_indel_shortform_loop_backward(const string &user_input, vector<string> marker_split) {
	string back = marker_split[marker_split.size() - 1];
	marker_split.pop_back();
	marker_split.insert(marker_split.begin(), back);
	vector<string> allele_split = str_split(user_input, '/');
	if (marker_split.size() < 2 || allele_split.size() < 2 || allele_split[0] != "-") return false;
	for (size_t i = 1; i < allele_split.size(); ++i) {
		int result = is_indel_shortform(allele_split[0] + "/" + allele_split[i], marker_split);
		if (result) return result;
	}
	return 0;
}

int check_indel_shortform(const string &user_input, const vector<string> &marker_split) {
	return max(is_indel_shortform_loop(user_input, marker_split), is_indel_shortform_loop_backward(user_input, marker_split));
}

int indel_shortform_normalize(const vector<string> &marker_split) {
	if (marker_split.size() != 2) return 0;
	if (marker_split[0].length() >= marker_split[1].length()) return 0;
	if (marker_split[0] == marker_split[1].substr(0, marker_split[0].length())) {
		return marker_split[0].length();
	}
	return 0;
}

int indel_shortform_normalize_backward(const vector<string> &marker_split) {
	if (marker_split.size() != 2) return 0;
	if (marker_split[1].length() >= marker_split[0].length()) return 0;
	if (marker_split[1] == marker_split[0].substr(0, marker_split[1].length())) {
		return marker_split[1].length();
	}
	return 0;
}

// Creates the reverse map from chromosome/position/allele to rsID
int create_cpa_table_pointer(const char *source_name, const char* rsid_table_name, bool include_all, ostream *log_file) {
	cout << "Creating table from position to rsID...\n";
	// Create diskhash file plus a text file with ".data" appended to filename
	// The diskhash will point to locations in the .data file.
	DiskHash<size_t> ht(rsid_table_name, key_maxlen_big, dht::DHOpenRW);
	string line;
	ifstream filebase;
	filebase.open(string(source_name));
	if (!filebase.is_open()) {
		cout << "Source file could not be opened" << endl;
		exit(1);
	}
	bool compressed = boost::filesystem::path(string(source_name)).extension() == ".gz";
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	if (compressed) {
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(filebase);
	}
	istream compressed_stream(&inbuf);
	istream &file = compressed ? compressed_stream : filebase;
	ofstream s_file;
	s_file.open(string(rsid_table_name) + ".data");
	int max_rsid_length = 0;
	int max_content_length = 0;
	size_t prev_pos = s_file.tellp();
	string prev_chromosome = "";
	string prev_position = "";
	string line_break = "";
	// Loop through each line of snp150Common.txt.
	// Only recognized markers will be added to the table unless -A is supplied.
	while (std::getline(file, line)) {
		vector<string> tsv = str_split(line, '\t');
		string chromosome = tsv[chromosome_col];
		string chromosome_key = chromosome_to_key(chromosome);
		string startpos = tsv[startpos_col];
		string rsid_str = tsv[rsid_col];
		string ref_allele = tsv[ref_allele_col];
		string alleles = tsv[alleles_col];
		string rsid_num = rsid_str.substr(2, rsid_str.length());
		string b_key = "";
		// tellp returns the position in the data file that will be pointed to with the hash code
		if (include_all || 
		(is_snp(alleles) || is_indel(alleles))) {
			if (prev_chromosome != chromosome || prev_position != startpos) {
				s_file << line_break;
				b_key = "sc" + chromosome_key + "_" + startpos;
				prev_pos = s_file.tellp();
			prev_position = startpos;
			prev_chromosome = chromosome;
				const bool inserted = ht.insert(b_key.c_str(), prev_pos);
				line_break = "\n";
			}
			else s_file << "\t";
			s_file << rsid_str << "\t" << alleles;
		}
		else {
			if (log_file)
				(*log_file) << "Excluded " << rsid_str << " " << chromosome << " " << startpos << " " << alleles << endl;
		}
	}
	s_file << line_break;
	filebase.close();
	s_file.close();
	return 0;
}

// Gets the rsID from the hash table and prints to standard out the chromosome, starting position, and alleles
int get_rsid(const char *rsid_table_name, const char* rsid_param, ostream *log_file) {
	// Format check rsNNNNNNNNN
	if (!(strlen(rsid_param) > 2 && tolower(rsid_param[0]) == 'r' && tolower(rsid_param[1]) == 's')) {
		cout << "Invalid RSID requested" << endl;
		return 1;
	}
	DiskHash<SNPData> ht(rsid_table_name, key_maxlen, dht::DHOpenRO);
	string rsid(rsid_param);
	const char *rsid_num_part = rsid.substr(2, rsid.length()).c_str();
	SNPData *member;
	if (ht.is_member(rsid_num_part)) {
		member = ht.lookup(rsid_num_part);
		string result = member->data;
		vector<string> pieces = str_split(result, '\t');
		string chromosome = pieces[0];
		string startpos = to_string(atoi(pieces[1].c_str()) + 1);
		string alleles = pieces[2];
		cout << rsid << "\t" << chromosome << "\t" << startpos << "\t" << alleles << endl;
	}
	else {
		if (log_file)
			(*log_file) << rsid << " not found" << endl;
		cout << "ERROR: The given RSID was not found" << endl;
		return 1;
	}
	return 0;
}

// Look up chromosome/position/alleles and print RSID(s) [pointer version]
int get_cpa_pointers(const char *data_file_name, const char *rsid_table_name, const char * chromosome_c, const char *position_c, const char *allele_c, ostream *log_file) {
	DiskHash<size_t> ht(rsid_table_name, key_maxlen_big, dht::DHOpenRO);
	ifstream file;
	file.open(string(rsid_table_name) + ".data");
	string chromosome(chromosome_c);
	string allele(allele_c);
	vector<string> alleles_in = str_split(allele, '/');
	int position = atoi(position_c) - 1 + max(indel_shortform_normalize(alleles_in), indel_shortform_normalize_backward(alleles_in));
	string snp_b_key = "sc" + chromosome_to_key(chromosome) + "_" + to_string(position);
	int found = 0;
	if (!ht.is_member(snp_b_key.c_str())) {
		cout << "No matches found" << endl;;
		return 1;
	}
	// Use seekg to go to pointed-to position, which is the beginning of a line read with getline
	size_t *item = ht.lookup(snp_b_key.c_str());
	string line;
	file.seekg(*item);
	getline(file, line);
	vector<string> pieces = str_split(line, '\t');
	// Loop through tab-separated values
	for (int i = 0; i < pieces.size() / 2; ++i) {
		vector<string> alleles = str_split(pieces[2*i + 1], '/');
		vector<string> complement = get_complement(alleles);
		string rsid_num = pieces[2*i];
		if (allele == "" || 
		pieces[2*i + 1] == allele || 
		(!is_subset({"-"}, alleles) && (is_subset(alleles_in, alleles) || is_subset(alleles_in, complement))) ||
		(alleles[0] == "-" && alleles_in[0] == "-" && is_subset(vector<string>(alleles_in.begin() + 1, alleles_in.end()), vector<string>(alleles.begin() + 1, alleles.end()))) ||
		(alleles[0] == "-" && alleles_in[alleles_in.size() - 1] == "-" && is_subset(vector<string>(alleles_in.begin(), alleles_in.end() - 1), vector<string>(alleles.begin() + 1, alleles.end()))) || 
		check_indel_shortform(pieces[2*i + 1], alleles_in)) {
			string print_allele = allele;
			if (allele == "") print_allele = pieces[2*i + 1];
			cout << rsid_num << "\t" << chromosome << "\t" << (position + 1 - check_indel_shortform(pieces[2*i + 1], alleles_in)) << "\t" << print_allele << endl;
			++found;
		}
		else {
			if (log_file) *(log_file) << chromosome << " " << (position + 1) << " " << pieces[2*i + 1] << " did not match with input " << allele << endl;
		}
	}
	file.close();
	if (!found) cout << "No matches found" << endl;
	return 0;
}

// Takes rsIDs from an input file, line-separated, and queries the hash table to output chromosome/position/allele data
void get_rsid_file(const string &filename, const string &source, ostream *log_file) {
	ifstream file;
	file.open(filename);
	string rsid;
	while (file >> rsid) {
		get_rsid(source.c_str(), rsid.c_str(), log_file);
	}
	file.close();
}

// Takes chromosome/position/allele data from an input file, tab-separated, one per line, and queries the hash table to output rsID data
void get_cpa_pointers_file(const string &filename, const string &source, const string &table, ostream *log_file) {
	ifstream file;
	file.open(filename);
	string chr, pos, allele_seq;
	while (file >> chr >> pos >> allele_seq) {
		get_cpa_pointers(source.c_str(), table.c_str(), chr.c_str(), pos.c_str(), allele_seq.c_str(), log_file);
	}
	file.close();
}

void load_positional_arguments(const vector<string> &args, const vector<string*> &target) {
	size_t nextpos = 0;
	for (size_t i = 0; i < target.size(); ++i) {
		string &deref = *(target[i]);
		if (deref.length() == 0) {
			if (nextpos == args.size()) break;
			deref = args[nextpos++];
		}
	}
}

// Main function
// Uses CLI11 for the command-line interface, which autogenerates a help screen
int main(int argc, char** argv) {
	CLI::App app{"rsLookup genomic hash table lookup program"};
	string type = "";
	string pos1 = "", pos2 = "", pos3 = "", pos4 = "";
	app.add_option("pos1", pos1, "pos1")->group("");
	app.add_option("pos2", pos2, "pos2")->group("");
	app.add_option("pos3", pos3, "pos3")->group("");
	app.add_option("pos4", pos4, "pos4")->group("");
	CLI::Option *set_option = app.add_set("-k,--kind", type, {"rsid", "cpa"}, "Kind of creation or lookup (rsID or chromosome/position/allele)");
	string table_name = "", source_name = "";
	CLI::Option *source_option = app.add_option("-s,--source", source_name, "Path to source data file used for hash table creation")->check(CLI::ExistingFile);
	CLI::Option *table_option = app.add_option("-t,--table", table_name, "Path to hash table to be created or examined");
	string folder_name = "";
	CLI::Option *folder_option = app.add_option("-d,--directory", folder_name, "Path to directory of hash table data to be created or examined, storing both types")->excludes(table_option);
	bool create = false, retrieve = false;
	CLI::Option *create_flag = app.add_flag("-C,--create", create, "Create hash table\nNOTE: In addition to command-line options, positional arguments can also be used for the lookups as follows:\nrsLookup -C [path_to_snp150Common.txt[.gz]] [target_directory_path]");
	CLI::Option *retrieve_flag = app.add_flag("-R,--retrieve", retrieve, "Retrieves a given rsID (starting with rs) or a chromosome/position/allele triplet separated by spaces\nNOTE: In addition to command-line options, positional arguments can also be used for the lookups as follows:\nrsLookup [data_directory_path_or_appropriate_table] rsNNNNNNNNN\nrsLookup [data_directory_path_or_appropriate_table] chrXX [position] [allele_sequence]\nrsLookup [data_directory_path_or_appropriate_table] [chr]:[pos] [a1/a2]");
	string file_path = "";
	CLI::Option *file_option = app.add_option("-f,--file", file_path, "Path to input file")->check(CLI::ExistingFile);
	string rsid = "", chromosome = "", position = "", alleles = "";
	CLI::Option *rsid_option = app.add_option("-r,--rsid", rsid, "rsID marker (format rsNNNNNNNNN)");
	CLI::Option *chromosome_option = app.add_option("-c,--chromosome", chromosome, "Chromosome (format chr{1,...,26,X,Y})");
	CLI::Option *allele_option = app.add_option("-a,--allele", alleles, "Allele to be looked up");
	CLI::Option *position_option = app.add_option("-p,--position", position, "Position to be looked up, indexed from 1");
	string log_name = "";
	bool quiet = false, verbose = false;
	CLI::Option *log_file_option = app.add_option("-l,--log-file", log_name, "Specifies the location of a log file to contain entries of a given chromosome and position that do not match the allele input");
	CLI::Option *quiet_option = app.add_flag("-q,--quiet", quiet, "Quiet mode");
	CLI::Option *verbose_option = app.add_flag("-v,--verbose", verbose, "Verbose mode")->excludes(quiet_option);
	bool include_all = false;
	CLI::Option *all_option = app.add_flag("-A,--all", include_all, "Include all markers in reverse map, not just recognized ones");
	CLI11_PARSE(app, argc, argv);
	ostream *log_file = nullptr;
	ofstream log_stream;
	if (log_name.length() > 0) {
		log_stream.open(log_name);
		log_file = &log_stream;
	}
	else if ((verbose && create) || (!create && !quiet)) log_file = &cerr;
	
	
	// Parse positional arguments
	vector<string> allArgs = {pos1, pos2, pos3, pos4};
	vector<string> args;
	for (size_t i = 0; i < allArgs.size(); ++i) {
		if (allArgs[i].length() > 0) args.push_back(allArgs[i]);
	}
	if (!create) retrieve = true;
	
	if (create) {
		bool targetFromPos = table_name.length() == 0 && folder_name.length() == 0;
		load_positional_arguments(args, {&source_name, &table_name});
		if (table_name.length() > 0 && folder_name.length() > 0) {
			cerr << "Error: Only one destination can be specified\n";
			return 1;
		}
		if (targetFromPos && boost::filesystem::path(table_name).extension() == "")
			folder_name = table_name;
		if (type.length() == 0 && !boost::filesystem::is_directory(folder_name))
			boost::filesystem::create_directories(folder_name);
		
		if (type == "rsid")
			create_rsid_table(source_name.c_str(), table_name.c_str(), log_file);
		else if (type == "cpa")
			create_cpa_table_pointer(source_name.c_str(), table_name.c_str(), include_all, log_file);
		else if (type.length() == 0 && folder_name.length() > 0) {
			create_rsid_table(source_name.c_str(), (folder_name + "/rsid.dht").c_str(), log_file);
			create_cpa_table_pointer(source_name.c_str(), (folder_name + "/cpa.dht").c_str(), include_all, log_file);
		}
		else {
			cerr << "The request was unrecognized. No action was taken.\n";
			return 1;
		}
	}
	else if (retrieve) {
		if (table_name.length() == 0 && folder_name.length() == 0) {
			if (args.size() > 0) {
				string location = args[0];
				args.erase(args.begin());
				if (boost::filesystem::path(location).extension() == "") folder_name = location;
				else table_name = location;
			}
			else {
				int help_argc = 2;
		const char *help_argv[2] = {"rsLookup", "-h"};
		CLI11_PARSE(app, help_argc, help_argv);
				return 1;
			}
		}
		if (type == "") {
			if (rsid.length() > 0) type = "rsid";
			else if (chromosome.length() > 0) type = "cpa";
			else {
				if (args.size() > 0) {
					if (args[0].substr(0, 2) == "rs") type = "rsid";
					else type = "cpa";
				}
				else {
					cerr << "Error: Lookup not specified\n";
					return 1;
				}
			}
		}
		if (type == "rsid") {
			load_positional_arguments(args, {&rsid});
			if (table_name.size() == 0) table_name = folder_name + "/rsid.dht";
			if (file_path == "") get_rsid(table_name.c_str(), rsid.c_str(), log_file);
			else get_rsid_file(file_path, table_name, log_file);
		}
		else if (type == "cpa") {
			load_positional_arguments(args, {&chromosome, &position, &alleles});
			if (table_name.size() == 0) table_name = folder_name + "/cpa.dht";
			if (chromosome.substr(0, 3) != "chr") chromosome = "chr" + chromosome;
			vector<string> chr_position = str_split(chromosome, ':');
			if (chr_position.size() > 1) {
				chromosome = chr_position[0];
				if (position.length() > 0 && alleles.length() == 0) alleles = position;
				position = chr_position[1];
				vector<string> position_alleles = str_split(position, '_');
				if (position_alleles.size() > 1) {
					position = position_alleles[0];
					alleles = position_alleles[1];
				}
			}
			if (file_path == "") get_cpa_pointers(source_name.c_str(), table_name.c_str(), chromosome.c_str(), position.c_str(), alleles.c_str(), log_file);
			else get_cpa_pointers_file(file_path, source_name, table_name, log_file);
		}
	}
	else {
		int help_argc = 2;
		const char *help_argv[2] = {"rsLookup", "-h"};
		CLI11_PARSE(app, help_argc, help_argv);
	}
}
