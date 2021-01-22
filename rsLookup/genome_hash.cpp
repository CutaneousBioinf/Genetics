#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <cmath>

#include <diskhash.hpp>
#include "CLI11.hpp"
using namespace std;
using namespace dht;

// Genome RSID hash table program

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

struct SNPDataBig {
	char data[MAX_BIG_DATA_LENGTH];
};

// Constants
const char *source_name = "../mini.txt";
const char *rsid_table_name = "rsid.dht";
const int chromosome_col = 1;
const int startpos_col = 2;
const int rsid_col = 4;
const int ref_allele_col = 8;
const int alleles_col = 9;
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

// Converts int to char array containing binary digits that make up the int
void convert_int_to_char_array(unsigned int num, char* arr) {
	arr[0] = (num >> 24);
	arr[1] = (num >> 16);
	arr[2] = (num >> 8);
	arr[3] = (num >> 0);
	int params = 1;
	for (int i = 0; i <= 3; ++i) {
		arr[i] = arr[i] / 2;
		if (arr[i] == 0) {
			arr[i] = 1;
			params = params | ((unsigned int) pow(2.0, i + 1));
		}
	}
	arr[4] = (char) params;
}

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
int create_rsid_table(const char *source_name, const char* rsid_table_name) {
	DiskHash<SNPData> ht(rsid_table_name, key_maxlen, dht::DHOpenRW);
	string line;
	ifstream file;
	file.open(string(source_name));
	if (!file.is_open()) {
		cout << "Source file could not be opened" << endl;
		return 1;
	}
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
		if (stored_data.length() > MAX_DATA_LENGTH) {
			cout << "ERROR: Maximum length exceeded, encountered string of length " << stored_data.length() << endl;
		}
		char stored_arr[MAX_DATA_LENGTH] = {};
		strcpy(stored_arr, stored_data.c_str());
		strcpy(item.data, stored_arr);
		const char *rsid_num_part = rsid_str.substr(2, rsid_str.length()).c_str();
		cout << "adding " << rsid_str << "..." << endl;
		const bool inserted = ht.insert(rsid_num_part, item);
		if (!inserted) {
			cout << "[already in table]" << endl;
		}
	}
	file.close();
	cout << "Max RSID length (excluding rs): " << max_rsid_length << endl;
	cout << "Max content length: " << max_content_length << endl;
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

// Returns whether string s is allele (A, C, T, G)
bool is_allele(string s) {
	return s == "A" || s == "C" || s == "G" || s == "T";
}

// Get complement of allele (A->T, C->G)
string get_complement(string allele) {
	if (allele == "A") return "T";
	else if (allele == "C") return "G";
	else if (allele == "G") return "C";
	else if (allele == "T") return "A";
	else return "";
}

vector<string> get_complement(const vector<string> &alleles) {
	vector<string> complement;
	for (size_t i = 0; i < alleles.size(); ++i) {
		complement.push_back(get_complement(alleles[i]));
	}
	return complement;
}

// Check if allele string is SNP on the condition that all slash-separated values are either A, C, G, or T
bool is_snp(const string &alleles_str) {
	if (alleles_str == "" || alleles_str.find("/") == string::npos) return false;
	vector<string> alleles = str_split(alleles_str, '/');
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (!is_allele(alleles[i])) return false;
	}
	// if (alleles.size() > 0) cout << "true " << join_str(alleles, "/");
	return alleles.size() > 0;
}

string canonical_alleles(const string &alleles_str) {
	if (alleles_str == "" || alleles_str.find("/") == string::npos) return "";
	vector<string> alleles = str_split(alleles_str, '/');
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (!is_allele(alleles[i])) return "";
	}
	if (alleles.size() == 0) return "";
	vector<string> c_alleles;
	if (in_list(alleles, "A") || in_list(alleles, "T")) c_alleles.push_back("A");
	if (in_list(alleles, "C") || in_list(alleles, "G")) c_alleles.push_back("C");
	return join_str(c_alleles, "/");
}

// Counts the number of members of the hash table of key c_x, where c is constant and x is a number incremented by 1
int count_other_members(DiskHash<SNPData> *ht, string b_key) {
	int i = -1;
	while (ht->is_member((b_key + "_" + to_string(i + 1)).c_str())) {
		++i;
	}
	return i + 1;
}

int count_other_members_big(DiskHash<SNPDataBig> *ht, string b_key) {
	int i = -1;
	while (ht->is_member((b_key + "_" + to_string(i + 1)).c_str())) {
		++i;
	}
	return i + 1;
}

// Count other members pointers
int count_other_members_pointers(DiskHash<size_t> *ht, string b_key) {
	int i = -1;
	while (ht->is_member((b_key + "_" + to_string(i + 1)).c_str())) {
		++i;
	}
	return i + 1;
}

// Create chromosome/position/alleles table
int create_cpa_table(const char *source_name, const char* rsid_table_name) {
	DiskHash<SNPDataBig> ht(rsid_table_name, key_maxlen_big, dht::DHOpenRW);
	string line;
	ifstream file;
	file.open(string(source_name));
	if (!file.is_open()) {
		cout << "Source file could not be opened" << endl;
		return 1;
	}
	int max_rsid_length = 0;
	int max_content_length = 0;
	while (std::getline(file, line)) {
		SNPDataBig item;
		vector<string> tsv = str_split(line, '\t');
		string chromosome = tsv[chromosome_col];
		string chromosome_key = chromosome_to_key(chromosome);
		string startpos = tsv[startpos_col];
		string rsid_str = tsv[rsid_col];
		if (rsid_str.length() - 2 > max_rsid_length) {
			max_rsid_length = rsid_str.length() - 2;
		}
		string ref_allele = tsv[ref_allele_col];
		string alleles = tsv[alleles_col];
		string rsid_num = rsid_str.substr(2, rsid_str.length());
		string stored_data = rsid_num + "\t" + ref_allele + "\t" + alleles;
		if (stored_data.length() > max_content_length) {
			max_content_length = stored_data.length();
		}
		if (stored_data.length() > MAX_BIG_DATA_LENGTH) {
			cout << "ERROR: Maximum length exceeded, encountered string of length " << stored_data.length() << endl;
			continue;
		}
		string b_key = "";
		if (is_snp(alleles)) {
			char stored_arr[MAX_DATA_LENGTH] = {};
			strcpy(stored_arr, stored_data.c_str());
			strcpy(item.data, stored_arr);
			b_key = "sc" + chromosome_key + "_" + startpos;
		}
		if (b_key != "") {
			cout << "adding " << b_key << "_" << count_other_members_big(&ht, b_key) << "..." << endl;
			cout << file.tellg() << endl;
			const bool inserted = ht.insert((b_key + "_" + to_string(count_other_members_big(&ht, b_key))).c_str(), item);
			if (!inserted) {
				cout << "[already in table]" << endl;
			}
		}
	}
	file.close();
	cout << "Max RSID length (excluding rs): " << max_rsid_length << endl;
	cout << "Max content length: " << max_content_length << endl;
	return 0;
}


int create_cpa_table_pointer(const char *source_name, const char* rsid_table_name) {
	DiskHash<size_t> ht(rsid_table_name, key_maxlen_big, dht::DHOpenRW);
	string line;
	ifstream file;
	file.open(string(source_name));
	ofstream s_file;
	s_file.open(string(rsid_table_name) + ".data");
	if (!file.is_open()) {
		cout << "Source file could not be opened" << endl;
		return 1;
	}
	int max_rsid_length = 0;
	int max_content_length = 0;
	size_t prev_pos = s_file.tellp();
	string prev_chromosome = "";
	string prev_position = "";
	string line_break = "";
	while (std::getline(file, line)) {
		vector<string> tsv = str_split(line, '\t');
		string chromosome = tsv[chromosome_col];
		string chromosome_key = chromosome_to_key(chromosome);
		string startpos = tsv[startpos_col];
		string rsid_str = tsv[rsid_col];
		if (rsid_str.length() - 2 > max_rsid_length) {
			max_rsid_length = rsid_str.length() - 2;
		}
		string ref_allele = tsv[ref_allele_col];
		string alleles = tsv[alleles_col];
		string c_alleles = canonical_alleles(alleles);
		string rsid_num = rsid_str.substr(2, rsid_str.length());
		string b_key = "";
		if (is_snp(alleles)) {
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
	}
	s_file << line_break;
	file.close();
	s_file.close();
	cout << "Max RSID length (excluding rs): " << max_rsid_length << endl;
	cout << "Max content length: " << max_content_length << endl;
	return 0;
}


// Gets the RSID from the hash table and prints to standard out the chromosome, starting position, and alleles
int get_rsid(const char *rsid_table_name, const char* rsid_param) {
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
		cout << member->data << endl;
	}
	else {
		cout << "ERROR: The given RSID was not found" << endl;
		return 1;
	}
	return 0;
}

// Lookup chromosome/position/alleles and print RSID(s)
int get_cpa(const char *rsid_table_name, const char * chromosome_c, const char *position_c, const char *allele_c) {
	DiskHash<SNPDataBig> ht(rsid_table_name, key_maxlen_big, dht::DHOpenRO);
	string chromosome(chromosome_c);
	int position = atoi(position_c) - 1;
	string allele(allele_c);
	string snp_b_key = "sc" + chromosome_to_key(chromosome) + "_" + to_string(position);
	SNPDataBig *member;
	int found = 0;
	for (int i = 0; i < count_other_members_big(&ht, snp_b_key); ++i) {
		SNPDataBig *item = ht.lookup((snp_b_key + "_" + to_string(i)).c_str());
		vector<string> pieces = str_split(item->data, '\t');
		string rsid_num = pieces[0];
		string ref_allele = pieces[1];
		vector<string> alleles = str_split(pieces[2], '/');
		if (in_list(alleles, ref_allele) && (in_list(alleles, allele) || in_list(alleles, get_complement(allele)))) {
			cout << "rs" << rsid_num << endl;
			++found;
		}
	}
	if (!found) cout << "No matches found" << endl;
	return 0;
}

// Look up chromosome/position/alleles and print RSID(s) [pointer version]
int get_cpa_pointers(const char *data_file_name, const char *rsid_table_name, const char * chromosome_c, const char *position_c, const char *allele_c) {
	DiskHash<size_t> ht(rsid_table_name, key_maxlen_big, dht::DHOpenRO);
	ifstream file;
	file.open(string(rsid_table_name) + ".data");
	string chromosome(chromosome_c);
	int position = atoi(position_c) - 1;
	string allele(allele_c);
	vector<string> alleles_in = str_split(allele, '/');
	string snp_b_key = "sc" + chromosome_to_key(chromosome) + "_" + to_string(position);
	int found = 0;
	if (!ht.is_member(snp_b_key.c_str())) {
		cout << "No matches found" << endl;;
		return 1;
	}
	size_t *item = ht.lookup(snp_b_key.c_str());
	string line;
	file.seekg(*item);
	getline(file, line);
	vector<string> pieces = str_split(line, '\t');
	for (int i = 0; i < pieces.size() / 2; ++i) {
		vector<string> alleles = str_split(pieces[2*i + 1], '/');
		vector<string> complement = get_complement(alleles);
		string rsid_num = pieces[2*i];
		if (is_subset(alleles_in, alleles) || is_subset(alleles_in, complement)) {
			cout << rsid_num << endl;
			++found;
		}
	}
	file.close();
	if (!found) cout << "No matches found" << endl;
	return 0;
}

// Command line interface
/*
int main(int argc, char** argv) {
	string main_command = "";
	if (argc >= 2) {
		main_command = string(argv[1]);
		if (main_command == "create_rsid_table") {
			if (argc >= 4) {
				return create_rsid_table(argv[2], argv[3]);
			}
			else {
				cout << "Please follow create_rsid_table by the file name of the source and destination" << endl;
				return 1;
			}
		}
		else if (main_command == "create_cpa_table") {
			if (argc >= 4) {
				return create_cpa_table(argv[2], argv[3]);
			}
			else {
				cout << "Please follow create_cpa_table by the file name of the source and destination" << endl;
				return 1;
			}
		}
		else if (main_command == "create_cpa_table_pointer") {
			if (argc >= 4) {
				return create_cpa_table_pointer(argv[2], argv[3]);
			}
			else {
				cout << "Please follow create_cpa_table_pointers by the file name of the source and destination" << endl;
				return 1;
			}
		}
		else if (main_command == "get_rsid") {
			if (argc >= 4) {
				return get_rsid(argv[2], argv[3]);
			}
			else {
				cout << "Please follow get_rsid by the file name of the table and the RSID in form RSxxxxxxxxxxxx" << endl;
				return 1;
			}
		}
		else if (main_command == "get_cpa") {
			if (argc >= 6) {
				return get_cpa(argv[2], argv[3], argv[4], argv[5]);
			}
			else {
				cout << "Please follow get_cpa by the file name of the table, chromosome (chr[1...22,x,y]), position number, and allele" << endl;
				return 1;
			}
		}
		else if (main_command == "get_cpa_pointer") {
			if (argc >= 7) {
				return get_cpa_pointers(argv[2], argv[3], argv[4], argv[5], argv[6]);
			}
			else {
				cout << "Please follow get_cpa_pointer by the file name of the dataset, file name of the table, chromosome (chr[1...22,x,y]), position number, and allele" << endl;
				return 1;
			}
		}
	}
	else {
		cout << "USAGE: [create_rsid_table|create_cpa_table] [source_tsv] [target_dht]" << endl;
		cout << "get_rsid [source_dht] [rsid]" << endl;
	}
	return 0;
}
*/


int main(int argc, char** argv) {
	CLI::App app{"Genomic hash table lookup program"};
	string type = "";
	CLI::Option *set_option = app.add_set("-k,--kind", type, {"rsid", "cpa"}, "Kind of creation or lookup (rsID or chromosome/position/allele)");
	string table_name = "", source_name = "";
	CLI::Option *source_option = app.add_option("-s,--source", source_name, "Path to source data file used for hash table creation")->check(CLI::ExistingFile);
	CLI::Option *table_option = app.add_option("-t,--table", table_name, "Path to hash table to be created or examined");
	bool create = false, retrieve = false;
	CLI::Option *create_flag = app.add_flag("-C,--create", create, "Create hash table")->needs(set_option)->needs(source_option)->needs(table_option);
	CLI::Option *retrieve_flag = app.add_flag("-R,--retrieve", retrieve, "Retrieves a given rsID (starting with rs) or a chromosome/position/allele triplet separated by spaces")->needs(set_option)->needs(source_option)->needs(table_option);
	string rsid = "", chromosome = "", position = "", alleles = "";
	CLI::Option *rsid_option = app.add_option("-r,--rsid", rsid, "rsID marker (format rsNNNNNNNNN)");
	CLI::Option *chromosome_option = app.add_option("-c,--chromosome", rsid, "Chromosome (format chr{1,...,26,X,Y})");
	CLI::Option *allele_option = app.add_option("-a,--allele", alleles, "Allele to be looked up");
	CLI::Option *position_option = app.add_option("-p,--position", position, "Position to be looked up, indexed from 1");
	CLI11_PARSE(app, argc, argv);
	if (create)
		if (type == "rsid")
			create_rsid_table(source_name.c_str(), table_name.c_str());
		else if (type == "cpa")
			create_cpa_table_pointer(source_name.c_str(), table_name.c_str());
	else if (retrieve)
		if (type == "rsid")
			get_rsid(table_name.c_str(), rsid.c_str());
		else if (type == "cpa")
			get_cpa_pointers(source_name.c_str(), table_name.c_str(), chromosome.c_str(), position.c_str(), alleles.c_str());
}
