#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <cmath>

#include <diskhash.hpp>
using namespace std;
using namespace dht;

// Genome RSID hash table program

const int MAX_DATA_LENGTH = 60;

// Defines the struct containing genome data struct
struct SNPData {
	char data[MAX_DATA_LENGTH];
};

// Constants
const char *source_name = "../mini.txt";
const char *rsid_table_name = "rsid.dht";
const int chromosome_col = 1;
const int startpos_col = 2;
const int rsid_col = 4;
const int alleles_col = 9;
const int key_maxlen = 15;

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

// Command line interface
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
		else if (main_command == "get_rsid") {
			if (argc >= 4) {
				return get_rsid(argv[2], argv[3]);
			}
			else {
				cout << "Please follow get_rsid by the file name of the table and the RSID in form RSxxxxxxxxxxxx" << endl;
				return 1;
			}
		}
	}
	else {
		cout << "USAGE: create_rsid_table [source_tsv] [target_dht]" << endl;
		cout << "get_rsid [source_dht] [rsid]" << endl;
	}
	return 0;
}
