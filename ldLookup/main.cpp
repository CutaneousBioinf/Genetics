#include "ldLookup.hpp"

int main (int argc, char** argv) {
    const std::string usage = "USAGE: ldlookupx create [table_name] [source_file]\n       ldlookupx get [table_name] [key]\n";
    if (argc != 4) {
        std::cout << usage;
        return 1;
    }

    try {
        const std::string command = std::string(argv[1]);
        const std::string table = std::string(argv[2]);
        const std::string arg = std::string(argv[3]);
        if (!command.compare("create")) {
            RecordParser p;
            LDTable::create_table(table, arg, p);
        } else if (!command.compare("get")) {
            LDTable ldt(table);
            for (std::string s : ldt.get(arg)) {
                std::cout << "Value: " << s << "\n";
            }
        } else {
            std::cout << usage;
            return 1;
        }
    } catch (std::exception &e) {
        std::cout << "ERROR: " << e.what() << "\n";
        return 1;
    }
}