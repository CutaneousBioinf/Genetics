#ifndef LDLOOKUP_GLOBAL_HPP
#define LDLOOKUP_GLOBAL_HPP
#define CHECK_FAIL(file, msg) if ((file).fail()) { throw std::runtime_error((msg)); }
#endif