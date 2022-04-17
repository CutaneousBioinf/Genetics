// #include <cassert>
// #include <string>

// #include "vdh.hpp"

// const std::string TEST_FILE = "exclude/test_file.vdhdat";
// const std::string TEST_TABLE = "exclude/test_table.vdhdht";
// const size_t TEST_KEY_SIZE = 16;

// const Options TEST_TABLE_OPTS {
//     TEST_FILE, TEST_TABLE, TEST_KEY_SIZE, false
// };

// const Options TEST_TABLE_CREATE_OPTS {
//     TEST_FILE, TEST_TABLE, TEST_KEY_SIZE, true
// };


// void check_append_empty_key(const VectorDiskHash& vdh) {

// }

// void check_append_max_key(const VectorDiskHash& vdh) {

// }

// void check_append_empty_value(const VectorDiskHash& vdh) {

// }

// void check_append_max_value_to_reserved_non_tail_key(const VectorDiskHash& vdh) {

// }

// void check_append_max_value_to_reserved_tail_key(const VectorDiskHash& vdh) {

// }

// void check_append_overlong_key(const VectorDiskHash& vdh) {
    
// }

// void check_append_overlong_value_to_reserved_non_tail_key(const VectorDiskHash& vdh) {

// }

// void check_append_overlong_value_to_reserved_tail_key(const VectorDiskHash& vdh) {

// }

// void check_non_tail_key(const VectorDiskHash& vdh) {

// }

// void check_reserve_empty_key(const VectorDiskHash& vdh) {
    
// }

// void check_reserve_max_key(const VectorDiskHash& vdh) {
    
// }

// void check_reserve_overlong_key(const VectorDiskHash& vdh) {
    
// }

// void check_reserved_tail_key(const VectorDiskHash& vdh) {

// }

// void check_reserved_non_tail_key(const VectorDiskHash& vdh) {

// }

// void check_tail_key(const VectorDiskHash& vdh) {

// }

// void check_create_duplicate_table() {
//     Options opts = TEST_TABLE_CREATE_OPTS;
//     try {
//         VectorDiskHash vdh(opts);
//         assert(false);
//     } catch (std::exception& ignore) {
//         ;
//     }

//     Options opts_2 = opts;
//     opts_2.file_path = "exclude/path_that_doesnt_exist_ae0a0c98";
//     try {
//         VectorDiskHash vdh(opts_2);
//         assert(false);
//     } catch (std::exception& ignore) {
//         ;
//     }

//     Options opts_3 = opts;
//     opts_3.table_path = "exclude/path_that_doesnt_exist_ae0a0c98";
//     try {
//         VectorDiskHash vdh(opts_3);
//         assert(false);
//     } catch (std::exception& ignore) {
//         ;
//     }
// }

// void check_open_existing_table_with_incorrect_key_size() {
//     Options opts = TEST_TABLE_OPTS;
//     opts.max_key_size = 32;

//     try {
//         VectorDiskHash vdh(opts);
//         assert(false);
//     } catch (std::exception& ignore) {
//         ;
//     }
// }

// void check_can_open_existing_table() {
//     Options opts = TEST_TABLE_OPTS;
//     VectorDiskHash test_vdh(TEST_TABLE_OPTS);

//     opts.max_key_size = 0;
//     VectorDiskHash test_vdh2(opts);
// }

// int main() {
//     // Destroy test_vdh before running table opening tests.
//     {
//         VectorDiskHash test_vdh(TEST_TABLE_CREATE_OPTS);
//     }
    
//     check_create_duplicate_table();
//     check_open_existing_table_with_incorrect_key_size();
//     check_can_open_existing_table();
// }

int main() {
    return 0;
}