import collections
import os.path as path
import shlex
import subprocess
import time
import typing


####################################
#####     GLOBAL CONSTANTS     #####
####################################

# ldLookup Directory
PARENT_DIR = path.dirname(path.dirname(path.abspath(__file__)))
# Fully-qualified path to the ldLookup executable
EXE = path.join(PARENT_DIR, 'ldLookup')
# Fully-qualified path to the test LD data.
TEST_DATA = path.join(PARENT_DIR, 'tst', 'test.ld')
# Used to measure execution time
GLOBAL_START_TIME = time.time()

#########################
#####     TYPES     #####
#########################

MARKER_ID = str
MARKER_LIST = typing.List[MARKER_ID]
MARKER_ASSOCIATIONS = typing.Mapping[MARKER_ID, MARKER_LIST]
DIFF = typing.Tuple[
    MARKER_LIST,
    MARKER_LIST,
    MARKER_ASSOCIATIONS,
    MARKER_ASSOCIATIONS
]

class LDPair:
    def __init__(self,
                 index_marker: MARKER_ID,
                 ld_marker: MARKER_ID,
                 index_maf: float,
                 r_squared: float):
        self.index_marker = index_marker
        self.ld_marker = ld_marker
        self.index_maf = index_maf
        self.r_squared = r_squared

    def __hash__(self):
        return hash((self.index_marker, self.ld_marker, self.index_maf, self.r_squared))

    @staticmethod
    def from_string(s: str,
                    index_snp_column:int=2,
                    index_snp_maf_column:int=3,
                    ld_snp_column:int=6,
                    r_squared_column:int=7,
                    delimiter:str=' ') -> "LDPair":
        fields = s.split(delimiter)
        try:
            return LDPair(fields[index_snp_column],
                          fields[ld_snp_column],
                          float(fields[index_snp_maf_column]),
                          float(fields[r_squared_column]))
        except (IndexError, ValueError):
            return None

    @staticmethod
    def validate(ld_pair: "LDPair", min_r_squared: float = 0):
        return min_r_squared <= ld_pair.r_squared < 1


#####################################
#####     HELPER FUNCTIONS     ######
#####################################


'''Prints `s` with timing information.'''
def print_info(s: str):
    print(f"[{round(time.time() - GLOBAL_START_TIME, 5)}] {s}")


'''Prints `d` with one key-value pair per line.'''
def pretty_print_dict(d: typing.Dict):
    for k, v in d.items():
        for v_i in v:
            print("\t", k, " : ", v_i)


'''Prints list elements with formatting.'''
def pretty_print_list(lst: typing.List):
    for i in lst:
        print("\t", i)


'''
EFFECTS: Executes `command.` If `silent` is False, prints timing information.
`command` can be a subprocess-executable list or string.
RETURNS: A tuple containing the command's standard output (str), standard
error (str), and return code (int).
'''
def call(command, silent: bool=False) -> typing.Tuple[str, str, int]:
    t0 = time.time()
    if not silent:
        print_info(f"Running '{command}'...")

    # Convert strings to lists so subprocess doesn't choke
    if type(command) == str:
        command = shlex.split(command)

    # Execute the command, capturing stderr and stdout as strings
    proc = subprocess.Popen(command,
                            encoding='utf-8',
                            stderr=subprocess.PIPE, 
                            stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    if not silent:
        delta = round(time.time() - t0, 5)
        print_info(f"Return Code: {proc.returncode} After {delta}s")

    return stdout, stderr, proc.returncode


'''Prints its inputs and exits the test suite.'''
def fail(cmd, out: str, err: str):
    print_info("TEST FAILED")
    print_info("COMMAND:")
    print(cmd)
    print_info("STDOUT:")
    print(out)
    print_info("STDERR:")
    print(err)
    exit(1)


'''
EFFECTS: Attempts to parse and validate each element of `string_iterable` as an
LDPair.
RETURNS: Tuple of well-formed valid LDPairs, well-formed invalid (i.e. excluded)
LDPairs, and ill-formed strings from `string_iterable`.
'''
def parse_ld_pairs(string_iterable: typing.Iterable[str],
                   pair_builder=LDPair.from_string,
                   pair_validator=LDPair.validate) -> typing.List[LDPair]:
    valid_pairs = []
    excluded_pairs = []
    malformed_lines = []
    for line in string_iterable:
        pair = pair_builder(line)
        if pair is None:
            malformed_lines.append(line)
        elif not pair_validator(pair):
            excluded_pairs.append(pair)
        else:
            valid_pairs.append(pair)
    return valid_pairs, excluded_pairs, malformed_lines


'''Parses strings output by ldLookup to associations genetic markers.'''
def parse_marker_assocs(ldLookup_output: str) -> MARKER_ASSOCIATIONS:
    # ldLookup output looks like
    # "key_marker1\tval_1_1\n\tkey_marker1\n\tvalue_1_2\n\tkey_marker2\tvalue_marker_2_1"...
    assocs = collections.defaultdict(list)
    for line in ldLookup_output.split('\n')[:-1]:
        try:
            key_marker, value_snp = line.split('\t')
            assocs[key_marker].append(value_snp)
        except ValueError:
            print_info(f"WARNING: Failed to Parse Marker Association '{line}'")
    
    # Convert assocs from defaultdict
    return dict(assocs)


'''
EFFECTS: Compares two marker associations. Ignores order of marker assocs.
RETURNS: A tuple containing:
    1) a list of keys appearing in the expected assoc, but not the actual
       assoc (i.e. a list of missing keys)
    2) a list of keys appearing in the actual assoc, but not the expected
       assoc (i.e. a list of extra keys)
    3) a mapping of keys in the expected assoc to markers that appear
       associated with those keys in the expected assoc, but not the actual
       assoc (i.e. missing associations)
    4) a mapping of keys in the actual assoc to markers that appear
       associated with those keys in the actual assoc, but not the expected
       assoc (i.e. extra associations)
'''
def compare_associations(expected: MARKER_ASSOCIATIONS,
                         actual: MARKER_ASSOCIATIONS) -> DIFF:
    # Copy actual so we don't modify it.
    actual = dict(actual)
    missing_keys = []
    extra_assoc = {}
    missing_assoc = collections.defaultdict(list)

    for marker, associations in expected.items():
        try:
            actual_assocs = actual[marker]
        except KeyError:
            missing_keys.append(marker)
            continue

        for assoc in associations:
            try:
                actual_assocs.remove(assoc)
            except ValueError:
                missing_assoc[marker].append(assoc)

        if actual_assocs:
            extra_assoc[marker] = actual_assocs
        del actual[marker]

    extra_keys = list(actual.keys())
    # Convert missing_assoc to a dict rather than a defaultdict.
    return missing_keys, extra_keys, dict(missing_assoc), extra_assoc


'''Prints a diff of two sets of marker associations.'''
def print_diff(diff: DIFF):
    missing_keys, extra_keys, missing_values, extra_values = diff
    print_info("Printing Diff...")
    print_info("MISSING KEY MARKERS")
    if missing_keys:
        pretty_print_list(missing_keys)
    else:
        print("No Missing Key Markers")

    print_info("EXTRA KEY MARKERS")
    if extra_keys:
        pretty_print_list(extra_keys)
    else:
        print("No Extra Key markers")

    print_info("MISSING ASSOCIATIONS")
    if missing_values:
        pretty_print_dict(missing_values)
    else:
        print("No Missing Associations")

    print_info("EXTRA ASSOCIATIONS")
    if extra_values:
        pretty_print_dict(extra_values)
    else:
        print("No Extra Associations")


########################################################
##########          TEST DRIVER CODE          ##########
########################################################


def load_test_data():
    try:
        f = open(TEST_DATA, 'r', encoding='utf-8')
        lines = f.readlines()
        f.close()
    except OSError:
        print_info(f"Failed to Open {TEST_DATA}")
        exit(1)
    return parse_ld_pairs(lines)


def test_get_ld(ld_surrogates: typing.Mapping[MARKER_ID, typing.Set[MARKER_ID]]):
    print_info("Testing get_ld Subcommand: test_get_ld")
    for index_snp, assoc_snps in ld_surrogates.items():
        cmd = [EXE, "get_ld", "test_dataset", "-m", index_snp]
        out, err, code = call(cmd, silent=True)
        if code:
            fail(cmd, out, err)
        
        parsed_out = parse_marker_assocs(out)
        expected = {index_snp : list(assoc_snps)}
        diff = compare_associations(expected, parsed_out)
        if any(diff):
            print_diff(diff)
            fail(cmd, out, err)
    
    print_info("Completed Test")


def test_get_ld_file(ld_surrogates: typing.Mapping[MARKER_ID, typing.Set[MARKER_ID]]):
    print_info("Testing get_ld Subcommand: test_get_ld")

    test_path = path.join(PARENT_DIR, 'tst', 'test.in')
    with open(test_path, 'w', encoding='utf-8') as of:
        for k in ld_surrogates.keys():
            of.write(k)
            of.write('\n')

    expected = {}
    for k, v in ld_surrogates.items():
        expected[k] = list(v)    

    cmd = [EXE, "get_ld", "test_dataset", test_path]
    out, err, code = call(cmd, silent=True)
    if code:
        fail(cmd, out, err)
        
    parsed_out = parse_marker_assocs(out)
    diff = compare_associations(expected, parsed_out)
    if any(diff):
        print_diff(diff)
        fail(cmd, out, err)
    
    print_info("Completed Test")


if __name__ == '__main__':
    valid, excluded, malformed = load_test_data()
    ld_surrogates = collections.defaultdict(set)
    for pair in valid:
        ld_surrogates[pair.index_marker].add(pair.ld_marker)
    
    print_info("Loaded Test Data")

    print_info("Testing create Subcommand")
    cmd = f"{EXE} create test_dataset {TEST_DATA}"
    out, err, code = call(cmd)
    if code:
        fail(cmd, out, err)

    test_get_ld(ld_surrogates)
    test_get_ld_file(ld_surrogates)
