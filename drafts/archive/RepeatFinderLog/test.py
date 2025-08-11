import timeit
import sys
from findRepeatsV4 import RepeatFinder, NoRepeats

rf = RepeatFinder()

str1 = "AGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCT"
str2 = "AGCTTAGCTAAGCTTAGCTTAGCTTAGCTAAGCTTAGCTGAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAGCTTAGCTAAGCTTAGCTAAGCTTAGCTAAGCTTAGCT"
sequence = "GGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGCCAAAGACCCGTGGTGTTGCCATGAATCCAGTTGATCACACTCACGGTGGTGGTAACCATCAACATATTGTTAAGGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGC"
quality = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,FF,:,F,,F:F,:F:FF,F,FFFFFFFFFF,,FF:F,,FFFF,:FFFFF,F,,F,FFFFF:FFFFFF,,FFFF::FFFFF:::F:FF:,FFFF:FFFFF:F,F:,FFF,,FF,FFF,:F:F::FFFFFFF,F"
try:
    repeat_distance = rf.find_repeat_dist(sequence, k=20, max_errors=2)
except NoRepeats:
    print("No Repeats Found")
    
errors = 2

def test_compare_string():
    rf.compare_string(str1, str2, errors)

def test_compare_string_v2():
    rf.compare_string_v2(str1, str2, errors)
    
def test_compare_string_hd():
    rf.compare_string_hd(str1, str2, errors)

def test_consensus_v0():
    rf.consensus_from_long_read(sequence, quality, repeat_distance)

def test_consensus_v1():
    rf.consensus_from_long_read_fast(sequence, quality, repeat_distance)

def test_consensus_np():
    rf.consensus_from_long_read_numpy(sequence, quality, repeat_distance)

# Number of runs from command line
repeat_count = int(sys.argv[1]) if len(sys.argv) > 1 else 1
number = 10000

for run in range(repeat_count):
    # time1 = timeit.timeit(test_compare_string, number=number)
    # time2 = timeit.timeit(test_compare_string_v2, number=number)
    # time3 = timeit.timeit(test_compare_string_hd, number=number)
    # print(f"Run {run+1}:")
    # print(f"  compare_string:    {time1:.5f} seconds for {number} calls")
    # print(f"  compare_string_v2: {time2:.5f} seconds for {number} calls")
    # print(f"  compare_string_hd: {time3:.5f} seconds for {number} calls")
    time1 = timeit.timeit(test_consensus_v0, number=number)
    time2 = timeit.timeit(test_consensus_v1, number=number)
    time3 = timeit.timeit(test_consensus_np, number=number)
    print(f"Run {run+1}:")
    print(f"  consensus_v0: {time1:.5f} seconds for {number} calls")
    print(f"  consensus_v1: {time2:.5f} seconds for {number} calls")
    print(f"  consensus_np: {time3:.5f} seconds for {number} calls")
