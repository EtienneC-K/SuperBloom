#!/usr/bin/env python3

######################
# This module is here to make plenty of benchmarks with many parameters, and is meant to be modified 
# often
#####################

from time import perf_counter
import os

def main():
    all_results = []
    #testing min_size vs nb blocks
    for minimizer_size in [7, 9, 11, 13]:
    #for minimizer_size in [7, 9, 11, 13, 15]:
        #for nb_blocks in range(2*minimizer_size-2, 2*minimizer_size+3):
        for nb_blocks in range(2*minimizer_size, 2*minimizer_size+1):
            times = []
            for _ in range (5): #to get better averages
                block_size = 37-nb_blocks
                command = f"\\time ./target/release/bloomybloom -t 16 --input-type 1 -m {minimizer_size} --block-size {block_size} --size 37 target/release/SRR_first_trad.fasta"
                start = perf_counter()
                os.system(command)
                execution_time = perf_counter()-start
                times.append(execution_time)
            average_time = sum(times)/len(times)
            all_results.append(f"m : {minimizer_size}, nb_blocks : {nb_blocks}; executed in {average_time}")
            print(all_results[-1])

    print()
    print()
    for (i, result) in enumerate(all_results):
        if i%5 == 0:
            print()
        print(result)

    

if __name__ == "__main__":
    main()
