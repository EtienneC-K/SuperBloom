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
    for minimizer_size in [7, 9, 11, 13, 15]:
        for nb_blocks in range(2*minimizer_size-2, 2*minimizer_size+3):
            block_size = 33-nb_blocks
            command = f"./target/release/bloomybloom -t 16 --input-type 1 -m {minimizer_size} --block-size {block_size} target/release/SRR_first_trad.fasta"
            start = perf_counter()
            os.system(command)
            execution_time = perf_counter()-start
            all_results.apppend(f"m : {minimizer_size}, nb_blocks : {nb_blocks}; executed in {execution_time}")
            print(all_results[-1])

    print()
    print()
    for result in all_results:
        print(result)

    

if __name__ == "__main__":
    main()
