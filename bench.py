#!/usr/bin/env python3

from sys import argv
import csv
import subprocess

#function that launches all the tests and writes in the output file
def main():
    #start by reading all arguments
    input_file, threads, max_ram = parse_arguments()

    #"header" of the csv
    data = write_csv_header(input_file, threads, max_ram)

    #launching kmc with time
    data.append(run_kmc(input_file, threads))
    print("KMC done")

    #launching gerbil with time
    data.append(run_gerbil(input_file, threads, max_ram))
    print("Gerbil done")

    #writing a "header" transition before running bloomybloom
    data = write_transitional_header(data)

    #launching all different bloomybloom with many different options
    data = launch_bloomys(data, input_file, threads, max_ram)

    #writing to the csv
    write_data_to_csv(data)


def parse_arguments():
    """parses the 3 arguments we are looking for"""
    if len(argv) != 4 :
        #incorrect use of the programm, print help message
        print_help_message()

    #if we do have the correct number of arguments, try to convert each one to the desired type
    input_file = argv[1]
    if not (input_file.endswith("fasta") or input_file.endswith("fastq")):
        print("WARNING : the file pointed to doesn't end with fasta or fastq, are you sure it's the right one ?")

    threads = 0
    try :
        threads = int(argv[2])
    except (TypeError, ValueError) :
        print("ERROR : number of threads isn't an integer")
        print()
        print_help_message()

    max_ram = 0
    try :
        max_ram = int(argv[3])
        assert (max_ram >= 12)
    except (TypeError, ValueError, AssertionError) :
        print("ERROR : amount of available RAM (in GB) isn't an integer greater or equal than 12")
        print()
        print_help_message()

    return (input_file, threads, max_ram)


def print_help_message():
    """simple function to print what the programm excpects"""
    print("To use this benchmark programm use the following 3 arguments :")
    print("    path to fasta file to read")
    print("    number of threads to use")
    print("    number of available GB of RAM, must be an integer >= 12")
    exit(1)


def write_csv_header(input_file, threads, max_ram):
    """returns the "data" list for the first couple of cells of the csv"""
    #as a reminder, each item of data represents a line of the final file
    file_name = input_file.split("/")[-1]
    data = [[f"Source file name : {file_name}"]]
    data.append([f"threads : {threads}, available RAM : {max_ram}GB"])
    data.append([])
    data.append([])
    data.append(["", "wall clock", "user time", "system time", "%CPU", "", "major page faults", "minor page faults", "swaps"])

    return data


def run_kmc(input_file, threads):
    """rusn kmc on the specified file, and returns all metrics and text wanted
    in the correct order to be appended to data
    Has a failsafe in case kmc (like always) fails to parse the file"""
    #TODO : i'll also have to run the ram only kmc sometime
    executable_path = "./KMC3.2.4.linux.x64/bin/kmc"
    completed_run = None
    command = f"\\time {executable_path} {input_file} results/31mers temps -t {threads}"
    try :
        completed_run = subprocess.run(command, shell = True, capture_output = True)
    except Exception as e :
        print("KMC problem")
        print(e)
        return ["KMC", "died"]

    #apparently i have to catch kmc exceptions this way
    if completed_run.stderr.decode("utf-8").startswith("Error") or \
        completed_run.stderr.decode("utf-8").startswith("\nError"):
        print("KMC problem")
        print(completed_run.stderr.decode("utf-8"))
        return ["KMC", "died"]

    print("Run data de KMC :")
    print(completed_run.stderr)
    print("#######################")

    run_data = parse_backslash_time(completed_run.stderr)

    #now adding "" to align the colums properly
    run_data.insert(0, "KMC")
    run_data.insert(5, "")

    return run_data


def run_gerbil(input_file, threads, max_ram):
    #runs gerbil, with 31mers, restricts max ram in GB
    executable_path = "./gerbil/build/gerbil"
    completed_run = None
    command = f"\\time {executable_path} {input_file} -k 31 -t {threads} -e {max_ram}GB temps/ results/gerbil_res"

    try :
        completed_run = subprocess.run(command, shell = True, capture_output = True)
    except Exception as e :
        print("Gerbil problem")
        print(e)
        return ["Gerbil", "died"]

    #Gerbil Errors might be like kmc's ones, caught like this
    if completed_run.stderr.decode("utf-8").startswith("ERROR") or \
        completed_run.stderr.decode("utf-8").startswith("\nError"):
        print("Gerbil problem")
        print(completed_run.stderr.decode("utf-8"))
        return ["Gerbil", "died"]

    run_data = parse_backslash_time(completed_run.stderr)

    #now adding "" to align the colums properly
    run_data.insert(0, "Gerbil")
    run_data.insert(5, "")

    return run_data


def write_transitional_header(data):
    """writes a couple of lines to make the transition between kmc/gerbil and all the bloomybloom runs"""
    data.append(["#"*100])
    data.append(["Bloomybloom"])
    header_columns = [
        "", 
        "wall_clock", 
        "user time", 
        "system time", 
        "%CPU",
        "",
        "major page faults",
        "minor page faults",
        "swaps",
        "",
        "non zero block rate",
        "max block fill",
        "average fill",
        "median fill",
        "",
        "hash table non zeros",
        "max block fill",
        "average fill"
        "median fill"
    ]
    data.append(header_columns)

    return data


def parse_backslash_time(output: str):
    """parses the output of \time, returning a list with :
        wall clock time,
        user time,
        system time,
        %CPU,
        major page faults,
        minor page faults,
        swaps,
    """
    #first split everything
    splited = output.split()

    #then extract string containing the time stamps and such
    wall_clock = splited[2].decode("utf-8").replace("elapsed", "")
    wall_clock = wall_clock.replace(":", "m") #because can auto format to a date otherwise
    #y'aura ptet un probleme la dedans si j'ai des tests qui depasse les 1h
    user = splited[0].decode("utf-8").replace("user", "")
    system = splited[1].decode("utf-8").replace("system", "")
    cpu = splited[3].decode("utf-8").replace("%CPU", "")
    swaps = splited[8].decode("utf-8").replace("swaps", "")

    #pagefualts need a bit more work to parse
    pagefaults=splited[7].decode("utf-8").replace("(", "").split("+")
    major = pagefaults[0].replace("major", "")
    minor = pagefaults[1].replace("minor)pagefaults", "")

    return [wall_clock, user, system, cpu, major, minor, swaps]


def write_data_to_csv(data):
    """writes the list of list, data, to bench_results.csv"""
    output_path = "bench_results.csv"

    #clear out any None that might be left by failing functions
    for (i, sub_list) in enumerate(data):
        if sub_list is None:
            data[i] = [""]

    #actual writing part
    with open(output_path, mode='w', newline='') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(data)
    print(f"wrote all results in {output_path}")


def launch_bloomys(data, input_file, threads, max_ram):
    """
    launches all different options of bloomybloom we want to test
    always launches them twice, one without counting of fill rates to have the \time results
    and one with, to have the fillings

    no_options, default values
    """
    #function structure follows : write the options down in data and in a variable
    #run and add the results to "data"

    #base values of the options first
    size = 33
    block_size = 14
    ht_size = 28
    ht_block_size = 14
    no_ht = False
    no_bloom = False
    minimizer_size = 11
    n_hashes = 7

    #no_options default values
    data, options = update_options(data, threads, size, block_size, ht_size, ht_block_size, no_ht, no_bloom, minimizer_size, n_hashes)
    data.append(launch_and_collect(input_file, options))
    print("Finished a bloomybloom option set")

    #dotn forget the return
    return (data)


def update_options(data, threads, size, block_size, ht_size, ht_block_size, no_ht, no_bloom, minimizer_size, n_hashes):
    """updates the options and data variable with all the specified options values"""
    options = f"-t {threads} --input-type 1 --size {size} --block-size {block_size} --table-size {ht_size} --table-block-size {ht_block_size}"
    if no_bloom:
        options += " --no-bloom"
    elif no_ht:
        options += " --no-hashtable"

    data.append(write_options(size, block_size, ht_size, ht_block_size, no_ht, no_bloom, minimizer_size))

    return (data, options)


def launch_and_collect(input_file, options):
    """launches bloomybloom twice, first time collects the times and such, second time collets fillings"""
    executable_path = "./target/release/bloomybloom"
    timed_results = None
    counted_results = None

    command = f"\\time {executable_path} {options} {input_file}"
    try :
        completed_timed_run = subprocess.run(command, shell = True, capture_output = True)
        timed_results = parse_backslash_time(completed_timed_run.stderr)
        timed_results.insert(0, "")
        timed_results.insert(5, "")
    except Exception as e :
        print("Bloomybloom problem")
        print(f"Option set : {options}")
        print(f"full command : {command}")
        print(e)
        timed_results = ["Failed"]*8

    #now for the counting version
    options += " --counting --auto-bench"
    command = f"{executable_path} {options} {input_file}"
    try :
        compelted_counted_run = subprocess.run(command, shell = True, capture_output = True)
        counted_results = parse_counted(compelted_counted_run.stdout)
    except Exception as e :
        print("Bloomybloom problem :")
        print(e)
        counted_results = ["Failed"]*8

    return (timed_results+counted_results)


def write_options(size, block_size, ht_size, ht_block_size, no_ht, no_bloom, minimizer_size):
    """
    function to write a human readable string with all options, this string will be appended to data
    important note : this writes all sizes fully, but the given variables are powers of 2
    """
    
    return_text = f"minimizer_size : {minimizer_size}"
    return_text += f", size {2**size}"
    return_text += f", block-size : {2**block_size}"
    return_text += f", nb-blocks{2**(size-block_size)}"
    return_text += f", ht-size : {2**ht_size}"
    return_text += f", ht-block-size : {2**ht_block_size}"
    return_text += f", nb-ht_blocks : {2**(ht_size-ht_block_size)}"
    return_text += f", no_hashtable : {no_ht}"
    return_text += f", no_bloom : {no_bloom}"

    return ([return_text])


def parse_counted(count_output):
    """reads the first line of the counted output, the one with all the slashes
    also writes down the data truncated to not overwelm output with too many decimal places"""
    data_line = count_output.decode("utf-8").splitlines()[0]
    data_line = data_line.split("|")
    return_list = []
    for i in range(8):
        #main use of this loop is to check the length of data_line, ie if no errors in bloomybloom
        data_bit = data_line[i]
        data_bit = data_bit[0:6]
        return_list.append(data_bit)

    #no to insert spaces in the list to align the columns
    return_list.insert(0, "")
    return_list.insert(5, "")

    return return_list




if __name__ == "__main__":
    main()
