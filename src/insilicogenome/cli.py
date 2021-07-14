import getopt
import sys
import argparse
from insilicogenome import insilicogenome

#def parse_args():
#    """
#    This is a help message made to make user understand how to use the script
#    """

help_message = """
usage: insilicogenome -o|--output <string> [-s|--size] <integer>

A Python utility to create artificial genomes, this tool is designed to be used in bioinformatics benchmarking programs.

Arguments :
    -h, --help    show this help message and exit
    -o, --output  Name of the genome and header [string]
    -s, --size    Size of the genome [base paire]

 Examples
    --------
    insilicogenome -o genome.fasta -s 100
"""
try:
    options, args = getopt.getopt(sys.argv[1:], "ho:s:", ["help", "output=", "size="])
except getopt.GetoptError as err:
    print(str(err))
    print(help_message)
    sys.exit(2)

for opt, arg in options:
    if opt in ("-h", "--help"):
        print(help_message)
        sys.exit()
    elif opt in ("-o", "--output"):
        output = arg
    elif opt in ("-s", "--size"):
        size = int(arg)
    else:
        print(help_message)
        sys.exit(2)

def main():
    sequence = insilicogenome.random_dnasequence(size)
    insilicogenome.write_fasta_genome(output, sequence)
    print(f"A genome of {size}bp have been write in {output}")

if __name__ == "__main__":
    main()