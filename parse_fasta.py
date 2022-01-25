#!/usr/bin/env python3

"""
Kate Gerhardt (kgerhar4)
parse_fasta.py

Parses a fasta file with the overall focus of obtaining information relevant to the
creation of a fragment recruitment plot. Collects headers, lengths of individual
entries, and total length. Prints output to a file in CSV format.

Accepts FASTA formatted files, can be gz compressed.

Assumptions:
    -Input files are in FASTA
    -Files have headers beginning with >
    -File has extension .fa , .fna ,  or .fasta (with .gz at end if compressed)


Run script like this:
python3 parse_fasta.py -i genome.fna -o **** output file name

"""

import re
import os
import gzip
import argparse
import fragment_recruit
import csv

def main():
    parser = argparse.ArgumentParser(description='xxx')
    parser.add_argument('-r', '--ref_genome', type=str, required=True, help='Path to file representing an input reference genome')
    parser.add_argument('-a', '--align_data', type=str, required=True, help='Path to SAM formatted file containing alignment data')
    args = parser.parse_args()

    filename = args.ref_genome

    # Checks if file is gzip compressed
    if args.ref_genome.endswith('.gz'):
        # If gzip compressed, use gzip.open to read without uncompressing
        inp = gzip.open(args.ref_genome, 'rt')

        # Store name of input file as string with .gz removed
        filename = args.ref_genome[0:-3]

    else:
        # If not compressed, use open
        inp = open(args.ref_genome)

    # Parses FASTA files
    number_molecules = 0
    seq_length = 0
    reference_details = []

    if filename.endswith('.fa') or filename.endswith('.fasta') or filename.endswith('.fna'):
        for line in inp:
            line = line.rstrip()

            if line.startswith('>'):
                if (seq_length > 0):
                    reference_details.append(seq_length)
                    seq_length = 0
                # Increments total sequence counter if fasta header character detected
                number_molecules += 1
                reference_details.append(line)
            else:
                seq_length += len(line)

    reference_details.append(seq_length)

    # Calculate total reference sequence length
    total_length = 0
    for x in reference_details:
        if (isinstance(x, int)):
            total_length += x

    reference_details.append(total_length)
    reference_details.append(number_molecules)


    with open("reference_details.csv", "w", newline="\n") as file:
        writer = csv.writer(file)
        writer.writerow(reference_details)

    # Closes files
    inp.close()


if __name__ == '__main__':
    main()
