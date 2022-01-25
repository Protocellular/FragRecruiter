#!/usr/bin/env python3

"""
Kate Gerhardt (kgerhar4)
parse_sam

Parses a SAM file with the overall focus of obtaining information relevant to the
creation of a fragment recruitment plot. Prints output to a file in CSV format.


Accepts SAM formatted files

Assumptions:
    -Input files are in SAM format
    -Files have header (please see ********* samtools if file does not contain a header)


Run script like this:
python3 parse_sam.py -i alignment_data.sam -o **** output file name

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

    filename = args.align_data

    inp = open(args.align_data)

    accessions = []          # [acc_1, acc_1_len, ... acc_n, acc_n_len]
    all_reads = []           # each read [flag, rname, pos, mapq, cigar, segment_seq, align_len, percent_identity]
    tool_info = []           # [tool name, tool version]

    for line in inp:
        line = line.rstrip()

        # process headers
        if line.startswith('@'):

            # get accessions / program info
            if re.search(r'SQ', line):
                seq_info = []

                # Get accessions
                if re.search(r'SN:', line):
                    accession_obj = re.search(r'(?<=SN:).+(?=\t)', line)
                    accession = accession_obj.group()
                    seq_info.append(accession)

                # Get Length
                if re.search(r'LN:', line):
                    len_obj = re.search(r'(?<=LN:).+(?=\W)', line)
                    seq_info.append(len_obj)
                    accessions.append(seq_info)

            # get program info
            if re.search(r'PG', line):
                # get name of alignment program
                tool_name_obj = re.search(r'(?<=ID:).+(?=PN)', line)
                tool_name = tool_name_obj.group()
                tool_info.append(str(tool_name))

                # get version of alignment program
                tool_ver_obj = re.search(r'(?<=VN:).+(?=\t)', line)
                tool_ver = tool_ver_obj.group()
                tool_info.append(tool_ver)


        else:
            # process aligned reads

            # split line by tab
            data_line = re.split('\t', line)
            # print(data_line)

            # parse relevant info
            read_name = data_line[0]
            flag = data_line[1]
            ref_name = data_line[2]
            pos = data_line[3]
            mapq = data_line[4]
            cigar = data_line[5]
            segment_seq = data_line[9]


            align_len = 0
            percent_identity = 0.00

            # process cigar string

            # parse matches
            matches_list = re.findall(r'\d+(?=M)', cigar)
            match_count = 0

            for m in matches_list:
                m = int(m)
                match_count = match_count + m

            # parse mismatches
            mismatches_list = re.findall(r'\d+(?=D|I|N|S)', cigar)
            mismatch_count = 0.0

            for mm in mismatches_list:
                mm = int(mm)
                mismatch_count = mismatch_count + mm


            # determine percent identity
            percent_identity = 100 * (match_count / (match_count + mismatch_count))


            # determine alignment length
            align_len = match_count - mismatch_count


            if (align_len > 30):
                align_read = [read_name, flag, ref_name, pos, mapq, cigar, segment_seq, align_len, percent_identity]
                all_reads.append(align_read)



    alignment_info = [tool_info]


    # Prints output to CSV
    with open("alignment_details.csv", "w", newline="\n") as file:
        writer = csv.writer(file)
        writer.writerows(all_reads)

    with open("alignment_info.csv", "w", newline="\n") as file:
        writer = csv.writer(file)
        writer.writerows(accessions)
        writer.writerows(alignment_info)



    # Closes files
    inp.close()


if __name__ == '__main__':
    main()
