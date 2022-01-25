#!/usr/bin/env python3

"""
Kate Gerhardt (kgerhar4)
Fragment Recruitment Plotting Tool

Generates a fragment recruitment plot for an input SAM alignment file and cognate
reference genome.

Plot Details:
    Y-axis = percent identity
    X-axis = position in genome

Accepts:
    -Reference genome: FASTA formatted file (gzip compressed or uncompressed)
    -Alignment data: SAM

Assumptions:
    -Input files are in SAM format (for alignment data) and FASTA (for reference
        genome)
    -Unmapped reads have been removed from SAM file
    -


Run script like this:
python3 fragment_recruit.py -r ref_genome.fna -a alignment.sam

"""

import re
import os
import argparse
import parse_sam
import parse_fasta
import csv
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import numpy as np

def main():

    parser = argparse.ArgumentParser(description='Generates a fragment recruitment plot for alignment data')

    parser.add_argument('-r', '--ref_genome', type=str, required=True, help='Path to file representing an input reference genome')
    parser.add_argument('-a', '--align_data', type=str, required=True, help='Path to SAM formatted file containing alignment data')
    args = parser.parse_args()


    parse_fasta.main()
    parse_sam.main()

    # get reference data
    reference_data = []
    fasta_headers = []
    with open('reference_details.csv') as input:
        rows = csv.reader(input)
        reference_data = next(rows)


    total_ref_len = int(reference_data[-1])
    main_ref_header = reference_data[0]
    main_ref_chrom_len = int(reference_data[1])
    fasta_header = [main_ref_header, main_ref_chrom_len, total_ref_len]
    fasta_headers.append(fasta_header)

    ref_data_len = len(reference_data)
    i = 3
    j = 0

    while (i < (ref_data_len - 1)):
        fasta_head = []
        fasta_head.append(reference_data[i-1])
        fasta_head.append(int(reference_data[i]))
        fasta_headers.append(fasta_head)
        i+=2


    # get tool name and version
    alignment_header = []
    with open('alignment_info.csv') as input:
        rows = csv.reader(input)
        for row in rows:
            alignment_header.append(row)


    main_accession = alignment_header[0][0]


    # get x coordinate, get y coordinate, get percent identity
    alignment_data = []
    with open('alignment_details.csv') as input:
        rows = csv.reader(input)
        for row in rows:
            reference = row[2]
            x_coordinate = row[3]
            y_coordinate = row[8]

            data_point = [reference, x_coordinate, y_coordinate]
            alignment_data.append(data_point)

    # make plot
    # main_chr length
    x_coord = []
    y_coord = []
    # accs = []
    # ii = 0

    for data_point in alignment_data:
        if (data_point[0] == main_accession):
            x_coord.append(int(data_point[1]))
            y_coord.append(float(data_point[2]))
            # acc.append(data_point[0])
            # acc.append(main_accession)
        # else:
        #     ii += 1
        # accs.append(acc)

    # for x in accs:
    #     print(x)

    # print(ii)

    x = np.array(x_coord)
    y = np.array(y_coord)

    mplstyle.use('fast')


    start_tick = 1
    mid_tick = int(main_ref_chrom_len/2)
    quart_tick = int(mid_tick/2)
    three_quart_tick = int(quart_tick * 3)
    end_tick = main_ref_chrom_len

    tick_2 = int(quart_tick/2)
    tick_5 = int(mid_tick+tick_2)
    tick_7 = int(three_quart_tick+tick_2)
    tick_4 = int(quart_tick + tick_2)

    xticks = [start_tick, tick_2, quart_tick, tick_4, mid_tick, tick_5, three_quart_tick, tick_7, end_tick]
    xtick_labels = []

    for tick in xticks:
        xtick_labels.append(str(tick))


    yticks = [80, 85, 90, 95, 100]
    ytick_labels = []

    for tick in yticks:
        ytick_labels.append(str(tick))


    organism = re.search(r'(?<=>).+', main_ref_header)
    organism = organism.group()

    plot = plt.scatter(x,y, marker="_", s=.09, alpha=.8)
    plt.xlim(0, main_ref_chrom_len)
    plt.ylim(80,100.5)
    plt.rc('xtick', labelsize=8)
    plt.xticks(xticks, xtick_labels, fontsize=5)
    plt.yticks(yticks, ytick_labels, fontsize=6)
    plt.ylabel('Sequence Identity (%)', fontsize=8)
    plt.xlabel('Genome Position (bp)', fontsize=7)
    plt.suptitle("Feline Fecal Microbiome")
    plt.title(organism, fontsize=9)

    plot_file_name = organism + "_figure"

    plt.savefig(plot_file_name, format="png", dpi=600)



    # delete files that were made
        # delete alignment_details.csv
        # delete alignment_info.csv
        # delete reference_details.csv




if __name__ == '__main__':
    main()
