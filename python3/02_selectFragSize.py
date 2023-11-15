#!/usr/bin/env python
import pysam
import sys
import argparse

def open_bam_file(infile):
    try:
        samfile = pysam.AlignmentFile(infile, 'rb')
    except FileNotFoundError:
        sys.exit(f"The file {infile} does not exist")
    except Exception:
        sys.exit(f"The file {infile} does not have BAM format")

    try:
        samfile.check_index()
    except Exception:
        sys.exit(f"{infile} is not indexed. You MUST index the file first!")

    return samfile

def select_read(input_file, output_file, select_mode, select_size):
    samfile = open_bam_file(input_file)
    samout = pysam.AlignmentFile(output_file, "wb", template=samfile)

    for read in samfile:
        if (select_mode == 1 and abs(read.template_length) < select_size) or \
           (select_mode == 2 and abs(read.template_length) > select_size):
            samout.write(read)

    samout.close()
    samfile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam')
    parser.add_argument('output_bam')
    parser.add_argument('-f', '--filter', type=int, default=150, help='fragment length selection, default=150')
    parser.add_argument('-m', '--mode', type=int, default=1, help='Select fragments smaller [1] or larger [2] than filter size. Default=1')
    args = parser.parse_args()

    select_read(args.input_bam, args.output_bam, args.mode, args.filter)

if __name__ == '__main__':
    main()
