#!/usr/bin/env python
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Peaks files, separated by commas')
    parser.add_argument('output')
    parser.add_argument('bam')

    args = parser.parse_args()

    peak = args.input
    outname = args.output
    bam = args.bam

    subprocess.call(f'''awk 'OFS="\t" {{print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}}' {peak} > {outname}.SAF''', shell=True)
    subprocess.call(f"featureCounts -p -a {outname}.SAF -F SAF -o {outname}-frip.txt {bam} >>{outname}-frip.log 2>&1", shell=True)

if __name__ == '__main__':
    main()
