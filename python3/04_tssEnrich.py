#!/usr/bin/env python
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Peaks files, separated by commas')
    parser.add_argument('gene', help='gtf name')
    parser.add_argument('output')
    parser.add_argument('bam')

    args = parser.parse_args()

    bw = args.input
    gene = args.gene
    outname = args.output
    bam = args.bam

    subprocess.call(f"computeMatrix reference-point  --referencePoint TSS -p 15 -b 1000 -a 1000 -S {bw} -R {gene}_gene_body_bed6.bed -o {outname}-TSS.gz --outFileSortedRegions {outname}-TSS.bed", shell=True)
    subprocess.call(f'''awk 'OFS="\t" {{print $4, $1, $2+1, $3, "."}}' {outname}-TSS.bed > {outname}-TSS_tmp.SAF''', shell=True)
    subprocess.call(f'''awk '(NR>1)' {outname}-TSS_tmp.SAF > {outname}-TSS.SAF''', shell=True)
    subprocess.call(f"rm {outname}-TSS_tmp.SAF ", shell=True)
    subprocess.call(f"featureCounts -T 10 -p -a {outname}-TSS.SAF  -F SAF -o {outname}-TSS.txt {bam} >>{outname}-TSS.log 2>&1", shell=True)

if __name__ == '__main__':
    main()

