import sys
import pandas as pd
import os
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='separate', help="1: integration site; 2: full-extend fragment", type=int, default=2)
    parser.add_argument("-shift", help="shift size from integration site(bp), default: 50", dest='shift', default=50, type=int)
    parser.add_argument("-ES", help="extend size from integration site (bp), default: 100", dest='extend', default=100, type=int)
    parser.add_argument("-bs", help="bin size for bigwig (bp), default: 10", dest='binsize', default=10, type=int)
    parser.add_argument("-c", dest="control_bam", help='input control bam file, default: none', default='', type=str)
    parser.add_argument('-p', help='the path to macs3', default=" ", dest='path', type=str)
    parser.add_argument('input_bam', help='input ATAC-seq bam file', type=str)
    parser.add_argument('output_name', help='name for output files', type=str)
    parser.add_argument('gene_bed', help='Gene or promoter annotation bed file, either gene_body_bed6.bed or gene_promoter_bed6.bed', type=str)

    args = parser.parse_args()

    input_bam = args.input_bam
    Outname = args.output_name
    control_bam = args.control_bam
    control_para = '' if not control_bam else '-c ' + control_bam
    shiftSize = args.shift
    extSize = args.extend
    bs = args.binsize
    gene = args.gene_bed
    path = args.path

    subprocess.call(f'''samtools index {input_bam}''', shell=True)

    if args.separate == 1:
        subprocess.call(f'''{path}macs3 callpeak {control_para} -t {input_bam} --nomodel --shift -{shiftSize} --extsize {extSize} -n {Outname} 2>&1''', shell=True)
        subprocess.call(f'''bamCoverage -b {input_bam} -bs {bs} --normalizeUsing RPKM --Offset 1 20 -o {Outname}_coverage.bw 2>&1''', shell=True)
        subprocess.call(f'''bedtools intersect -wo -a {gene} -b {Outname}_peaks.narrowPeak > {Outname}_peak_gene_list.txt''', shell=True)
    elif args.separate == 2:
        subprocess.call(f'''{path}macs3 callpeak {control_para} -t {input_bam} --format BAMPE -n {Outname} 2>&1''', shell=True)
        subprocess.call(f'''bamCoverage -b {input_bam} -bs {bs} --normalizeUsing RPKM -e -o {Outname}_coverage.bw 2>&1''', shell=True)
        subprocess.call(f'''bedtools intersect -wo -a {gene} -b {Outname}_peaks.narrowPeak > {Outname}_peak_gene_list.txt''', shell=True)
    else:
        print("Please enter separation number, 1:")
        subprocess.call(f'''bedtools intersect -wo -a {gene} -b {Outname}_peaks.narrowPeak > {Outname}_peak_gene_list.txt''', shell=True)

        g1 = pd.read_csv(f"{Outname}_peak_gene_list.txt", sep="\t", header=None)
        g1.columns = ['chr', 'sta', 'end', 'gene', 'score', 'dir', 'chrp', 'sta_p', 'end_p', 'peakname', 'score_p', 'strand', 'signalValue', 'pval', 'qValue', 'peak', 'overlap']

        g1['Peak'] = "(" + g1['chrp'] + ":" + g1['sta_p'].astype(str) + "-" + g1['end_p'].astype(str) + ')'

        g2 = g1[['chr', 'sta', 'end', 'gene', 'Peak']]
        g2.to_csv(f"{Outname}_peak_gene_list.txt", sep="\t", index=None)


if __name__ == '__main__':
    main()
