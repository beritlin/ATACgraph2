#!/usr/bin/env python

import pandas as pd
import subprocess
import math
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

def genome_size(samplename):
    return pd.read_csv(f"{samplename}_genome.bed", sep='\t', header=None)[1].sum() * 1.0

def peak_bp(peak):
    peak = pd.read_csv("/work1/home/peiyu/2023_atac/demo_rmM_peakcall_peaks.narrowPeak", header=None, sep="\t").iloc[:, 0:6].drop_duplicates(subset=[1, 2], keep='first')
    peak.columns = ['chr', 'peak_str', 'peak_end', 'peak_id', 'peak_value', 'peak_dir']
    peak['bp'] = (peak['peak_end'] - peak['peak_str'])
    return peak['bp'].sum()

def ano_bp(anno, samplename):
    return pd.read_csv(f"{samplename}_{anno}_merge.bed", header=None, sep="\t").drop_duplicates(subset=[1, 2], keep='first')[2].sum()

def annopeakbp(peak, anno):
    return float(pd.read_csv(f"{peak}_{anno}.txt", header=None, sep="\t").iloc[:, -1].sum())

def enrichment_num(peak, anno, genome_bp, samplename):
    try:
        return math.log((annopeakbp(peak, anno) / peak_bp(peak)) / (ano_bp(anno, samplename) / genome_bp), 2)
    except pd.errors.EmptyDataError:
        return 0

def associate(input_peak, samplename, promoter):
    peak = pd.read_csv(input_peak, header=None, sep="\t").iloc[:, 0:6]
    gene_body = pd.read_csv(f"{samplename}_gene_body_merge.bed", header=None, sep="\t")
    peak.columns = ['chr', 'peak_str', 'peak_end', 'peak_name', 'peak_value', 'peak_dir']
    gene_body.columns = ['chr', 'gbed_str', 'gbed_end', 'gene_id', 'gene_value', "gene_dir"]
    gene_body['pro_str'] = np.where(gene_body['gene_dir'] == '+', gene_body['gbed_str'] - promoter, gene_body['gbed_end'] - 0)
    gene_body['pro_end'] = np.where(gene_body['gene_dir'] == '+', gene_body['gbed_str'] + 0, gene_body['gbed_end'] + promoter)
    combined = pd.merge(gene_body, peak, on='chr')
    combined['Genebody'] = np.where((combined['gbed_str'] < combined['peak_end']) & (combined['gbed_end'] > combined['peak_str']), 1, 0)
    combined['Promoter'] = np.where((combined['pro_str'] < combined['peak_end']) & (combined['pro_end'] > combined['peak_str']), 1, 0)
    summary = combined[(combined['Promoter'] > 0) | (combined['Genebody'] > 0)]
    s1_group = summary.groupby(['gene_id', 'chr', 'gbed_str', 'gbed_end', 'gene_dir']).agg({'Genebody': 'sum', 'Promoter': 'sum'}).reset_index().sort_values(["chr", "gbed_str"])
    s1_group.to_csv(f"{input_peak}_gene_summary_table.xls", sep='\t', index=False)

def coverage_heatmap(coverage, samplename, input_peak):
    subprocess.call(
        f"computeMatrix scale-regions -S {coverage} -R {samplename}_gene_body_merge.bed --missingDataAsZero -bs 10 -a 1000 -b 1000 -out {coverage}gene_body.matrix.gz --outFileNameMatrix {coverage}gene_body.matrix.txt 2>&1",
        shell=True)
    subprocess.call(
        f"plotHeatmap -m {coverage}gene_body.matrix.gz -out {coverage}_gene_body_heatmap.pdf --legendLocation none",
        shell=True)
    subprocess.call(
        f"computeMatrix reference-point --referencePoint center -S {coverage} -R {input_peak} --missingDataAsZero -bs 10 -a 1000 -b 1000 -out {coverage}_peak.matrix.gz --outFileNameMatrix {coverage}_peak.matrix.txt 2>&1",
        shell=True)
    subprocess.call(
        f"plotHeatmap --refPointLabel peak  --regionsLabel peaks --xAxisLabel 'peak distance(bp)' -m {coverage}_peak.matrix.gz -out {coverage}_Peak_heatmap.pdf --legendLocation none",
        shell=True)

def annotation_peak(input_peak, samplename):
    annotation_name = ["gene_promoter", "gene_body", "exons", "introns", "5utr", "cds", "3utr", "gene_igr"]
    for anno in annotation_name:
        subprocess.call(
            f"bedtools intersect -nonamecheck -a {samplename}_{anno}_merge.bed -b {input_peak} -wo > {input_peak}_{anno}.txt",
            shell=True)

def enrichment(input_peak, genome_bp, samplename):
    annotation_name = ["gene_promoter", "gene_body", "exons", "introns", "5utr", "cds", "3utr", "gene_igr"]
    annotationname = ['Promoter', 'Genebody', 'Exon', 'Intron', '5UTR', 'CDS', '3UTR', 'IGR']
    enrichment_data = [enrichment_num(input_peak, anno, genome_bp, samplename) for anno in annotation_name]

    fig, ax1 = plt.subplots()
    ax1.bar(range(len(annotationname)), enrichment_data, align='center', color="gray")
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    # plt.ylim(ymin=-13, ymax=6)
    plt.xticks(range(len(annotationname)), annotationname, fontsize='small')
    plt.ylabel("Fold Enrichment (log2)")
    plt.savefig(f"{input_peak}_Fold_Enrichment.pdf", dpi=400, bbox_inches='tight')
    plt.close()

    fe_table = pd.DataFrame([enrichment_data], columns=annotationname)
    fe_table.to_csv(f"{input_peak}_Fold_Enrichment_Table.txt", index=None, sep="\t")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_peak', type=str)
    parser.add_argument('input_bigwig', type=str)
    parser.add_argument('gtf_name', type=str)
    parser.add_argument('-p', '--promoter', type=int, default=2000)
    args = parser.parse_args()

    promoter = args.promoter
    input_peak = args.input_peak
    samplename = args.gtf_name
    bigwig = args.input_bigwig

    associate(input_peak, samplename, promoter)
    coverage_heatmap(bigwig, samplename, input_peak)
    annotation_peak(input_peak, samplename)
    genome_bp = genome_size(samplename)
    enrichment(input_peak, genome_bp, samplename)

if __name__ == '__main__':
    main()
