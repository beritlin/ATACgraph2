import subprocess
import pandas as pd 
import time 
import subprocess
import argparse
import csv
import numpy as np
import os

ATACgraph_dir = os.path.dirname(os.path.abspath(__file__))

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-promoter", help="promoter region from gene transcript start site (TSS), default 2000", type=int, default=2000)
    parser.add_argument('input_gtf')
    parser.add_argument('gtf_name')
    args = parser.parse_args()

    input_gtf_file = args.input_gtf
    outname = args.gtf_name

    # annotation name
    gene = "gene_body"
    exon = "exons"
    intron = "introns"
    utr3 = "3utr"
    utr5 = "5utr"
    cds = "cds"
    promoter = "gene_promoter"
    igr = "gene_igr"
    bed12 = [exon, intron, utr5, cds, utr3]
    annotation_name = [promoter, gene, exon, intron, utr5, cds, utr3, igr]

    # Extract UTR, exon, cds from gene
    filename = os.path.join(ATACgraph_dir, "extract_transcript_regions.py")
    subprocess.call(f"python {filename} -i {input_gtf_file} -o {outname} --gtf", shell=True)
    subprocess.call("rm *_coding*.bed *noncoding*.bed *5utr_start.bed", shell=True)

    # Convert this blockbed (bed12) to bed6
    for i in bed12:
        subprocess.call(f"cat {outname}_{i}.bed | bed12ToBed6 -i stdin -n > {outname}_{i}_bed6.bed", shell=True)

    # Find gene_body.bed
    genes = pd.read_csv(input_gtf_file, header=None, sep="\t", dtype={0: str})
    genes.columns = ['chr', 'unknow', 'exon', 'g_str', 'g_end', 'g_score', 'g_dir', '.', 'gene_name']
    genes = genes[genes.exon == 'exon']

    genes.gene_name = ' ' + genes.gene_name

    gene_col = genes['gene_name'].str.split(';', expand=True)
    gene_col.columns = gene_col.iloc[1, :]
    gene_id = gene_col.filter(regex='gene_id')
    gene_id = gene_id.iloc[:, 0].str.split(' ', expand=True)
    gene_id[2] = gene_id[2].map(lambda x: x.lstrip('"').rstrip('"'))
    gene_id.columns = ['num', 'g_name', 'gene_id']
    gene_bed = genes[['chr', 'g_str', 'g_end', 'g_score', 'g_dir']].join(gene_id['gene_id'])
    gene_bed = gene_bed[['chr', 'g_str', 'g_end', 'gene_id', 'g_score', 'g_dir']]
    gene_bed = gene_bed.drop_duplicates(subset=['g_str', 'g_end'], keep='first')
    gene_bed = gene_bed.sort_values(['chr', 'g_str'], ascending=[True, True])
    gene_bed = gene_bed.drop_duplicates(subset=['g_str'], keep='last')

    gene_group = gene_bed.groupby(['chr', 'gene_id', 'g_score', 'g_dir']).agg({'g_str': 'min', 'g_end': 'max'}).reset_index()
    gene_group = gene_group.drop_duplicates(subset=['g_str', 'g_end'], keep='first')

    gene_body = gene_group.sort_values(['chr', 'g_str'], ascending=[True, True])
    gene_body = gene_body.drop_duplicates(subset=['g_str'], keep='last')
    gene_body = gene_body[['chr', 'g_str', 'g_end', 'gene_id', 'g_score', 'g_dir']]

    gene_body.to_csv(f"{outname}_gene_body_bed6.bed", sep='\t', index=False, header=None)

    # Calculate genome size
    genome = gene_body.groupby(['chr']).agg({'g_end': 'max'}).reset_index()
    genome.to_csv(f"{outname}_genome.bed", sep="\t", index=False, header=None)
    genome_bp = genome['g_end'].sum()
    genome_bp = genome_bp * 1.0

    # Find promoter.bed
    gene_body['pro_str'] = np.where(gene_body.g_dir == '+', gene_body.g_str - args.promoter, gene_body.g_end - 0)
    gene_body['pro_end'] = np.where(gene_body.g_dir == '+', gene_body.g_str + 0, gene_body.g_end + args.promoter)
    num = gene_body._get_numeric_data()
    num[num < 0] = 0
    gene_promoter = gene_body[['chr', 'pro_str', 'pro_end', 'gene_id', 'g_score', 'g_dir']]
    gene_promoter.columns = ['chr', 'g_str', 'g_end', 'gene_id', 'g_score', 'g_dir']

    gene_promoter.to_csv(f"{outname}_gene_promoter_bed6.bed", sep='\t', index=False, header=None)

    # Find iGR bed
    gbname = [promoter, gene, exon, intron, utr5, cds, utr3]

    for i in gbname:
        subprocess.call(f"bedtools sort -i {outname}_{i}_bed6.bed|bedtools merge -c 4,5,6 -o collapse,collapse,collapse > {outname}_{i}_merge.bed", shell=True)

    subprocess.call(f"bedtools complement -i {outname}_gene_body_merge.bed -g {outname}_genome.bed > {outname}_gene_igr.bed", shell=True)

    geneigr = pd.read_csv(f"{outname}_gene_igr.bed", header=None, sep="\t")
    geneigr['.'], geneigr['...'], geneigr['..'] = ['.', '.', '.']
    geneigr.to_csv(f"{outname}_gene_igr_bed6.bed", header=None, index=False, sep="\t")

    subprocess.call(f"bedtools sort -i {outname}_gene_igr_bed6.bed|bedtools merge -c 4,5,6 -o collapse,collapse,collapse >{outname}_gene_igr_merge.bed", shell=True)

    subprocess.call(f"rm *_bed6.bed |rm *cds_bed6.bed | rm *5utr_bed6.bed | rm *exons_bed6.bed | rm *gene_igr_bed6.bed | rm *introns_bed6.bed |rm *3utr.bed|rm *5utr.bed |rm *_cds.bed|rm *_exons.bed|rm *_igr.bed|rm *_introns.bed ", shell=True)

if __name__ == '__main__':
    main()
