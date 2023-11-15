#!/usr/bin/env python

import pysam
import pandas as pd
import numpy as np
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam')
    parser.add_argument('output_bed')
    parser.add_argument('-s', '--separate', type=int, default=150, help="border length of long and short fragment(bp), default: 150 ")
    parser.add_argument('-b', '--bin', type=int, default=10, help="bin size of a group on genome(bp), default: 10")
    parser.add_argument('-f', '--filter', type=int, default=2, help="minimal  fragment junction tracks, default:2")
    args = parser.parse_args()

    input_bam, outbed, selectSize, filterFunc, bin_size = args.input_bam, args.output_bed, args.separate, args.filter, args.bin
    leftshift, rightshift = 0, 0

    with open(outbed, 'w') as c:
        c.write(f'track name=junctions description="TopHat junctions"\n')

    bam = pysam.AlignmentFile(input_bam)
    data = [(bam.references[b.rname], b.pos + leftshift, b.pos + b.tlen - rightshift) for b in bam if 0 < b.tlen < 2000]
    c = pd.DataFrame(data, columns=['chr', 'str', 'end'])
    c[['str', 'end']] = c[['str', 'end']] // bin_size * bin_size
    c1_score = c.groupby(['chr', 'str', 'end']).size().reset_index(name='count2')
    c1_score = c1_score[c1_score['count2'] >= filterFunc]
    c1_score['location'] = c1_score['end'] - c1_score['str'] - 1
    c1_score['name'] = c1_score['chr'] + '_' + c1_score['str'].astype(str) + '_' + c1_score['end'].astype(str)
    c1_score['dir'] = np.where(c1_score['location'] > selectSize, '+', '-')
    c1_score[['thickstart', 'thickend']] = c1_score[['str', 'end']]
    c1_score['block_count'] = 2
    c1_score['block_size'] = "1,1"
    c1_score['rgb'] = np.where(c1_score['location'] > selectSize, '255,0,0', '0,0,255')
    c1_score['zero'] = 0
    c1_score['block_location'] = c1_score['zero'].astype(str) + "," + c1_score['location'].astype(str)
    c1_junction_pd = c1_score[['chr', 'str', 'end', 'name', 'count2', 'dir', 'thickstart', 'thickend', 'rgb', 'block_count',
                                'block_size', 'block_location']]
    c1_junction_pd.to_csv(outbed, mode='a', header=None, index=None, sep="\t")

if __name__ == '__main__':
    main()
