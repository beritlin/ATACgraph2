#!/usr/bin/env python
# 2019.12.9 added FC and pvalue
# 2029.03.09 added parameter for FC
import pybedtools
import argparse
import csv
import warnings
import pyBigWig
import numpy as np
from scipy.stats import ttest_ind

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_peakAs', help='Peaks files, separate by comma')
    parser.add_argument('input_peakBs', help='Peaks files, separate by comma')
    parser.add_argument('input_bwAs', help='BigWig files, separate by comma')
    parser.add_argument('input_bwBs', help='BigWig files, separate by comma')
    parser.add_argument('output_bedA')
    parser.add_argument('output_bedB')
    parser.add_argument('-c', '--fold_change', type=float, default=2.0,
                        help='Fold change cutoff. Default: 2')
    parser.add_argument('-p', '--p_value', type=float, default=0.05,
                        help='P-value cutoff. Default: 0.05')
    args = parser.parse_args()

    inPeakA, inPeakB, inBwA, inBwB = map(lambda x: x.split(','), [args.input_peakAs, args.input_peakBs, args.input_bwAs, args.input_bwBs])
    fold_change, p_value = args.fold_change, args.p_value

    groupBw = inBwA + inBwB
    outBedA, outBedB = args.output_bedA, args.output_bedB

    x = pybedtools.BedTool()
    result = x.multi_intersect(i=inPeakA + inPeakB, cluster='T')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        df = result.to_dataframe()

    df = df[df.name >= min(len(inPeakA), len(inPeakB))]

    for bwFile in groupBw:
        bw = pyBigWig.open(bwFile)
        df[bwFile] = [np.mean(bw.values(peak['chrom'], peak['start'] - 1, peak['end'])) for _, peak in df.iterrows()]

    df = df[df[groupBw].min(axis=1) > 1]
    df.loc[:, 'fc'] = np.log2(df[inBwA].mean(axis=1) / df[inBwB].mean(axis=1))
    df.loc[:, 'pvalue'] = np.where(
        min(len(inPeakA), len(inPeakB)) >= 2,
        [ttest_ind(
            peak[inBwA].astype(float),
            peak[inBwB].astype(float),
            equal_var=False
        ).pvalue for _, peak in df.iterrows()], 0
    )

    peaksA = df[(df.fc >= np.log2(fold_change)) & (df.pvalue < p_value)].copy()
    peaksA.loc[:, 'id'] = [f'PeakA_{i}' for i in range(1, len(peaksA) + 1)]
    peaksA.loc[:, 'fc2'] = peaksA.fc.map("{:.2f}".format)
    peaksA.loc[:, 'pvalue2'] = peaksA.pvalue.map("{:.2f}".format)
    peaksA = peaksA[['chrom', 'start', 'end', 'id',  'pvalue2', 'fc2']]

    peaksB = df[(df.fc <= -np.log2(fold_change)) & (df.pvalue < p_value)].copy()
    peaksB.loc[:, 'id'] = [f'PeakB_{i}' for i in range(1, len(peaksB) + 1)]
    peaksB.loc[:, 'fc2'] = peaksB.fc.map("{:.2f}".format)
    peaksB.loc[:, 'pvalue2'] = peaksB.pvalue.map("{:.2f}".format)
    peaksB = peaksB[['chrom', 'start', 'end', 'id', 'pvalue2','fc2']]

    peaksA.to_csv(outBedA, sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
    peaksB.to_csv(outBedB, sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)

    dar =  df[(df.fc <= -np.log2(fold_change)) or (df.fc <= np.log2(fold_change)) & (df.pvalue < p_value)].copy()
    dar.loc[:, 'id'] = [f'dar_{i}' for i in range(1, len(dar) + 1)]
    dar.loc[:, 'fc2'] = dar.fc.map("{:.2f}".format)
    dar.loc[:, 'pvalue2'] = dar.pvalue.map("{:.2f}".format)
    dar = dar[['chrom', 'start', 'end', 'id',  'pvalue2', 'fc2']]

    dar.to_csv("dar.txt", sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
