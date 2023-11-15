#!/usr/bin/env python3

import pysam
import sys
import subprocess
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
import scipy.fftpack
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam')
    parser.add_argument('fragment_distribution_outname')
    parser.add_argument('fragment_fft_outname')
    args = parser.parse_args()
    input_bam = args.input_bam
    output1 = args.fragment_distribution_outname
    output2 = args.fragment_fft_outname
    fragment_distribution(input_bam, output1, output2)

def fragment_distribution(input_bam, output1, output2):
    # Read bam file
    bam = pysam.AlignmentFile(input_bam)
    
    # Get fragment lengths
    fragment_len = np.array([b.template_length for b in bam if 0 < b.template_length < 2000], dtype=int)

    # Making histogram
    plt.style.use('classic')
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))
    n, bins, patches = ax1.hist(fragment_len, bins=100, color='gray', range=(0, 1000))
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel('Fragment length', fontsize=20)
    plt.ylabel('Fragment count', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax1.set_title('Fragment distribution', fontsize=20)
    outfile1 = output1
    plt.savefig(outfile1, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Calculating fft
    d2 = pd.concat([ffttable(n[15:96]),
                    ffttable(n[15:95]),
                    ffttable(n[15:94]),
                    ffttable(n[15:93]),
                    ffttable(n[15:92]),
                    ffttable(n[15:91]),
                    ffttable(n[15:90]),
                    ffttable(n[15:89]),
                    ffttable(n[15:88]),
                    ffttable(n[15:87])], ignore_index=True)
    
    d2 = d2.sort_values(by=['freq'], ascending=False)
    max_y_pos = 1 / d2.loc[d2.value.idxmax(), 'freq']
    max_y_pos_1 = round(max_y_pos,1)
    
    plt.style.use('classic')
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(1 / d2.freq, 10 * np.log10(d2.value + 1), color="royalblue")
    ax.axvline(max_y_pos, color="crimson", linestyle='--',linewidth=3) 
    ax.set_xlim(0, 400)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.set_xlabel('Period (bp)', fontsize=20)
    ax.set_ylabel('Power', fontsize=20)
    ymin, ymax = plt.ylim()
    xpeak = np.interp(max_y_pos, 1 / d2.freq, 10 * np.log10(d2.value + 1))
    xy = "(" + str(xpeak)[0:4] + "," + str(max_y_pos_1) + ")"
    plt.text(250, ymax - 1.5, xy, fontsize=20)
    ax.set_title('Period of fragment distribution', fontsize=20)
    outfile2 = output2
    plt.savefig(outfile2, dpi=300, bbox_inches='tight')

def ffttable(selected):
    selected = np.log2(selected + 1)
    selected2 = np.diff(selected)
    fragment_fft = sp.fftpack.fft(selected2)
    fragment_psd = np.abs(fragment_fft) ** 2
    fftfreq = sp.fftpack.fftfreq(len(fragment_psd), 10.0)
    d = pd.DataFrame({'freq': fftfreq, 'value': fragment_psd})
    d2 = d[d.freq > 0]
    return d2

if __name__ == '__main__':
    main()
