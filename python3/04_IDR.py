### need to be fixed
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_peakAs', help='Peaks files, separated by commas')
    parser.add_argument('input_peakBs', help='Peaks files, separated by commas')
    parser.add_argument('output')
    parser.add_argument('-p', help='path to idr', default=" ",dest='path')

    args = parser.parse_args()

    groupPeakA = args.input_peakAs
    groupPeakB = args.input_peakBs
    outBed = args.output
    path = args.path

    # idrAFile, idrBFile = "IDR_A.txt", "IDR_B.txt"

    subprocess.call(f"{path}idr --sample {groupPeakA} {groupPeakB} --input-file-type narrowPeak --rank p.value  --output-file {outBed}  --plot --log-output-file  {outBed}.log", shell=True)
    # subprocess.call(f"bedtools subtract -a {idrAFile} -b {idrBFile} -A > {outBed}_A", shell=True)
    # subprocess.call(f"bedtools subtract -a {idrBFile} -b {idrAFile} -A > {outBed}_B", shell=True)

if __name__ == '__main__':
    main()

