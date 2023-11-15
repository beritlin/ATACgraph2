#!/usr/bin/env python3

import pysam
import sys

def processInput(target, chrs):
    split = target.split(',')
    a = [x for x in chrs if x not in split]
    if any(eachTarget not in chrs for eachTarget in split):
        print("One or more chromosomes are not in the list. Please check!")
        sys.exit()
    return (a, split)

def openBamFile(infile):
    try:
        samfile = pysam.AlignmentFile(infile, 'rb')
    except IOError:
        sys.exit("The file {} does not exist".format(infile))
    except Exception as e:
        sys.exit("Error opening {}: {}".format(infile, e))

    try:
        samfile.check_index()
    except Exception as e:
        sys.exit("{} is not indexed. You MUST index the file first! Error: {}".format(infile, e))

    return samfile

def rmChr(infile, outfile, targetChr):
    samfile = openBamFile(infile)
    samout = pysam.AlignmentFile(outfile, "wb", template=samfile)
    chrs = samfile.references
    chrs_keep, chrs_rm = processInput(targetChr, chrs)

    keepread = 0
    dupread = 0

    for chrom in chrs_keep:
        fetchread = samfile.fetch(chrom)
        for read in fetchread:
            samout.write(read)
            keepread += 1
            if read.is_duplicate:
                dupread += 1

    totalread = sum(stats.total for stats in samfile.get_index_statistics())
    removeread = totalread - keepread
    removeratio = removeread / totalread if totalread != 0 else 0  # Avoid division by zero

    for chrom in samfile.get_index_statistics():
        if chrom.contig in chrs_rm:
            rmcount = chrom.total
            print("Removed {} {} reads".format(chrom.contig, rmcount))

    print("Removed total {} out of {} ({:.3%})".format(removeread, totalread, removeratio))

    samout.close()
    samfile.close()

def main():
    if len(sys.argv) != 4:
        print("Usage: python 00_rmChr <input.bam> <output.bam> <chrM>")
        print("")
        print("If you need to remove multiple chromosomes, use a comma to separate. For example chrM,chrPt")
        sys.exit()

    infile = sys.argv[1]
    outfile = sys.argv[2]
    targetChr = sys.argv[3]
    rmChr(infile, outfile, targetChr)

if __name__ == '__main__':
    main()
