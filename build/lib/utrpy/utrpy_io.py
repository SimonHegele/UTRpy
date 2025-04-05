"""
Script Name:    utrpy_io.py
Description:    File reading / Writing
Author:         Simon Hegele
Date:           2025-04-01
Version:        0.2
License:        GPL-3
"""

from logging import info

from .utrpy_gff_utils import load_gff, seqname_split_sort, get_seqnames

def load_data(args):

    info("Reading ...")
    gff_prediction = load_gff(args.gff_prediction)
    gff_assembly   = load_gff(args.gff_assembly)

    seqnames       = get_seqnames(gff_prediction)

    gff_prediction = seqname_split_sort(gff_prediction, seqnames)
    gff_assembly   = seqname_split_sort(gff_assembly,   seqnames)

    return gff_prediction, gff_assembly, seqnames

def write_data(args, gffs, scaffolds):

    info("Writing ...")
    gff_utrpy = sorted(gffs, key=lambda x: scaffolds.index(x.iloc[0,0]))
    for df in gff_utrpy:
        df.to_csv(args.gff_utrpy, sep="\t", mode="a", index=False, header=None)