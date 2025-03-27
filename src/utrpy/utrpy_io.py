"""
Functions for file reading and writing
"""

from logging import info
from pandas  import read_csv

def load_data(args):

    info("Reading ...")
    gff_prediction = read_csv(args.gff_prediction, sep="\t", header=None, comment="#")
    gff_assembly   = read_csv(args.gff_assembly,   sep="\t", header=None, comment="#")

    return gff_prediction, gff_assembly

def write_data(args, gffs, scaffolds):

    info("Writing ...")
    gff_utrpy = sorted(gffs, key=lambda x: scaffolds.index(x.iloc[0][0]))
    for df in gff_utrpy:
        df.to_csv(args.gff_utrpy, sep="\t", mode="a", index=False, header=None)
