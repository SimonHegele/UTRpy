"""
Functions for file reading and writing
"""

from pandas import read_csv

def load_data(args):

    gff_prediction = read_csv(args.gff_prediction, sep="\t", header=None, comment="#")
    gff_assembly   = read_csv(args.gff_assembly,   sep="\t", header=None, comment="#")

    return gff_prediction, gff_assembly

def write_data(args, gffs, scaffolds):
    
    gff_utrpy = sorted(gffs, key=lambda x: scaffolds.index(x.iloc[0][0]))
    for df in gff_utrpy:
        df.to_csv(args.gff_utrpy, sep="\t", mode="a", index=False, header=None)
