from logging         import info
from pandas          import DataFrame, read_csv
from multiprocessing import Pool
from numpy           import sum

from .utrpy_argumentparser import UTRpyArgparser
from .utrpy_exon_extend    import exon_extend_threaded, exon_extend
from .utrpy_logging        import logging_setup
from .utrpy_scaffold_split import scaffold_split

def main():

    args = UTRpyArgparser().parse_args()

    logging_setup("info", "UTRpy.log")

    info("Reading ...")
    gff_prediction = read_csv(args.gff_prediction, sep="\t", header=None, comment="#")
    gff_assembly   = read_csv(args.gff_assembly,   sep="\t", header=None, comment="#")

    info("Scaffold splitting ...")
    scaffolds      = list(gff_prediction[0].unique())
    gff_prediction = scaffold_split(gff_prediction)
    gff_assembly   = scaffold_split(gff_assembly)

    info("Exon extension ...")
    mp_args = zip([gff_prediction[s] for s in scaffolds],
                  [gff_assembly[s] if s in gff_assembly.keys() else DataFrame() for s in scaffolds],
                  [args.minimum_exon_overlap for _ in scaffolds],
                  [args.conservative for _ in scaffolds],
                  [args.maximum_exon_length for _ in scaffolds])
    
    with Pool(processes=args.threads) as pool:
        results = pool.map(exon_extend_threaded, mp_args)

    info(f"\tDone all  (Extended exons: {sum([r[1] for r in results])})")

    info("Writing ...")
    gff_utrpy = sorted([r[0] for r in results], key=lambda x: scaffolds.index(x.iloc[0][0]))
    for df in gff_utrpy:
        df.to_csv(args.gff_utrpy, sep="\t", mode="a", index=False, header=None)

    info("++++++++++++++++++++++++++++++++++")
    info("Simon says: thanks for using UTRpy")
    info("++++++++++++++++++++++++++++++++++")

if __name__ == '__main__':
    main()