from logging         import info
from pandas          import DataFrame
from multiprocessing import Pool
from numpy           import sum

from .utrpy_argumentparser import UTRpyArgparser
from .utrpy_utr_add        import add_utrs
from .utrpy_io             import load_data, write_data
from .utrpy_logging        import logging_setup
from .utrpy_exon_extend    import exon_extend_threaded
from .utrpy_scaffold_split import scaffold_split

def main():

    args  = UTRpyArgparser().parse_args()

    logging_setup("info", f"{args.gff_utrpy}.log")

    info("Reading ...")
    gff_prediction, gff_assembly = load_data(args)

    info("Scaffold splitting ...")
    scaffolds      = list(gff_prediction[0].unique())
    gff_prediction = scaffold_split(gff_prediction)
    gff_assembly   = scaffold_split(gff_assembly)

    info("Extending exons ...")
    mp_args = zip([gff_prediction[s] for s in scaffolds],
                  [gff_assembly[s] if s in gff_assembly.keys() else DataFrame() for s in scaffolds],
                  [args.strict_strandedness for _ in scaffolds],
                  [args.maximum_exon_length for _ in scaffolds])
    
    with Pool(processes=args.threads) as pool:
        results = pool.map(exon_extend_threaded, mp_args)
        gffs    = [r[0] for r in results]

    info(f"Done all  (Extended exons: {sum([r[1] for r in results])})")

    if args.explicit:
        info("Adding UTRs")
        with Pool(processes=args.threads) as pool:
            results = pool.map(add_utrs, gffs)
            gffs    = [r[0] for r in results]
            info(f"Done all  (Annotated UTRs: {sum([r[1] for r in results])})")

    info("Writing ...")
    write_data(args, gffs, scaffolds)

    info("++++++++++++++++++++++++++++++++++")
    info("Simon says: thanks for using UTRpy")
    info("++++++++++++++++++++++++++++++++++")

if __name__ == '__main__':
    main()
