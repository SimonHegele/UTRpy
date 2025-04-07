"""
Module Name:    utrpy_exon_extend.py
Description:    Extending exons (and their corresponding genes and transcripts) using
                suitable exons from the gene prediction to obtain UTRs
Author:         Simon Hegele
Date:           2025-04-01
Version:        0.2
License:        GPL-3
"""

from logging            import getLogger, info
from multiprocessing    import Pool
from pandas             import DataFrame, Series
from typing             import Generator

from .utrpy_check       import check
from .utrpy_gff_utils   import next_feature_index, get_feature_ancestor
from .utrpy_logging     import logging_setup

def overlapping_ta_exons(gff_ta: DataFrame, gp_exon: Series) -> DataFrame:
    
    mask1 = ((gp_exon["start"] <= gff_ta["end"])  & (gp_exon["end"] >= gff_ta["start"]))
    mask2 = ((gff_ta["start"]  <= gp_exon["end"]) & (gff_ta["end"]  >= gp_exon["start"]))

    return gff_ta.loc[(mask1 | mask2) & (gff_ta["feature"] == "exon")]

def select_ta_exon(ta_exons: DataFrame, ambiguities: str) -> Series | None:

    if ta_exons.shape[0] == 0:
        return None
    
    lengths = (ta_exons["end"]-ta_exons["start"])
    
    if len(lengths.unique()) == 1:
        return ta_exons.iloc[0]
    
    match ambiguities:
        case "none":
            return None
        case "smallest":
            return ta_exons.loc[lengths.idxmin()]
        case "longest":
            return ta_exons.loc[lengths.idxmax()]
     
def exon_matches(gff_gp: DataFrame,
                 gff_ta: DataFrame,
                 known_strand: bool,
                 max_exon_length: int,
                 ambiguities: str) -> Generator:

    i = next_feature_index(gff_gp, -1, "exon")

    while not i is None:

        gp_exon  = gff_gp.iloc[i]
        gp_tran  = get_feature_ancestor(gff_gp, gp_exon, "transcript")
        ta_exons = overlapping_ta_exons(gff_ta, gp_exon)
        ta_exons = ta_exons.loc[ta_exons.apply(lambda ta_exon: check(ta_exon,
                                                                     gp_exon,
                                                                     gp_tran,
                                                                     known_strand,
                                                                     max_exon_length),
                                               axis = 1)]
        ta_exon  = select_ta_exon(ta_exons, ambiguities)

        if not ta_exon is None:
            yield gp_exon, ta_exon

        i = next_feature_index(gff_gp, i, "exon")

def update(gff: DataFrame, gp_exon: Series, ta_exon: Series, type: str):
    
    ancestor      = get_feature_ancestor(gff, gp_exon, type)
    i             = (gff == ancestor).all(axis=1).idxmax()

    if not "(UTRpy)" in gff.iloc[i,1]:
        gff.iloc[i,1] = f"{ancestor["source"]} + {ta_exon["source"]} (UTRpy)"
    gff.iloc[i,3] = min(gff.iloc[i, 3], ta_exon["start"])
    gff.iloc[i,4] = max(gff.iloc[i, 4], ta_exon["end"])

def exon_extend(gff_gp: DataFrame,
                gff_ta: DataFrame,
                known_strand: bool,
                max_ex_len: int,
                log_file: str,
                ambiguities: str) -> tuple[DataFrame,int]:

    n = 0

    logging_setup("info", log_file)
    logger = getLogger()

    if gff_gp.shape[0]==0 or gff_ta.shape[0]==0:
        return gff_gp, 0
    
    for gp_exon, ta_exon in exon_matches(gff_gp, gff_ta, known_strand, max_ex_len, ambiguities):

        for type in ["exon", "transcript", "gene"]:
            update(gff_gp, gp_exon, ta_exon, type)
        
        n += 1

    logger.info("\tDone {0:<20} (Extended exons: {1:>5})".format(gff_gp.iloc[0,0],n))
    return gff_gp, n

def exon_extend_threaded(args: tuple[DataFrame,DataFrame,int,bool,int]) -> tuple[DataFrame,int]:

    gff_gp, gff_ta, known_strand, max_ex_len, log_file, ambiguities = args

    return exon_extend(gff_gp, gff_ta, known_strand, max_ex_len, log_file, ambiguities)

def exon_extend_multithreaded(gff_prediction: DataFrame,
                              gff_assembly: DataFrame,
                              seqnames: list[str],
                              args) -> list[DataFrame]:
    
    info("Extending exons ...")
    mp_args = zip([gff_prediction[s]        for s in seqnames],
                  [gff_assembly[s]          for s in seqnames],
                  [args.know_strand         for _ in seqnames],
                  [args.max_ex_len          for _ in seqnames],
                  [f"{args.gff_utrpy}.log"  for _ in seqnames],
                  [args.ambiguities         for _ in seqnames],)
    
    with Pool(processes=args.threads) as pool:
        results = pool.map(exon_extend_threaded, mp_args)
        gffs    = [r[0] for r in results]

    info(f"Done all  (Extended exons: {sum([r[1] for r in results])})")
    return gffs
