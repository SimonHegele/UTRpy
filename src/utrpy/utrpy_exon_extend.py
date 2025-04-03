"""
Script Name:    utrpy_exon_extend.py
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

def overlapping_ta_exons(gff_ta: DataFrame,
                         gp_exon: Series) -> DataFrame:
    """
    Args:
        gff_ta (DataFrame): Part of the transcriptome assembly GFF with features of
                            the same sequence
        gp_exon (Series):   An exon from the gene prediction of the same sequence

    Returns:
        DataFrame: gff_ta filtered for exons overlapping with gp_exon
    """
    
    # Case 1:
    # -> gp_exon["start"] <= ta_exon["end"] and gp_exon["end"] >= ta_exon["start"]
    # gp_exon starts at position before ta_exon ends and also ends before ta_exon ends
    # Examples:
    #   gp_exon -----       -----     ---
    #                   or        or 
    #   ta_exon   -----     -----    -----
    mask1 = ((gp_exon["start"] <= gff_ta["end"]) & (gp_exon["end"] >= gff_ta["start"]))

    # Case 2:
    # -> ta_exon["start"] <= gp_exon["end"] and ta_exon["end"] >= gp_exon["start"]
    # Examples:
    #   gp_exon   -----    -----     -----
    #                   or        or 
    #   ta_exon -----      -----      ---
    mask2 = ((gff_ta["start"] <= gp_exon["end"]) & (gff_ta["end"] >= gp_exon["start"]))

    return gff_ta.loc[(mask1 | mask2) & (gff_ta["feature"] == "exon")]
     
def exon_matches_sensitive(gff_gp: DataFrame,
                           gff_ta: DataFrame,
                           known_strand: bool,
                           max_exon_length: int) -> Generator:
    """
    Retrieving indices of matching exons from the gene prediction and transcriptome
    assembly.

    Args:
        gff_gp (DataFrame):     Part of the gene prediction GFF with features of the same
                                sequence
        gff_ta (DataFrame):     Part of the transcriptome assembly GFF with features of
                                the same sequence
        known_strand (bool):    Only use exons if the strands of both exons are known.
        max_exon_length (int):  Length limit for exons from the transcriptome assembly

    Yields:
        Generator: 2-tuples (i,j) with i referring to an exon gp_exon from the gene
                   gene prediction and j referring to an exon ta_exon from the
                   transcriptome assembly assembly where ta_exon is an UTR extension of
                   gp_exon.
    """
    i = next_feature_index(gff_gp, -1, "exon")
    while i != None:
        gp_exon  = gff_gp.iloc[i]
        gp_tran = get_feature_ancestor(gff_gp, gp_exon, "transcript")
        for j, ta_exon in overlapping_ta_exons(gff_ta, gp_exon).iterrows():
        #for j, ta_exon in overlapping_ta_exons(gff_ta, gp_exon).iterrows():
            if check(ta_exon, gp_exon, gp_tran, known_strand, max_exon_length):
                yield i, j
                break
        i = next_feature_index(gff_gp, i, "exon")

def update(gff: DataFrame, gp_exon: Series, ta_exon: Series, type: str):
    """
    Updating feature in GFF with UTR-extension

    Args:
        gff (DataFrame):    Part of the gene prediction GFF with features of the same
                            sequence
        gp_exon (Series):   An exon from the gene prediction
        ta_exon (Series):   An exon from the transcripome assembly with UTR-extension for
                            exon_gp
        type (str):         "exon", to update the gp_exon itself
                            <type>, to update an ancestor of gp_exon of the specified type                                  
    """
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
                log_file: str) -> tuple[DataFrame,int]:
    """
    Finds matching exons from the transcriptome assembly that pose UTR-extensions for
    exons from the gene prediction and updates the "source", "start" and "end" fields of
    the exons, their corresponding transcripts accordingly.

    Args:
        gff_gp (DataFrame):     Part of the gene prediction GFF with features of the same
                                sequence
        gff_ta (DataFrame):     Part of the transcriptome assembly GFF with features of
                                the same sequence
        known_strand (bool):    Only use exons if the strands of both exons are known.
        max_ex_len (int):       Length limit for exons from the transcriptome assembly
        log_file (str):         Path to log-file

    Returns:
        tuple[DataFrame,int]:   Updated GFF with extended features and number of extensions
    """

    n = 0

    logging_setup("info", log_file)
    logger = getLogger()

    if gff_gp.shape[0]==0 or gff_ta.shape[0]==0:
        return gff_gp, 0
    
    for i, j in exon_matches_sensitive(gff_gp, gff_ta, known_strand, max_ex_len):

        gp_exon = gff_gp.iloc[i]
        ta_exon = gff_ta.iloc[j]

        for type in ["exon", "transcript", "gene"]:
            update(gff_gp, gp_exon, ta_exon, type)
        
        n += 1

    logger.info("\tDone {0:<20} (Extended exons: {1:>5})".format(gff_gp.iloc[0,0],n))
    return gff_gp, n

def exon_extend_threaded(args: tuple[DataFrame,DataFrame,int,bool,int]) -> tuple[DataFrame,int]:
    """
    Wrapper function for exon_extend() allowing it to be used with multiprocessing
    """
    gff_gp, gff_ta, known_strand, max_ex_len, log_file = args

    return exon_extend(gff_gp, gff_ta, known_strand, max_ex_len, log_file)

def exon_extend_multithreaded(gff_prediction: DataFrame,
                              gff_assembly: DataFrame,
                              seqnames: list[str],
                              args) -> list[DataFrame]:
    
    info("Extending exons ...")
    mp_args = zip([gff_prediction[s]        for s in seqnames],
                  [gff_assembly[s]          for s in seqnames],
                  [args.know_strand         for _ in seqnames],
                  [args.max_ex_len          for _ in seqnames],
                  [f"{args.gff_utrpy}.log"  for _ in seqnames])
    
    with Pool(processes=args.threads) as pool:
        results = pool.map(exon_extend_threaded, mp_args)
        gffs    = [r[0] for r in results]

    info(f"Done all  (Extended exons: {sum([r[1] for r in results])})")
    return gffs