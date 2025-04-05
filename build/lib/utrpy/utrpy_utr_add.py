"""
Script Name:    utrpy_utr_add.py
Description:    Determine UTRs from differences of transcript starts / ends relative to
                to start codons and stop codons and explicitly adding them as features.
Author:         Simon Hegele
Date:           2025-04-01
Version:        1.0
License:        GPL-3
"""

from logging            import info
from multiprocessing    import Pool
from pandas             import concat, DataFrame

from .utrpy_gff_utils   import get_attributes_dict

##########################################################################################
#                                                                                        #
# Orientation help:                                                                      #
#                                                                                        #
#   5'---UTR---start----------stop---UTR---3' + strand                                   #
#                                                                                        #
#   3'---UTR---stop----------start---UTR---5' - strand                                   #
#                                                                                        #
##########################################################################################

def add_utr(gff: DataFrame,
            transcript: dict,
            type: str,
            utr_start: int,
            utr_end: int,
            i: int) -> DataFrame:
    
    transcript_id = get_attributes_dict(transcript)["ID"]
    seqname       =  transcript["seqname"]
    
    return concat([gff, DataFrame({
        "seqname":      [seqname],
        "source":       [transcript["source"]],
        "feature":      [type],
        "start":        [min(utr_start, utr_end)],
        "end":          [max(utr_start, utr_end)],
        "score":        ["."],
        "strand":       [transcript["strand"]],
        "frame":        ["."],
        "attributes":   [f"ID={seqname}_{type}_{i};Parent={transcript_id}"]
    })])

def add_utrs(gff: DataFrame) -> tuple[DataFrame, int]:
    """
    Args:
        gff (DataFrame):        GFF-file (or a part of it) represented as pandas DataFrame

    Returns:
        tuple[DataFrame, int]:  The input DataFrame with added UTRs

    Yes this method is actually three or more methods under a single trenchcoat.
    """

    n = gff.shape[0]

    transcripts  = gff.loc[gff["feature"]=="transcript"]
    start_codons = gff.loc[gff["feature"]=="start_codon"]
    stop_codons  = gff.loc[gff["feature"]=="stop_codon"]

    for i, transcript in transcripts.iterrows():

        start_codon = start_codons.loc[start_codons["attributes"]
                                       .str.contains(get_attributes_dict(transcript)["ID"])]
        stop_codon  = stop_codons.loc[stop_codons["attributes"]
                                      .str.contains(get_attributes_dict(transcript)["ID"])]
        
        if transcript["strand"] == "+":
            if len(start_codon) and (transcript["start"] < start_codon.iloc[0]["start"]):
                gff = add_utr(gff,
                              transcript,
                              "five_prime_UTR",
                              transcript["start"],
                              start_codon.iloc[0]["start"]-1,
                              i)
            if len(stop_codon) and (stop_codon.iloc[0]["end"] < transcript["end"]):
                gff = add_utr(gff,
                              transcript,
                              "three_prime_UTR",
                              stop_codon.iloc[0]["end"]+1,
                              transcript["end"],
                              i)
        if transcript["strand"] == "-":
            if len(stop_codon) and (transcript["start"] < stop_codon.iloc[0]["start"]):
                gff = add_utr(gff,
                              transcript,
                              "three_prime_UTR",
                              transcript["start"],
                              stop_codon.iloc[0]["start"]-1,
                              i)
            if len(start_codon) and (start_codon.iloc[0]["end"] < transcript["end"]):
                gff = add_utr(gff,
                              transcript,
                              "five_prime_UTR",
                              start_codon.iloc[0]["end"]+1,
                              transcript["end"],
                              i)
                
    n = gff.shape[0]-n
    info("\tDone {0:<20} (Annotated UTRs: {1:>5})".format(gff.iloc[0,0],n))

    return gff, n

def add_utrs_multithreaded(gff: list[DataFrame], threads):

    info("Adding UTRs...")

    with Pool(processes=threads) as pool:
        results = pool.map(add_utrs, gff)
        gff    = [r[0] for r in results]

    info(f"Done all  (Annotated UTRs: {sum([r[1] for r in results])})")

    return gff