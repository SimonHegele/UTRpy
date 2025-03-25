"""
Functions to merge scaffold separated parts of the input GFF-files
"""

from logging import info
from pandas  import DataFrame
from typing  import Generator

from .utrpy_gff_ops import next_feature_index, get_feature_ancestor, utr_extend
from .utrpy_check   import check

def exon_matches(gff_gp: DataFrame,
                 gff_ta: DataFrame,
                 strict: bool,
                 max_exon_length: int) -> Generator:
    """
    Args:
        gff_gp (DataFrame):     A pandas DataFrame representing a GFF file (or part of it)
        gff_ta (DataFrame):     A pandas DataFrame representing a GFF file (or part of it)
        strict (bool):          Only use exons if the strands of both exons are known.
        max_exon_length (int):  Length limit for exons from the transcriptome assembly

    Yields:
        Generator: 2-tuples (i,j) with i referring to an exon gp_exon from the gene
                   gene prediction and j referring to an exon ta_exon from the transcriptome
                   assembly assembly where ta_exon is an UTR extension of gp_exon.
    """
    i = next_feature_index(gff_gp, -1, "exon")
    j = next_feature_index(gff_ta, -1, "exon")
    while i != None and j != None:
        gp_exon = gff_gp.iloc[i]
        ta_exon = gff_ta.iloc[j]
        gp_tran = get_feature_ancestor(gff_gp, gp_exon, "transcript")
        if check(ta_exon, gp_exon, gp_tran, strict, max_exon_length):
            yield i, j
            i = next_feature_index(gff_gp, i, "exon")
            j = next_feature_index(gff_ta, j, "exon")
        else:
            if gp_exon[3] < ta_exon[3]:
                i = next_feature_index(gff_gp, i, "exon")
            else:
                j = next_feature_index(gff_ta, j, "exon")

def exon_extend(gff_gp: DataFrame,
                gff_ta: DataFrame,
                strict: bool,
                max_exon_length: int) -> tuple[DataFrame,int]:
    """
    Extends exons from the gene prediction using suitable exons from the transcriptome
    assembly 

    Args:
        gff_gp (DataFrame): A pandas DataFrame representing the GFF file (or part of it)
                            from the gene prediction
        gff_ta (DataFrame): A pandas DataFrame representing the GFF file (or part of it)
                            from the transcriptome assembly
        strict (bool)
        min_overlap (int):  Minimum overlap of pairs of exons to be considered

    Returns:
        tuple[DataFrame, int]: Updated GFF with extended features and number of extensions
    """
    n = 0

    if gff_gp.shape[0]==0 or gff_ta.shape[0]==0:
        return gff_gp, 0
                  
    # Sorting features by start positions
    gff_gp = gff_gp.sort_values(by=gff_gp.columns[3])
    gff_ta = gff_ta.sort_values(by=gff_ta.columns[3])
    
    for i, j in exon_matches(gff_gp, gff_ta, strict, max_exon_length):

        utr_extend(gff_gp, gff_ta, i, j)
        n += 1
           
    info("\tDone {0:<20} (Extended exons: {1:>5})".format(gff_gp.iloc[0,0],n))
    return gff_gp, n

def exon_extend_threaded(args: tuple[DataFrame,DataFrame,int,bool,int]) -> tuple[DataFrame,int]:
    """
    Wrapper function for process allowing it to be used with multiprocessing
    """
    gff_gp, gff, strict, max_exon_length = args

    return exon_extend(gff_gp, gff, strict, max_exon_length)