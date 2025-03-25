"""
Functions to operate on pandas DataFrames representing GFF-files
"""

from copy    import copy
from logging import error
from pandas  import DataFrame, Series
from re      import search

def next_feature_index(gff: DataFrame, i: int, feature_type:str) -> int | None:
    """
    Determines the next position in the GFF were a row with the specified feature type can
    be found.

    Args:
        gff (DataFrame):    A pandas DataFrame representing a GFF file (or part of it)
        i (int):            An index refering to a position in the GFF file (or part of it)
        feature_type (str): A feature type

    Returns:
        int | None: Next index of a row with the specified feature type or None if there
                    is no such index.
    """
    while True:
        i += 1
        if i == gff.shape[0]:
            return
        if gff.iloc[i,2]==feature_type:
            return i
        
def get_feature_ancestor(gff: DataFrame, feature: Series, ancestor_type: str)->str:
    """
    Returns the ancestor of the specified type for the input feature either directly if 
    <ancestor_type>_id=ancestor is in the attributes column or by recursive upwards
    traversal GFF-hierarchy
    """
    if ancestor_type==(feature[2]):
        return feature

    attrs = feature[8].split(";")
    try:
        ancestor_id = [a for a in attrs if a.startswith(f"{ancestor_type}_id=")][0].split("=")[1]
        ancestor    = gff.loc[gff[8].str.contains(f"ID={ancestor_id}", na=False)].iloc[0]
        return ancestor
    except:
        parent_id   = [a for a in attrs if a.startswith("Parent=")][0].split("=")[1]
        parent      = gff.loc[gff[8].str.contains(f"ID={parent_id}", na=False)].iloc[0]
        ancestor    = get_feature_ancestor(gff, parent, ancestor_type)
        if not ancestor is None:
            return ancestor
    # I'm not paid enough (anything) for proper error handling.
    # If the GFF is fine, then this is too.

def merge_exons(gff_gp: DataFrame,
                ta_exon: Series,
                i: int):
    
    gff_gp.iloc[i,1] = f"{gff_gp.iloc[i,1]} + {ta_exon[1]} (UTRpy)"
    gff_gp.iloc[i,3] = ta_exon[3]
    gff_gp.iloc[i,4] = ta_exon[4]

def update_length(feature: Series, new_exon: Series, gff: DataFrame):

    row_idx = (gff == feature).all(axis=1).idxmax()
    gff.iloc[row_idx, 3] = min(gff.iloc[row_idx, 3], new_exon[3])
    gff.iloc[row_idx, 4] = max(gff.iloc[row_idx, 4], new_exon[4])

def utr_extend(gff_gp: DataFrame, gff_ta:DataFrame, i: int, j: int):

    gp_exon = gff_gp.iloc[i]
    ta_exon = gff_ta.iloc[j]

    tran = get_feature_ancestor(gff_gp, gp_exon, "transcript")
    gene = get_feature_ancestor(gff_gp, tran, "gene")

    merge_exons(gff_gp, ta_exon, i)
    update_length(tran, ta_exon, gff_gp)
    update_length(gene, ta_exon, gff_gp)
    
    tran = get_feature_ancestor(gff_gp, gp_exon, "transcript")
    gene = get_feature_ancestor(gff_gp, tran, "gene")
