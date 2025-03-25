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
    Recursive upwards traversal of the GFF-file to find the ancestor of the specified type
    for the input feature
    """
    if (feature[2])==ancestor_type:
        return feature
    
    parent_id = search(r'Parent=([^;]+)', feature[8])
    if parent_id:
        parent_id = parent_id.group(0).split("=")[1]
        parent    = gff.loc[gff[8].str.contains(f"ID={parent_id}", na=False)].iloc[0]
        ancestor  = get_feature_ancestor(gff, parent, ancestor_type)
        if not ancestor is None:
            return ancestor
    error(f"Failed finding ancestor of type {ancestor_type} for feature {list(feature)}")

def merge_exons(gff_gp: DataFrame,
                ta_exon: Series,
                i: int):
    
    gff_gp.iloc[i,1] = f"{gff_gp.iloc[i,1]} + {ta_exon[1]} (UTRpy)"
    gff_gp.iloc[i,3] = ta_exon[3]
    gff_gp.iloc[i,4] = ta_exon[4]

def update_length(feature: Series, new_exon: Series, gff: DataFrame):

    try:
        row_idx = feature.index[0]
        gff.at[row_idx, 3] = min(gff.at[row_idx, 3], new_exon[3])
        gff.at[row_idx, 4] = max(gff.at[row_idx, 4], new_exon[4])
    except:
        feature

def utr_extend(gff_gp: DataFrame, gff_ta:DataFrame, i: int, j: int):

    gp_exon = gff_gp.iloc[i]
    ta_exon = gff_ta.iloc[j]

    tran = get_feature_ancestor(gff_gp, gp_exon, "transcript")
    gene = get_feature_ancestor(gff_gp, gp_exon, "gene")

    merge_exons(gff_gp, ta_exon, i)
    update_length(gene, ta_exon, gff_gp)
    update_length(tran, ta_exon, gff_gp)