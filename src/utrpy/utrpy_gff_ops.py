"""
Functions to operate on pandas DataFrames representing GFF-files
"""

from logging import info
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