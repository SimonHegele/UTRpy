"""
Module Name:    utrpy_gff_utils.py
Description:    Provides functions to retrive information from GFF-files represented as
                pandas.DataFrame
Author:         Simon Hegele
Date:           2025-04-01
Version:        0.2
License:        GPL-3
"""

from pandas import DataFrame, read_csv, Series

gff_columns = ["seqname",
               "source",
               "feature",
               "start",
               "end",
               "score",
               "strand",
               "frame",
               "attributes"]

def load_gff(file_path: str) -> DataFrame:
    """
    Args:
        file_path (str): location of the GFF-file in the file system

    Returns:
        DataFrame: GFF-file
    """
    return read_csv(file_path, sep="\t", header=None, comment="#", names=gff_columns)

def next_feature_index(gff: DataFrame, i: int, feature_type:str) -> int | None:
    """
    Args:
        gff (DataFrame):    A pandas DataFrame representing a GFF file (or part of it)
        i (int):            An index refering to a position in the GFF file (or part of it)
        feature_type (str): A feature type

    Returns:
        int | None: Next position of a row with the specified feature type or None if there
                    is no such index.
    """
    while True:
        i += 1
        if i == gff.shape[0]:
            return
        if gff.iloc[i,2]==feature_type:
            return i
        
def get_attributes_dict(feature: Series) -> dict:
    """
    Parses the attributes field of a feature to dictionary
    """
    return {a.split("=")[0]: a.split("=")[1] for a in feature["attributes"].split(";")}

def get_feature_string(feature: dict) -> str:
    """
    Parses feature to tab separated string
    """
    return "\t".join(str(feature[c]) for c in gff_columns)

def get_seqnames(gff: DataFrame) -> list[str]:
    
    return gff["seqname"].unique().tolist()

def seqname_split_sort(gff: DataFrame, seqnames=None) -> dict[str, DataFrame]:
    """
    Args:
        gff (DataFrame):                GFF-file
        seqnames (list[str], optional): List of seqnames. Defaults to None.

    Returns:
        dict[str, DataFrame]: key:   seqname
                              value: DataFrame with features from seqname sorted by
                                     start 
    """
    if seqnames is None:
        seqnames = get_seqnames(gff)

    return {seqname: gff.loc[gff["seqname"]==seqname]
            .sort_values("start")
            .reset_index(drop=True)
            for seqname in seqnames}

def get_feature_ancestor(gff: DataFrame, feature: dict, ancestor_type: str) -> Series:
    """
    Returns the ancestor of the specified type for the input feature
    a) Directly, if <ancestor_type>_id=ancestor is in the attributes column or
    b) By recursively traversing the GFF-hierarchy following the child-parent relationship

    Args:
        gff (DataFrame):        
        feature (dict):        
        ancestor_type (str):    

    Raises:
        Exception: When an ancestor of the specified type could not be found.
                   -> Wrong use (e.g. searching ancestor of type exon for feature of type
                      gene)
                   -> ill-formatted GFF

    Returns:
        Series: ancestor feature of the specified type for the input feature
    """
    
    if ancestor_type==(feature["feature"]):
        return feature

    attributes = get_attributes_dict(feature)

    if f"{ancestor_type}_id" in attributes.keys():
        ancestor_id = attributes[f"{ancestor_type}_id"]
        ancestor    = gff.loc[gff["attributes"].str.contains(f"ID={ancestor_id}", na=False)].iloc[0]
        return ancestor
    
    elif "Parent" in attributes.keys():
        parent_id = attributes["Parent"]
        parent    = gff.loc[gff["attributes"].str.contains(f"ID={parent_id}", na=False)].iloc[0]
        return get_feature_ancestor(gff, parent, ancestor_type)
    
    fs = get_feature_string(feature)
    raise Exception(f"Unable to find ancestor of type {ancestor_type} for feature\n{fs}")
