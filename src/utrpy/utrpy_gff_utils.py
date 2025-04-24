"""
Module Name:    utrpy_gff_utils.py
Description:    Functions for pandas.Dataframe represented GFF-files
                - attributes_dict(feature: pandas.Series) -> dict[str, str]
                - attributes_str(attributes: dict[str, str]) -> str
                - check_strands(feature1: pandas.Series, feature2: pandas.Series, know_strand=False) -> bool
                - empty_gff() -> pandas.DataFrame
                - get_ancestor(gff: pandas.DataFrame, feature: pandas.Series, type: str) -> pandas.Series:
                - features_overlap(feature_1, feature_2) -> bool
                - load_gff(file_path: str) -> pandas.DataFrame
                - seqname_split(gff: pandas.DataFrame, seqnames=None) -> dict[str, pandas.DataFrame]
                - write_gff(gff: pandas.DataFrame, file_path: str) -> None
Author:         Simon Hegele
Date:           2025-04-01
Version:        1.0
License:        GPL-3
"""

import logging
import pandas

gff_columns = ["seqname",
               "source",
               "type",
               "start",
               "end",
               "score",
               "strand",
               "frame",
               "attributes"]

def attributes_dict(feature: pandas.Series) -> dict[str, str]:
    """
    Parsing the key=value pairs from a features attributes fields into a hashmap
    """
    
    return {a.split("=")[0]: a.split("=")[1] for a in feature["attributes"].split(";")}

def attributes_str(attributes: dict[str, str]) -> str:
    
    return ";".join([f"{key}={attributes[key]}" for key in attributes.keys()])

def check_strands(feature1: pandas.Series,
                  feature2: pandas.Series,
                  know_strand=False) -> bool:

    if know_strand:
        if not feature1["strand"]=="." or feature2["strand"]==".":
            if feature1["strand"]==feature2["strand"]:
                return True
    else:
        if feature1["strand"]=="." or feature2["strand"]==".":
            return True
        if feature1["strand"]==feature2["strand"]:
            return True
    return False

def empty_gff() -> pandas.DataFrame:
    """
    Creates and returns an empty GFF-file
    """

    return pandas.DataFrame({}, columns = gff_columns)

def get_ancestor(gff: pandas.DataFrame,
                 feature: pandas.Series,
                 type: str) -> pandas.Series:
    """
    Returns the ancestor of the specified type for the input feature
    a) Directly, if <ancestor_type>_id=ancestor is in the attributes column or
    b) By recursively following the GFF-hierarchy following the child-parent relationship
    """

    ancestor = feature

    try:
        for _ in range(5): 
            if type==(ancestor["type"]):
                return ancestor
            attrs = attributes_dict(ancestor)
            if f"{type}_id" in attrs.keys():
                id = attrs[f"{type}_id"]
            elif "Parent" in attrs.keys():
                id = attrs["Parent"]
            ancestor = gff[gff.apply(lambda f: attributes_dict(f)["ID"]==id, axis=1)].iloc[0]

        logging.error(f"Failed to find ancestor of type {type} for {list(feature)}")

    except Exception as e:

        logging.error(f"Failed to find ancestor of type {type} for {list(feature)}\n{e}")

def feature_includes(feature_1, feature_2):
    """
    Checks if feature_1 includes feature_2
    """

    if not feature_1["seqname"] == feature_2["seqname"]:
        return False
    if not feature_1["start"] <= feature_2["start"]:
        return False
    if not feature_1["end"] >= feature_2["end"]:
        return False
    return True

def features_overlap(feature_1, feature_2) -> bool:
    """
    Checks if two features overlap
    """

    if feature_1["seqname"] == feature_2["seqname"]:
        if (feature_1["start"]<=feature_2["end"]) and (feature_1["end"]>=feature_2["start"]):
            return True
        if (feature_2["start"]<=feature_1["end"]) and (feature_2["end"]>=feature_1["start"]):
            return True
    return False

def load_gff(file_path: str) -> pandas.DataFrame:
    """
    Loading a GFF-file from the file-system
    """
    return pandas.read_csv(file_path, sep="\t", header=None, comment="#", names=gff_columns)

def included_features(gff: pandas.DataFrame,
                      feature: pandas.Series,
                      type=""):
    """
    Returns features of the specified type from the input GFF included by the input feature 
    """
    mask_seqname  = gff["seqname"] == feature["seqname"]
    mask_type     = gff["type"].str.contains(type, regex=False)
    prefiltered   = gff[mask_seqname & mask_type]
    mask_includes = prefiltered.apply(lambda f: feature_includes(f, feature), axis=1)

    return prefiltered[mask_includes]

def overlapping_features(gff: pandas.DataFrame,
                         feature: pandas.Series,
                         type="") -> pandas.DataFrame:
    """
    Returns features of the specified type from the input GFF overlapping the input feature 
    """
    
    mask_seqname = gff["seqname"] == feature["seqname"]
    mask_type    = gff["type"].str.contains(type, regex=False)
    prefiltered  = gff[mask_seqname & mask_type]
    mask_overlap = prefiltered.apply(lambda f: features_overlap(feature, f), axis=1)

    return prefiltered[mask_overlap]

def seqname_split(gff: pandas.DataFrame, seqnames=None) -> dict[str, pandas.DataFrame]:

    if seqnames is None:
        seqnames = gff["seqname"].unique()

    return {seqname: gff.loc[gff["seqname"]==seqname]
            .sort_values("start")
            .reset_index(drop=True)
            for seqname in seqnames}

def write_gff(gff: pandas.DataFrame, file_path: str) -> None:

    gff.to_csv(file_path, sep="\t", index=False, header=None)
