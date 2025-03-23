from pandas import DataFrame, Series
from typing import Generator

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
        
def features_overlap(feature1: Series, feature2: Series, min_overlap: int) -> bool:
    """
    Determines if two features overlap

    Args:
        feature1 (Series): A pandas Series representing a feature from a GFF file (or part of it)
        feature2 (Series): A pandas Series representing a feature from a GFF file (or part of it)
        min_overlap (int): The minimum required overlap length

    Returns:
        bool: True if the features overlap, False if not
    """
    return max(feature1[3], feature2[3]) + min_overlap <= min(feature1[4], feature2[4])

def feature_matches(gff1: DataFrame,
                    gff2: DataFrame,
                    feature_type1: str,
                    feature_type2: str,
                    min_overlap) -> Generator:
    """
    Args:
        gff1 (DataFrame):            A pandas DataFrame representing a GFF file (or part of it)
        gff2 (DataFrame):            A pandas DataFrame representing a GFF file (or part of it)
        feature_type1 (str):         A feature type
        feature_type2 (str):         A feature type
        min_overlap (int, optional): The minimum required overlap length

    Yields:
        Generator: 2-tuples (i,j) refering to indices of overlapping features in gff1 and gff2.
    """
    i = next_feature_index(gff1, -1, feature_type1)
    j = next_feature_index(gff2, -1, feature_type2)
    while i != None and j != None:
        feature1 = gff1.iloc[i]
        feature2 = gff2.iloc[j]
        if features_overlap(feature1, feature2, min_overlap):
            yield i, j
            i = next_feature_index(gff1, i, feature_type1)
            j = next_feature_index(gff2, j, feature_type2)
        else:
            if feature1[3] < feature2[3]:
                i = next_feature_index(gff1, i, feature_type1)
            else:
                j = next_feature_index(gff2, j, feature_type2)

