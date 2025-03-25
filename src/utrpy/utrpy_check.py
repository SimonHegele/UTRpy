"""
Functions to check if an exon from the assembly provides an UTR-extension
"""

from pandas import DataFrame, Series

def check_feature_length(feature: Series, max_feature_length: int) -> bool:
    return feature[4]- feature[3] <= max_feature_length

def check_strandedness(feature1: Series, feature2: Series, strict=False):
    """
    Args:
        feature1 (Series): A pandas Series representing a feature from a GFF file (or part of it)
        feature2 (Series): A pandas Series representing a feature from a GFF file (or part of it)

    Returns:
        bool: If strict: Only returns True if strands are known to be of the same strand
              Else:      Also returns True if strands aren not known to be different
    """
    if strict:
        if not feature1[6]=="." or feature2[6]==".":
            if feature1[6]==feature2[6]:
                return True
    else:
        if feature1[6]=="." or feature2[6]==".":
            return True
        if feature1[6]==feature2[6]:
            return True
    return False

def check_extends(ta_exon: Series, gp_exon: Series, tran: Series) -> bool:
    """
    Args:
        ta_exon (Series):    A pandas Series, an exon from the transcriptome assembly
        gp_exon (Series):    A pandas Series, an exon from the gene prediction
        tran (Series):       A pandas Series, a transcript from the gene prediction
        conservative (bool): Only use exons if the strands of both exons are known.

    Returns:
        bool: True if ta_exon extends tran to one side without extending gp_exon on the
              other side
              False else
    """
    # Checks if the ta_exon extends the transcript to the left
    if (int(ta_exon[3])<int(tran[3])):
        # Checks if the transcript end at the same position
        if (int(ta_exon[4])==int(gp_exon[4])):
            return True
    # Checks if the ta_exon extends the transcript to the left
    if (int(ta_exon[4])>int(tran[4])):
        # Checks if the transcript starts at the same position     
        if (int(ta_exon[3])==int(gp_exon[3])):    
            return True
    return False

def check(ta_exon: Series,
          gp_exon: Series,
          tran: Series,
          strict: bool,
          max_exon_length: int) -> bool:
    """
    Checks if the exon from the transcript assembly is a legitimate extension for the
    given transcript

    Args:
        ta_exon (Series):   A pandas Series, an exon from the transcriptome assembly
        gp_exon (Series):   A pandas Series, an exon from the gene prediction
        tran (Series):      A pandas Series, a transcript from the gene prediction
        strict (bool):      Only use exons if the strands of both exons are known.
        max_exon_length:    Length limit for ta_exon  
    Returns:
        bool: True if ta_exon extends the transcript
    """
    if not check_feature_length(ta_exon, max_exon_length):
        return False
    if not check_strandedness(ta_exon, gp_exon, strict):
        return False
    if not check_extends(ta_exon, gp_exon, tran):
        return False
    return True