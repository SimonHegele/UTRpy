from logging    import info, warning
from pandas     import DataFrame, Series
from re         import search

from .utrpy_feature_matches import feature_matches

def extends_transcript(ta_exon: Series,
                       gp_exon: Series,
                       tran: Series,
                       conservative: bool,
                       max_exon_length: int) -> bool:
    """
    Checks if the exon from the transcript assembly is a legitimate extension for the
    given transcript

    Args:
        ta_exon (Series):    A pandas Series, an exon from thetranscriptome assembly
        gp_exon (Series):    A pandas Series, an exon from the gene prediction
        tran (Series):       A pandas Series, a transcript from the gene prediction
        conservative (bool): Only use exons if the strands of both exons are known.

    Returns:
        bool: True if ta_exon extends the transcript
    """
    # Checks if the exon length exceeds what we consider to be reasonable
    if ta_exon[4]- ta_exon[3] <= max_exon_length:
        # Checks strandedness
        if (ta_exon[6]==tran[6]) or (not conservative and (ta_exon[6]=="." or gp_exon[6]==".")):
            # Checks if the ta_exon extends the transcript to the left
            if (int(ta_exon[3])<int(tran[3])):
                # Checks if the transcript end at the same position
                if (int(ta_exon[4])==int(gp_exon[4])):
                    return True
            # Extends transcript to the right
            if (int(ta_exon[4])>int(tran[4])):
                # Starts at same position as old exon     
                if (int(ta_exon[3])==int(gp_exon[3])):    
                    return True
            
    return False

def replace_feature(gff: DataFrame,
                    new_feature: Series,
                    i: int):
    """
    Replaces the old feature at postion i with the new feature while keeping the
    attributes

    Args:
        gff (DataFrame):   A pandas DataFrame representing a GFF file (or part of it)
        old_feature (Series): _description_
        new_feature (Series): _description_
        i (int): _description_
    """
    old_feature   = gff.iloc[i].copy()
    gff.iloc[i]   = Series(new_feature, index=gff.columns)
    gff.iloc[i,8] = old_feature[8]
    gff.iloc[i,1] = f"{old_feature[1]} + {new_feature[1]} (UTRpy)"

    print(f"GP:  {list(old_feature)}")
    print(f"TA:  {list(new_feature)}")
    print(f"New: {list(gff.iloc[i])}")
    print()

def update_length(feature: Series, new_exon: Series, gff: DataFrame):
    try:
        row_idx = feature.index[0]
        gff.at[row_idx, 3] = min(gff.at[row_idx, 3], new_exon[3])
        gff.at[row_idx, 4] = max(gff.at[row_idx, 4], new_exon[4])
    except:
        feature

def update_gff(gff, new_exon, i, trans, gene):

    replace_feature(gff, new_exon, i)
    update_length(gene, new_exon, gff)
    update_length(trans, new_exon, gff)

def exon_extend(gff_gp: DataFrame,
                gff_ta: DataFrame,
                min_overlap: int,
                conservative: bool,
                max_exon_length: int) -> tuple[DataFrame,int]:
    """
    Extends exons from the gene prediction using suitable exons from the transcriptome
    assembly 

    Args:
        gff_gp (DataFrame): A pandas DataFrame representing the GFF file (or part of it)
                            from the gene prediction
        gff_ta (DataFrame): A pandas DataFrame representing the GFF file (or part of it)
                            from the transcriptome assembly
        min_overlap (int):  Minimum overlap of pairs of exons to be considered

    Returns:
        tuple[DataFrame, int]: Updated GFF with extended features and number of extensions
    """
    n = 0
    
    for i, j in feature_matches(gff_gp, gff_ta, "exon", "exon", min_overlap):

        #try:
        gp_exon = gff_gp.iloc[i]
        ta_exon = gff_ta.iloc[j]
        tran_id  = search(r'transcript_id=([^;]+)', gp_exon[8]).group(0).split("=")[1]
        gene_id  = search(r'gene_id=([^;]+)',       gp_exon[8]).group(0).split("=")[1]
        tran     = gff_gp.loc[gff_gp[8].str.contains(f"ID={tran_id}", na=False)]
        gene     = gff_gp.loc[gff_gp[8].str.contains(f"ID={gene_id}", na=False)]

        if extends_transcript(ta_exon,
                                gp_exon,
                                tran.iloc[0],
                                conservative,
                                max_exon_length):
            n += 1
            update_gff(gff_gp, ta_exon, i, tran, gene)
            

        #except:
        #warning(f"failed at {gff_gp.iloc[0,0]}\ngff_gp row {list(gp_exon)}\n gff_ta row {list(ta_exon)}")

    info("\tDone {0:<20} (Extended exons: {1:>5})".format(gff_gp.iloc[0,0],n))
    return gff_gp, n

def exon_extend_threaded(args: tuple[DataFrame,DataFrame,int,bool,int]) -> tuple[DataFrame,int]:
    """
    Wrapper function for exon_extend allowing it to be used with multiprocessing
    """
    gff_gp, gff, min_overlap, conservative, max_exon_length = args

    return exon_extend(gff_gp, gff, min_overlap, conservative, max_exon_length)