import logging
import pandas
import numpy

from .utrpy_gff_utils           import attributes_dict
from .utrpy_transcript          import Transcript
from .utrpy_transcript_matching import transcript_matches
from .utrpy_utr_variant         import utr_variant

def select_from_variants(variants: list[dict], select: str):

    if select == "all":
        return variants
    
    exons   = [v["transcript"].loc[v["type"]=="exon"] for v in variants]
    lengths = [numpy.sum(e["end"]-e["start"]-1) for e in exons]
        
    match select:
        case "shortest":
            return [variants[lengths.index(numpy.min(lengths))]]
        case "longest":
            return [variants[lengths.index(numpy.min(lengths))]]
        
def update_gene_lengths(utr_variants: list, p_gff: pandas.DataFrame):

    for variant in utr_variants:

        i = variant["gene"].name
        p_gff.iloc[i,3] = min(p_gff.iloc[i,3], variant["transcript"]["start"].min())
        p_gff.iloc[i,4] = min(p_gff.iloc[i,4], variant["transcript"]["end"].max())

def delete_original_transcripts(p_gff: pandas.DataFrame,
                                to_delete: list[Transcript]):
    
    for transcript in to_delete:
            p_gff = transcript.delete(p_gff)

    return p_gff
    
def utr_extend(p_gff: pandas.DataFrame,
               a_gff: pandas.DataFrame,
               know_strand: bool,
               match_middle_exons: bool,
               keep: bool,
               select: str,
               max_exon_length: int) -> pandas.DataFrame:
    
    p_transcripts    = p_gff.loc[p_gff["type"]=="transcript"]
    all_utr_variants = []
    to_delete        = []

    for i, p_transcript in p_transcripts.iterrows():

        transcript_id = attributes_dict(p_transcript)["ID"]

        matches = list(transcript_matches(p_transcript,
                                          p_gff,
                                          a_gff,
                                          know_strand,
                                          match_middle_exons,
                                          max_exon_length))
        
        if any(matches):

            utr_variants = [utr_variant(match, p_gff, i)
                            for i, match in enumerate(matches)]
            utr_variants = [v for v in utr_variants if not v is None]
            
            if any(utr_variants):

                utr_variants = select_from_variants(utr_variants, select)
                to_delete.append(matches[0]["p_transcript"])
                all_utr_variants += utr_variants
        else:
            utr_variants = []
            
        logging.info(f"{transcript_id:<70} {len(utr_variants)} UTR-variants")

    update_gene_lengths(all_utr_variants, p_gff)

    if not keep:
        p_gff = delete_original_transcripts(p_gff, to_delete)
        
    return pandas.concat([p_gff] + [v["transcript"] for v in all_utr_variants])

def utr_extend_threaded(args):

    p_gff, a_gff, know_strand, match_middle_exons, keep, select, max_exon_length = args

    return utr_extend(p_gff, a_gff, know_strand, match_middle_exons, keep, select, max_exon_length)