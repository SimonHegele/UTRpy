"""
Module Name:    utrpy_variant_insertion.py
                Provides methods for the creation of a DataFrame representing an UTR -
                variant and its features merged from the assembly and prediction 
Description:    
Author:         Simon Hegele
Date:           2025-04-01
Version:        1.0
License:        GPL-3
"""

import pandas

from .utrpy_gff_utils  import attributes_dict, attributes_str, get_ancestor, features_overlap

def purely_assembled_exons(transcript_match):

    a_transcript = transcript_match["a_transcript"]

    return pandas.concat([a_transcript.exons.iloc[:transcript_match["start"]+1],
                          a_transcript.exons.iloc[transcript_match["end"]:]])

def assembled_and_predicted_exons(transcript_match):

    return transcript_match["p_transcript"].exons.iloc[1:-1]

def purely_predicted_features(transcript_match):

    p_transcript = transcript_match["p_transcript"]

    return p_transcript.features.loc[p_transcript.features["type"] != "exon"]

def create_transcript_id(transcript_match: dict,
                         variant) -> str:

    return attributes_dict(transcript_match["p_transcript"].data)["ID"] + f"_utr_{variant}"

def combined_features(transcript_match: dict) -> pandas.DataFrame:

    return pandas.concat([purely_assembled_exons(transcript_match),
                          assembled_and_predicted_exons(transcript_match),
                          purely_predicted_features(transcript_match)]).reset_index(drop=True)

def annotate_features(features: pandas.DataFrame,
                      p_transcript: pandas.Series,
                      transcript_id: str,
                      gene_id: str,
                      assembler: str,
                      predictor: str):

    for i, feature in features.iterrows():

        attributes = attributes_dict(feature)

        attributes["ID"]            = f"{transcript_id}_feature_{i}"
        attributes["Parent"]        = transcript_id
        attributes["gene_id"]       = gene_id
        attributes["transcript_id"] = transcript_id
        attributes["assembler"]     = assembler
        
        if features_overlap(p_transcript.data, feature):
            attributes["predictor"] = predictor

        features.iloc[i,1] = "UTRpy"
        features.iloc[i,8] = attributes_str(attributes)

def build_transcript_row(p_transcript,
                         tran_id,
                         gene_id,
                         assembler,
                         predictor,
                         gene,
                         features) -> pandas.DataFrame:

    attributes = attributes_dict(p_transcript.data)
    attributes["ID"]            = tran_id
    attributes["Parent"]        = gene_id
    attributes["gene_id"]       = gene_id
    attributes["assembler"]     = assembler
    attributes["predictor"]     = predictor
         
    return pandas.DataFrame({"seqname":   [gene["seqname"]],
                             "source":    ["UTRpy"],
                             "type":      ["transcript"],
                             "start":     [features["start"].min()],
                             "end":       [features["end"].max()],
                             "score":     ["."],
                             "strand":    [gene["strand"]],
                             "frame":     ["."],
                             "attributes": attributes_str(attributes)})

def utr_variant(transcript_match: dict,
                p_gff: pandas.DataFrame,
                variant: int) -> dict[str, pandas.DataFrame | pandas.Series] | None:
    
    p_transcript  = transcript_match["p_transcript"]
    a_transcript  = transcript_match["a_transcript"]
    gene          = get_ancestor(p_gff, transcript_match["p_transcript"].data, "gene")

    if gene is None:
        return None
    
    gene_id       = "???" if gene is None else attributes_dict(gene)["ID"]
    transcript_id = create_transcript_id(transcript_match, variant)
    predictor     = p_transcript.exons.iloc[0]["source"]
    assembler     = a_transcript.exons.iloc[0]["source"]
    features      = combined_features(transcript_match)

    annotate_features(features,p_transcript,transcript_id,gene_id,assembler,predictor)

    return {"gene":       gene,
            "transcript": pandas.concat([features,
                                         build_transcript_row(p_transcript,
                                                              transcript_id,
                                                              gene_id,
                                                              assembler,
                                                              predictor,
                                                              gene,
                                                              features)])
           } 
    