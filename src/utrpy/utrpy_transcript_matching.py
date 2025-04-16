"""
Module Name:    utrpy_argumentparser.py
Description:    Provides methods to match exons of a predicted transcript to the exons of
                an assembled transcript.
Author:         Simon Hegele
Date:           2025-04-01
Version:        1.0
License:        GPL-3
"""

import pandas
import typing
import logging

from .utrpy_gff_utils  import check_strands, included_features
from .utrpy_transcript import Transcript

def first_matching_exon(p_transcript: Transcript,
                        a_transcript: Transcript) -> int | None:
    
    for i in range(len(a_transcript.exons)):
        if a_transcript.exons.iloc[i]["start"] >= p_transcript.exons.iloc[0]["start"]:
            break
        if a_transcript.exons.iloc[i]["end"] == p_transcript.exons.iloc[0]["end"]:
            return i
        
def last_matching_exon(p_transcript: Transcript,
                       a_transcript: Transcript) -> int | None:
    
    for i in reversed(range(len(a_transcript.exons))):
        if a_transcript.exons.iloc[i]["end"] <= p_transcript.exons.iloc[-1]["end"]:
            break
        if a_transcript.exons.iloc[i]["start"] == p_transcript.exons.iloc[-1]["start"]:
            return i
        
def all_middle_exons_match(p_transcript: Transcript,
                           a_transcript: Transcript,
                           i: int,
                           k: int) -> bool:
    
    for j in range(i+1,k-1):
        p_exon = p_transcript.exons.iloc[j-i]
        a_exon = a_transcript.exons.iloc[j]
        if a_exon["start"] != p_exon["start"]:
            return False
        if a_exon["end"] != p_exon["end"]:
            return False
    return True

def match_single_exon(p_transcript: Transcript,
                      a_transcript: Transcript) -> dict | None:
    
    p_exon = p_transcript.exons.iloc[0]
    for i in range(len(a_transcript.exons)):
        a_exon = a_transcript.exons.iloc[i]
        if a_exon["start"] < p_exon["start"]:
            if a_exon["end"] > p_exon["end"]:
                return {"p_transcript": p_transcript,
                        "a_transcript": a_transcript,
                        "start"       : i,
                        "end"         : i}
    return

def match_multiple_exons(p_transcript: Transcript,
                         a_transcript: Transcript,
                         match_middle_exons) -> dict | None:
    
    i = first_matching_exon(p_transcript, a_transcript)
    k = last_matching_exon(p_transcript, a_transcript)

    if (i is None) or (k is None):
        return
    if match_middle_exons and not all_middle_exons_match(p_transcript, a_transcript, i, k):
        return
    return {"p_transcript": p_transcript,
            "a_transcript": a_transcript,
            "start"       : i,
            "end"         : k}

def transcript_match(p_transcript: Transcript,
                     a_transcript: Transcript,
                     match_middle_exons) -> dict | None:
    
    if len(p_transcript.exons) == 1:
        return match_single_exon(p_transcript, a_transcript)
    else:
        return match_multiple_exons(p_transcript, a_transcript, match_middle_exons)
    
def check_exon_lengths(transcript: Transcript, max_exon_length) -> bool:

    exon_lenghts = transcript.exons["end"] - transcript.exons["start"] + 1

    return max_exon_length > exon_lenghts.max()
 
def transcript_matches(p_transcript: Transcript,
                       p_gff: pandas.DataFrame,
                       a_gff: pandas.DataFrame,
                       know_strand: bool,
                       match_middle_exons: bool,
                       max_exon_length: int) -> typing.Generator:
    
    a_transcripts = included_features(a_gff, p_transcript, "RNA")

    a_transcripts = a_transcripts[
        a_transcripts.apply(lambda t: check_strands(t,
                                                    p_transcript,
                                                    know_strand),
                            axis=1)]
    
    p_transcript = Transcript(p_transcript, p_gff)

    if len(p_transcript.exons) == 0:
        logging.error(f"No exons found for transcript {p_transcript.id}")
        return

    for i, a_transcript in a_transcripts.iterrows():

        a_transcript = Transcript(a_transcript, a_gff)

        if check_exon_lengths(a_transcript, max_exon_length):
            m = transcript_match(p_transcript, a_transcript, match_middle_exons)
            if not m is None:
                yield m