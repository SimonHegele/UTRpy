""" 
Functions to explicitly add UTRs as features
"""

from logging import info
from pandas  import DataFrame
from re      import search

def add_utr(gff, type, start, stop, strand, chr, i, transcript_id):
    if not start == stop:
        gff.loc[len(gff)] = [chr, "UTRpy", type, start, stop, '.', strand, '.',
                            f"ID=3pUTR_{chr}_{i};Parent={transcript_id}"]

def add_3_prime_utrs(gff: DataFrame):

    transcripts  = gff.loc[gff[2]=="transcript"]
    stop_codons  = gff.loc[gff[2]=="stop_codon"]
    type         = "3'-UTR"

    for i, transcript in transcripts.iterrows():

        transcript_id = search(r'ID=([^;]+)', transcript[8]).group(0).split("=")[1]
        stop_codon    = stop_codons.loc[stop_codons[8].str.contains(transcript_id)]
        chr           = transcript[0]
        strand        = transcript[6]

        if stop_codon.shape[0] > 0 and strand != '.':
            stop_codon = stop_codon.iloc[0]
            if (strand == '+') and (stop_codon[4] < transcript[4]):
                start  = stop_codon[4] + 1
                stop   = transcript[4]
                add_utr(gff, type, start, stop, strand, chr, i, transcript_id)
            if (strand == '-') and (stop_codon[3] > transcript[3]):
                start  = transcript[3]
                stop   = stop_codon[3] -1
                add_utr(gff, type, start, stop, strand, chr, i, transcript_id)

def add_5_prime_utrs(gff: DataFrame) -> DataFrame:

    transcripts  = gff.loc[gff[2]=="transcript"]
    start_codons = gff.loc[gff[2]=="start_codon"]
    type         = "5'-UTR"

    for i, transcript in transcripts.iterrows():

        transcript_id = search(r'ID=([^;]+)', transcript[8]).group(0).split("=")[1]
        start_codon   = start_codons.loc[start_codons[8].str.contains(transcript_id)]
        chr           = transcript[0]
        strand        = transcript[6]

        if start_codon.shape[0] > 0 and strand != '.':
            start_codon = start_codon.iloc[0]
            if (strand == '+') and (start_codon[3] > transcript[3]):
                start  = transcript[3] 
                stop   = start_codon[3] - 1
                add_utr(gff, type, start, stop, strand, chr, i, transcript_id)
            if (strand == '-') and (start_codon[4] < transcript[4]):
                start  = start_codon[4] -1 
                stop   = transcript[4]
                add_utr(gff, type, start, stop, strand, chr, i, transcript_id)
    
def add_utrs(gff: DataFrame) -> tuple[DataFrame, int]:

    n = gff.shape[0]

    add_3_prime_utrs(gff)
    add_5_prime_utrs(gff)
    gff = gff.sort_values(by=gff.columns[3])
    gff = gff.reset_index(drop=True)

    n = gff.shape[0]-n
    info("\tDone {0:<20} (Annotated UTRs: {1:>5})".format(gff.iloc[0,0],n))
    return gff, n


    
