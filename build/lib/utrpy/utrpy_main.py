"""
Script Name:    utrpy_main.py
Description:    Main script of UTRpy
                Note that this script is not intended for direct execution,
                please refer to the README.md for instructions on installation and usage.
Author:         Simon Hegele
Date:           2025-04-01
Version:        0.2
License:        GPL-3
"""

from logging                import info
from time                   import time

from .utrpy_argumentparser  import UTRpyArgparser
from .utrpy_exon_extend     import exon_extend_multithreaded
from .utrpy_io              import load_data, write_data
from .utrpy_logging         import logging_setup
from .utrpy_utr_add         import add_utrs_multithreaded

def time_format(seconds):

    seconds = int(seconds)
    hours,   seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"{hours:02}h:{minutes:02}m:{seconds:02}s"

def main():

    start       = time()
    utrpy_args  = UTRpyArgparser().parse_args()

    logging_setup("info", f"{utrpy_args.gff_utrpy}.log")

    gff_prediction, gff_assembly, seqnames = load_data(utrpy_args)

    gff_utrpy = exon_extend_multithreaded(gff_prediction,
                                          gff_assembly,
                                          seqnames,
                                          utrpy_args)
    
    if utrpy_args.explicit:
        gff_utrpy = add_utrs_multithreaded(gff_utrpy, utrpy_args.threads)

    write_data(utrpy_args, gff_utrpy, seqnames)

    info(f"Completed after {time_format(time()-start)}")
    info("#############################################")
    info("#    Simon says: Thanks for using UTRpy!    #")
    info("#############################################")
