"""
Module Name:    utrpy_argumentparser.py
Description:    Provides class UTRpyArgparser(ArgumentParser)
                - Arguments for UTRpy added on initialization
                - Extended parse_args() to 
                    a) Check input 
                    b) Create the output directory and the temporary directory
                    c) Write a file with parameter selection to the output directory
Author:         Simon Hegele
Date:           2025-04-01
Version:        1.0
License:        GPL-3
"""

import argparse
import datetime
import logging
import os

class UTRpyArgparser(argparse.ArgumentParser):

    prog        =  "utrpy"

    description = ("UTR extension of transcript exons from protein orthology based gene "
                   "prediction using exons from reference based assembly")
                       
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        # Input files
        self.add_argument("prediction",
                          help="Annotation from gene prediction (GFF/GTF)")
        self.add_argument("assembly",
                          help="Annotation from transcriptome assembly (GFF/GTF)")
        self.add_argument("outdir",
                          help="Output directory (Must not exist already)")
        self.add_argument("-pp","--pinky_promise",
                          help="Pinky promise that gff_prediction is correct (Will fix it otherwise)",
                          metavar="")
        
        # Transcript matching
        self.add_argument("-m", "--match",
                          help="What exons of predicted transcripts to match"
                          " [choices: ends, all] [default: all]",
                          metavar="",
                          choices=["ends", "all"],
                          default="all")
        self.add_argument("-ks","--know_strand",
                          action="store_true",
                          help="Use only transcripts where the strand is known")
        self.add_argument("-me","--max_exon_length",
                          help="Don't use assembled transcripts with exons longer than this [default: 20000]",
                          metavar="",
                          type=int,
                          default=20000)
        
        # UTR-variant selection
        self.add_argument("-s", "--select",
                          help="How to select UTR-variants if there are multiple [choices: shortest, longest, all] [default: all]",
                          choices=["shortest", "longest", "all"],
                          default="all",
                          metavar="",)
        self.add_argument("-k","--keep",
                          help="Keep the original transcript instead of deleting them",
                          metavar="",)
    
        self.add_argument("-t","--tmpdir",
                          help="Temporary directory",
                          metavar="",
                          default=f"tmp_{datetime.datetime.now()}")
        self.add_argument("-p", "--processes",
                          help="Number of parallel processes to use [Default:4]",
                          type=int,
                          metavar="",
                          default=4)
        self.add_argument("-l","--log_level",
                          help="[default: info]",
                          default="info",
                          metavar="")
               
    def write_parameter_file(self):

        with open(os.path.join(self.args.outdir, "utrpy.param"), "w") as param_file:
            
            param_file.write(f"gff_prediction   {self.args.gff_prediction:<15}\n")
            param_file.write(f"gff_assembly     {self.args.gff_assembly:<15}\n")
            param_file.write(f"--match          {self.args.match:<15}\n")
            param_file.write(f"--select         {self.args.select:<15}\n")
            param_file.write(f"--know_strand    {self.args.know_strand:<15}\n")
            param_file.write(f"--processes      {self.args.processes:<15}\n")
            if self.args.pinky_promise:
                param_file.write("--pinky_promise")

    def check_input(self):

        if not os.path.isfile(self.args.gff_prediction):
            logging.error(f"{self.args.gff_prediction} is not a file")
            exit(1)
        if not os.path.isfile(self.args.gff_assembly):
            logging.error(f"{self.args.gff_assembly} is not a file")
            exit(1)
        if os.path.isdir(self.args.outdir):
            logging.error(f"{self.args.outdir} is exists")
            exit(1)

    def parse_args(self):

        self.args = super().parse_args()

        self.check_input()
        os.mkdir(self.args.outdir)
        os.mkdir(self.args.tmpdir)
        self.write_parameter_file()  
        
        return self.args
