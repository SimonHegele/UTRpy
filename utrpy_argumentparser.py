"""
Script Name:    utrpy_argumentparser.py
Description:    Provides UTRpyArgparser(ArgumentParser)
                - Required arguments added on initialization
                - Extended parse_args() to check input and write a file with the parameters
Author:         Simon Hegele
Date:           2025-04-01
Version:        0.2
License:        GPL-3
"""

import logging

from argparse import ArgumentParser
from os.path  import isfile

class UTRpyArgparser(ArgumentParser):

    prog        =  "utrpy"

    description = ("UTR extension of transcript exons from protein orthology based gene "
                   "prediction using exons from reference based assembly")
                       
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        # Input files
        self.add_argument("gff_prediction",
                          help="Annotation from gene prediction (GFF)")
        self.add_argument("gff_assembly",
                          help="Annotation from transcriptome assembly (GFF/GTF)")
        self.add_argument("gff_utrpy",
                          help="Output file")
        
        # Exon extension
        self.add_argument("-a","--ambiguities",
                          help="How to chose UTR extension when there are multiple [Default:smallest]",
                          choices=["none", "smallest", "longest"],
                          default="smallest")
        self.add_argument("-k","--know_strand",
                          action="store_true",
                          help="Exons with unknown strandedness are not used")
        self.add_argument("-m","--max_ex_len",
                          type=int,
                          metavar="",
                          default=20_000,
                          help="Maximum allowed exon length [Default:20000]")
        
        # UTR annotation
        self.add_argument("-e","--explicit",
                          action="store_true",
                          help="Explicitly add UTRs as features")
    
        # Others
        self.add_argument("-t", "--threads",
                          type=int,
                          metavar="",
                          default=4,
                          help="[Default:4]")
        
        
    def write_parameter_file(self):

        with open(f"{self.args.gff_utrpy}.param", "w") as param_file:
            
            # Files
            param_file.write("{0:<15} {1:<15}\n".format("gff_prediction",   self.args.gff_prediction))
            param_file.write("{0:<15} {1:<15}\n".format("gff_assembly",     self.args.gff_assembly))
            param_file.write("{0:<15} {1:<15}\n".format("gff_utrpy",        self.args.gff_utrpy))
            param_file.write("{0:<15} {1:<15}\n".format("--know_strand",    self.args.know_strand))
            param_file.write("{0:<15} {1:<15}\n".format("--explicit",       self.args.explicit))
            param_file.write("{0:<15} {1:<15}\n".format("--max_ex_len",     self.args.max_ex_len))
            param_file.write("{0:<15} {1:<15}\n".format("--threads",        self.args.threads))

    def check_input(self):

        # Checking give file-paths
        if not isfile(self.args.gff_prediction):
            raise FileNotFoundError(f"{self.args.gff_prediction}")
        if not isfile(self.args.gff_assembly):
            raise FileNotFoundError(f"{self.args.gff_assembly}")
        if isfile(self.args.gff_utrpy):
            print("File exists")
            raise FileExistsError(f"{self.args.gff_utrpy}")

    def parse_args(self):

        self.args = super().parse_args()

        self.check_input()
        self.write_parameter_file()  
        
        return self.args