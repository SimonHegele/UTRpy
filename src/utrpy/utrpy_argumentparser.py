from argparse import ArgumentParser

class UTRpyArgparser(ArgumentParser):

    prog        =   "utrpy"

    description = ("UTR extension of transcript exons from protein orthology based gene "
                   "prediction using exons from reference based assembly")
                       
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        # Input files
        self.add_argument("gff_prediction",
                          help="GFF-format genome annotation from gene prediction.")
        self.add_argument("gff_assembly",
                          help="GFF-format genome annotation from transcriptome assembly.")
        self.add_argument("gff_utrpy",
                          help="Output file path")
        
        # Options
        self.add_argument("-max","--maximum_exon_length",
                          type=int,
                          default=20_000,
                          help="Maximum exon length prevents the use of unreasonably long exons. [Default:20000]")
        self.add_argument("-t", "--threads",
                          type=int,
                          default=4,
                          help="Number of threads to use [Default:4]")
        self.add_argument("-s","--strict_strandedness",
                          action="store_true",
                          default=False,
                          help="Exons with unknown strandedness are not used")
        self.add_argument("-e","--explicit",
                          action="store_true",
                          default=False,
                          help="Explicitly add UTRs as features")
