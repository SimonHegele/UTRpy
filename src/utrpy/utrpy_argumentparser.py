from argparse import ArgumentParser

class UTRpyArgparser(ArgumentParser):

    prog        =   "utrpy"

    description = ( "Genome annotions from protein orthology based gene prediction tools"
                    "lack UTRs. Exons in these kind of gene annotations can be extended to"
                    "also cover UTRs using exons generated with reference based assembly"
                    "tools such as StringTie." )
                       
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
        self.add_argument("-mel","--maximum_exon_length",
                          type=int,
                          default=20_000,
                          help="Maximum exon length prevents the use of unreasonably long exons. [Default:20000]")
        self.add_argument("-meo","--minimum_exon_overlap",
                          type=int,
                          help="Minimum overlap of exons to be considered [Default:Auto]")
        self.add_argument("-t", "--threads",
                          type=int,
                          default=4,
                          help="Number of threads to use [Default:4]")
        self.add_argument("-c","--conservative",
                          action="store_true",
                          default=True,
                          help="Will not use exons if their strand is not known [Default:True]")
