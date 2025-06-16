import logging
import multiprocessing
import pandas
import os
import shutil
import subprocess

from .utrpy_agat_prepare     import agat_prepare
from .utrpy_argumentparser   import UTRpyArgparser
from .utrpy_gff_utils        import load_gff, seqname_split, write_gff
from .utrpy_logging          import logging_setup
from .utrpy_utr_extend       import utr_extend_threaded

def main():

    logging_setup(args)

    args = UTRpyArgparser().parse_args()

    logging.info("Preprocessing with AGAT")
    agat_prepare(args)

    p_gff = load_gff(os.path.join(args.tmpdir, "prediction.gff"))
    a_gff = load_gff(os.path.join(args.tmpdir, "assembly.gff"))

    match args.match:
        case "ends":
            match_middle_exons = False
        case "all":
            match_middle_exons = True

    seqnames = list(set(p_gff["seqname"].unique()) &
                    set(a_gff["seqname"].unique()))
    
    p_gff = seqname_split(p_gff)
    a_gff = seqname_split(a_gff)

    mp_args = zip([p_gff[s]             for s in seqnames],
                  [a_gff[s]             for s in seqnames],
                  [match_middle_exons   for _ in seqnames],
                  [args.know_strand     for _ in seqnames],
                  [args.keep            for _ in seqnames],
                  [args.select          for _ in seqnames],
                  [args.max_exon_length for _ in seqnames])
    
    with multiprocessing.Pool(args.processes) as pool:
        u_gff = pool.map(utr_extend_threaded, mp_args)

    u_gff = pandas.concat(u_gff)

    write_gff(u_gff, os.path.join(args.tmpdir, "utrpy.gff"))

    subprocess.run(["agat_convert_sp_gxf2gxf.pl",
                    "--gff", os.path.join(args.tmpdir, "utrpy.gff"),
                    "-o", os.path.join(args.outdir, "utrpy.gff")],
                    check=True)
    
    shutil.rmtree(args.tmpdir)

    logging.info("#############################################")
    logging.info("#    Simon says: Thanks for using UTRpy!    #")
    logging.info("#############################################")
