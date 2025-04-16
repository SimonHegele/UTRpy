"""
Module Name:    utrpy_agat_prepare
Description:    Provides method agat_prepare which calls AGAT to fix inconsistencies in
                the input files and adds missing features that are only implicitly given
                in the atrributes.
Author:         Simon Hegele
Date:           2025-04-01
Version:        1.0
License:        GPL-3
"""

import argparse
import os
import subprocess

def agat_prepare(args: argparse.Namespace) -> None:

    if args.gff_assembly.endswith(".gtf"):
        subprocess.run(["agat_convert_sp_gxf2gxf.pl",
                        "--gtf", args.gff_assembly,
                        "-o", os.path.join(args.tmpdir, "assembly.gff")],
                        check=True)
    if args.gff_assembly.endswith(".gff"):
        subprocess.run(["agat_convert_sp_gxf2gxf.pl",
                        "--gff", args.gff_assembly,
                        "-o", os.path.join(args.tmpdir, "assembly.gff")],
                        check=True)

    if args.pinky_promise:
        subprocess.run(["cp",
                        args.gff_prediction,
                        os.path.join(args.tmpdir, "prediction.gff")],
                        check=True)
    else:
        if args.gff_assembly.endswith(".gtf"):
            subprocess.run(["agat_convert_sp_gxf2gxf.pl",
                            "--gtf", args.gff_prediction,
                            "-o", os.path.join(args.tmpdir, "prediction.gff")],
                            check=True)
        if args.gff_assembly.endswith(".gff"):
            subprocess.run(["agat_convert_sp_gxf2gxf.pl",
                            "--gff", args.gff_prediction,
                            "-o", os.path.join(args.tmpdir, "prediction.gff")],
                            check=True)  