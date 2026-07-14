#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse

from . import profile as profile_mod
from .utils import run_cmd
from .paths import get_dbg_ganon


def _ganon():
    return os.environ.get("PANTAX_DBG_GANON_BIN") or get_dbg_ganon()


BUILD_CUSTOM_HELP = """\
usage: pantax-dbg build-custom [-h] [-i [...]] [-e] [-c] [-n] [-a] [-l] [-m [...]] [-z [...]] [--skip-genome-size] [-r [...]]
                          [-q [...]] -d DB_PREFIX [-x] [-t] [-p] [-k] [-w] [-s] [-f] [-j] [-y] [-v] [--restart]
                          [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -i [ ...], --input [ ...]
                        Input file(s) and/or folder(s). Mutually exclusive --input-file. (default: None)
  -e , --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: fna.gz)
  -c, --input-recursive
                        Look for files recursively in folder(s) provided with --input (default: False)
  -d DB_PREFIX, --db-prefix DB_PREFIX
                        Database output prefix (default: None)

custom arguments:
  -n , --input-file     Tab-separated file with all necessary file/sequence information. Fields: file [<tab> target
                        <tab> node <tab> specialization <tab> specialization name]. Mutually exclusive --input (default: None)
  -a , --input-target   Target to use [file, sequence]. Parse input by file or by sequence. Using 'file' is recommended
                        and will speed-up the building process (default: file)
  -l , --level          Max. level to build the database. By default, --level is the --input-target. Options: any
                        available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy). Further
                        specialization options [assembly, custom]. assembly will retrieve and use the assembly accession
                        and name. custom requires and uses the specialization field in the --input-file. (default: None)
  -m [ ...], --taxonomy-files [ ...]
                        Specific files for taxonomy - otherwise files will be downloaded (default: None)
  -z [ ...], --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)

ncbi arguments:
  -r [ ...], --ncbi-sequence-info [ ...]
                        Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information.
                        [eutils, nucl_gb, nucl_wgs, nucl_est, nucl_gss, pdb, prot, dead_nucl, dead_wgs, dead_prot or one
                        or more accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/].
                        By default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise. (default:
                        [])
  -q [ ...], --ncbi-file-info [ ...]
                        Downloads assembly_summary files to extract target information. [refseq, genbank,
                        refseq_historical, genbank_historical or one or more assembly_summary files from
                        https://ftp.ncbi.nlm.nih.gov/genomes/] (default: ['refseq', 'genbank'])

important arguments:
  -x , --taxonomy       Set taxonomy to enable taxonomic classification, lca and reports [ncbi, gtdb, skip] (default:
                        ncbi)
  -t , --threads 

advanced arguments:
  -p , --max-fp         Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k , --kmer-size      The k-mer size to split sequences. (default: 19)
  -w , --window-size    The window-size to build filter with minimizers. (default: 31)
  -s , --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f , --filter-size    Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j , --mode           Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive rate. By default, an average value is calculated
                        to balance classification speed and database size. Only valid for --filter-type ibf. (default:
                        avg)
  -y , --min-length     Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v , --filter-type    Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
                        path set with --raptor-path. --mode, --filter-size and --min-length will be ignored with hibf.
                        hibf will set --max-fp 0.001 as default. (default: hibf)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
"""

PROFILE_USAGE = """\
pantax-dbg profile [-h] -r READ1 -r READ2 -d PREFIX -i FILE -o DIR [options]
"""


def subcmd_build_custom(argv):
    if any(a in ("-h", "--help") for a in argv):
        print(BUILD_CUSTOM_HELP)
        sys.exit(0)
    cmd = [_ganon(), "build-custom"] + list(argv)
    run_cmd(cmd, echo=False)


def subcmd_profile(args):
    if len(args.reads) != 2:
        raise SystemExit(
            "[PanTax-DBG][error] Expected two paired-end short-read FASTQ files. "
            "Please provide: -r read1.fq.gz -r read2.fq.gz"
        )

    profile_mod.run(
        reads=args.reads,
        single=False,
        db_prefix=args.db_prefix,
        out_prefix=args.out,
        report_type="abundance",
        file_info=args.ref_info,
        threads=args.threads,
        k=args.kmer,
        strain_min_reads=args.strain_min_reads,
        strain_topk=args.strain_topk,
        strain_div=args.strain_div,
        r1_large_n=args.r1_large_n,
        r1_topk_large=args.r1_topk_large,
        r1_min_abundance=args.r1_min_abundance,
        species_min_abundance=args.species_min_abundance,
        ggcat_temp_dir=args.ggcat_temp_dir,
        ggcat_memory=None,
        ggcat_prefer_memory=False,
        ggcat_disk_optimization_level=None,
        ggcat_intermediate_compression_level=None,
    )


def build_parser():
    parser = argparse.ArgumentParser(
        prog="pantax-dbg",
        description=(
            "PanTax-DBG: accurate species- and strain-level profiling "
            "for paired-end short-read metagenomes."
        ),
    )
    subparsers = parser.add_subparsers(dest="subcmd", metavar="<command>")

    build_parser = subparsers.add_parser(
        "build-custom",
        help="Build a custom PanTax-DBG reference database.",
    )
    build_parser.add_argument(
        "ganon_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to the PanTax-DBG database builder.",
    )

    profile_parser = subparsers.add_parser(
        "profile",
        help="Profile paired-end reads at species and strain levels.",
        description=(
            "Estimate species- and strain-level relative abundances from "
            "paired-end short-read metagenomes using PanTax-DBG."
        ),
        usage=PROFILE_USAGE,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    required = profile_parser.add_argument_group("required arguments")
    required.add_argument(
        "-r", "--reads",
        action="append",
        required=True,
        metavar="FASTQ",
        help=(
            "Input paired-end FASTQ file. Supply exactly twice in mate order:\n"
            "  -r read1.fq.gz -r read2.fq.gz"
        ),
    )
    required.add_argument(
        "-d", "--db-prefix",
        required=True,
        metavar="PREFIX",
        help="Prefix of a database built with 'pantax-dbg build-custom'.",
    )
    required.add_argument(
        "-i", "--ref-info",
        required=True,
        metavar="FILE",
        help=(
            "Tab-delimited reference metadata table for DBG refinement.\n"
            "Columns: strain_name, strain_taxid, species_taxid, species_name, genome_path."
        ),
    )
    required.add_argument(
        "-o", "--output", "--out",
        dest="out",
        required=True,
        metavar="DIR",
        help="Output directory for final species- and strain-level profiles.",
    )

    profiling = profile_parser.add_argument_group("profiling arguments")
    profiling.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        metavar="INT",
        help="Number of worker threads. (default: 8)",
    )
    profiling.add_argument(
        "-k", "--kmer-size",
        dest="kmer",
        type=int,
        default=31,
        metavar="INT",
        help="k-mer size for colored DBG refinement. (default: 31)",
    )
    profiling.add_argument(
        "--species-min-abundance",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help=(
            "Minimum relative abundance retained in the final species profile;\n"
            "retained abundances are renormalized. (default: 0.0)"
        ),
    )

    strain = profile_parser.add_argument_group("strain refinement arguments")
    strain.add_argument(
        "--strain-min-reads",
        type=float,
        default=5.0,
        metavar="FLOAT",
        help=(
            "Minimum assigned-read support required before strain candidate\n"
            "ranking. (default: 5.0)"
        ),
    )
    strain.add_argument(
        "--strain-topk",
        type=int,
        default=5,
        metavar="INT",
        help=(
            "Maximum number of retained strain candidates per species after\n"
            "ranking. (default: 5)"
        ),
    )
    strain.add_argument(
        "--strain-div",
        type=float,
        default=5.0,
        metavar="FLOAT",
        help=(
            "Divisor used in secondary strain-support filtering.\n"
            "(default: 5.0)"
        ),
    )

    candidates = profile_parser.add_argument_group("candidate screening arguments")
    candidates.add_argument(
        "--r1-large-n",
        type=int,
        default=1000,
        metavar="INT",
        help=(
            "Candidate-set size above which large-set top-k filtering is\n"
            "applied. (default: 1000)"
        ),
    )
    candidates.add_argument(
        "--r1-topk-large",
        type=int,
        default=3,
        metavar="INT",
        help=(
            "Number of strain candidates retained per species in large-set\n"
            "filtering mode. (default: 3)"
        ),
    )
    candidates.add_argument(
        "--r1-min-abundance",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help=(
            "Minimum first-round species abundance retained for DBG\n"
            "refinement. (default: 0.0)"
        ),
    )

    technical = profile_parser.add_argument_group("technical arguments")
    technical.add_argument(
        "--ggcat-temp-dir",
        default=None,
        metavar="DIR",
        help=(
            "Directory for DBG-ggcat temporary build files. By default, a\n"
            "run-specific directory is created under $TMPDIR or /tmp."
        ),
    )
    return parser


def main():
    parser = build_parser()
    argv = sys.argv[1:]

    if not argv or argv[0] in ("-h", "--help"):
        parser.print_help()
        sys.exit(0)

    subcmd = argv[0]
    if subcmd == "build-custom":
        subcmd_build_custom(argv[1:])
    elif subcmd == "profile":
        args = parser.parse_args(argv)
        subcmd_profile(args)
    else:
        parser.print_help()
        sys.exit(1)


def build_custom_main():
    subcmd_build_custom(sys.argv[1:])


def profile_main():
    parser = build_parser()
    args = parser.parse_args(["profile"] + sys.argv[1:])
    subcmd_profile(args)
