#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from pathlib import Path

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
pantax-dbg profile [-h] -r [...] [--single]--db-prefix DB_PREFIX -i REF_INFO -o [-t] [-k]
"""

def subcmd_build_custom(argv):
    if any(a in ("-h", "--help") for a in argv):
        print(BUILD_CUSTOM_HELP)
        sys.exit(0)
    cmd = [_ganon(), "build-custom"] + list(argv)
    run_cmd(cmd, echo=False)


def subcmd_profile(args):
    if not args.reads:
        raise SystemExit("[PanTax-DBG] --reads/-r .Paired: -r R1 -r R2；Single: -r R.")

    profile_mod.run(
        reads=args.reads,
        single=args.single,
        db_prefix=args.db_prefix,
        out_prefix=args.out,
        report_type="abundance",
        file_info=args.ref_info,
        threads=args.threads,
        k=args.kmer,
        # new knobs
        strain_min_reads=args.strain_min_reads,
        strain_topk=args.strain_topk,
        #strain_keep_ties=args.strain_keep_ties,
        strain_div=args.strain_div,

        r1_large_n=args.r1_large_n,
        r1_topk_large=args.r1_topk_large,
        r1_min_abundance=args.r1_min_abundance,
        species_min_abundance=args.species_min_abundance,
        # ganon_singleton_min_abund=args.ganon_singleton_min_abund,
    )


def build_parser():
    p = argparse.ArgumentParser(
        prog="pantax-dbg",
        description="PanTax-DBG: a robust and accurate species-level metagenomic profiler."
    )
    sub = p.add_subparsers(dest="subcmd", metavar="<command>")

   
    p_build = sub.add_parser(
        "build-custom",
        help="Build custom PanTax-DBG databases."
    )
    p_build.add_argument("ganon_args", nargs=argparse.REMAINDER,
                         help="Arguments passed directly to 'ganon build-cutstom'.")


    #p_prof = sub.add_parser("profile", help="Profile reads against custom databases.")
    p_prof = sub.add_parser(
        "profile",
        help="Profile reads against custom databases.",
        usage=PROFILE_USAGE,                         
    )
    p_prof.add_argument("-r", "--reads", action="append", required=True,metavar="",
                        help=("For paired-end data, specify mates consecutively: -r R1.fq -r R2.fq. "
                        "For single-end data, use --single and give one -r per file. "))
    p_prof.add_argument("--single", action="store_true", help="Treat input as single-end reads. ")
    p_prof.add_argument("-d", "--db-prefix", required=True,metavar="",
                        help="Database input prefix.")
    p_prof.add_argument("-i", "--ref-info", required=True,metavar="",
                        help=("Tab-separated reference metadata file. Fields: "
                        "strain_name <tab> strain_taxid <tab> species_taxid "
                        "<tab> species_name <tab> genome_path. "
                        "strain_name and strain_taxid must be unique."))
    p_prof.add_argument("-o", "--output", required=True, metavar="", dest="out", help="Output directory for profiling results.")
    p_prof.add_argument("-t", "--threads", type=int, default=8, metavar="", help="Number of threads.")
    p_prof.add_argument("-k", "--kmer-size", dest="kmer",type=int, default=31, metavar="", help="k-mer size used in the ccDBG-based profiling step.")



    # ---- strain postprocess params ----
    p_prof.add_argument("--strain-min-reads", type=float, default=5.0,
                        help="Minimum assigned_reads to keep a strain-group row before TopK.")
    p_prof.add_argument("--strain-topk", type=int, default=5,
                        help="Keep top-K strain groups per species (by rel_abundance_within_species).")
    # p_prof.add_argument("--strain-keep-ties", action="store_true", default=False,
    #                     help="Keep all ties at K-th score when selecting TopK strain groups.")
    p_prof.add_argument("--strain-div", type=float, default=5.0,
                        help="cutoff = threshold / div for further strain filtering (default 5.0).")

    # ---- ganon mapping switch params ----
    p_prof.add_argument("--r1-large-n", type=int, default=1000,
                        help="If r1 species hits > this value, use large-n TopK+singleton filtering mapping.")
    p_prof.add_argument("--r1-topk-large", type=int, default=3,
                        help="TopK strains used in ccDBG.")
    # p_prof.add_argument("--ganon-singleton-min-abund", type=float, default=7e-7,
    #                     help="singleton_min_abund used in large-n mapping step.")
    p_prof.add_argument("--r1-min-abundance", type=float, default=0.0,
                        help="Minimum relative abundance (from Ganon) to keep a species for ccDBG construction (default: 0.0).")
    p_prof.add_argument("--species-min-abundance", type=float, default=0.0,
                        help="Minimum relative abundance for final species output. Species below this are removed and abundances are re-normalized. (default: 0.0).")
    return p


def main():
    p = build_parser()
    argv = sys.argv[1:]

    
    if not argv:
        p.print_help()
        sys.exit(0)

    
    if argv[0] in ("-h", "--help"):
        p.print_help()
        sys.exit(0)

    
    subcmd = argv[0]

    if subcmd == "build-custom":
        
        subcmd_build_custom(argv[1:])
    elif subcmd == "profile":
        
        args = p.parse_args(argv)
        subcmd_profile(args)
    else:
        
        p.print_help()
        sys.exit(1)


def build_custom_main():
    subcmd_build_custom(sys.argv[1:])

def profile_main():
    p = build_parser()
    fake_argv = ["themis", "profile"] + sys.argv[1:]
    args = p.parse_args(fake_argv[1:])
    subcmd_profile(args)
