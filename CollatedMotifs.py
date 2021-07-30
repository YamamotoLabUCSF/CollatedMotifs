#!/usr/local/bin/anaconda3/bin/python3
# Note: edit shebang line above as appropriate for your system
# AUTH: Kirk Ehmsen
# FILE: CollatedMotifs.py
# DATE: 09-04-2018/08-04-2019 & 3-02-2021/06-26-2021/07-01-2021
# DESC: This script accepts text to standard input, and returns TFBS motifs for samples
# from a demultiplexed NGS fastq dataset.  Provided with a reference sequence, the script returns TFBS motifs
# collated as 'new' or 'lost' relative to the reference sequence.
# USAGE: ./CollatedMotifs.py or python3 CollatedMotifs.py
# REPO: https://github.com/YamamotoLabUCSF/CollatedMotifs

#############################################################################
# Background notes:
# ==============================================
# This is CollatedMotifs.py v1.1
# ==============================================
# https://github.com/YamamotoLabUCSF/CollatedMotifs
# v1.1/Committed 7-01-2021
# ----------------------------------------------
# For usage details, please refer to README file at GitHub location and to the following manuscript:
#   Ehmsen, Knuesel, Asahina, Aridomi, Yamamoto (2021)
# Please cite usage as:
#   CollatedMotifs.py
#   Ehmsen, Knuesel, Asahina, Aridomi, Yamamoto (2021)

# Operation notes:
# ==============================================
# This script accepts text to standard input, and for samples from a demultiplexed NGS fastq dataset,
# (1) uses BLASTN to produce alignments with user-provided reference sequences, (2) uses FIMO
# to identify matches to TFBS motifs present in read and reference sequences, and (3) presents a summary
# comparison of TFBS shared between read and reference, or gained/lost in read vs. reference.
# Python3, BLASTN, MAKEBLASTDB, FIMO, and FASTA-GET-MARKOV are required.
#
# What does this script do?
# 1. counts reads per well/sample (sample defined by fastq file); fastq file name provides the sample name
# 2. identifies the top 5 reads per well (in terms of read abundance), and calculates representation among reads within
# the well at four levels:
#   (a) raw frequency (% read type in question, relative to total reads)
#   (b) percentile (% of other read types that fall below the frequency of the read type in question)
#   (c) adjusted frequency @ 1% (% read type in question, relative to reads that occur at >1% frequency)
#   (d) adjusted frequency @ 10% (% read type in question, relative to reads that occur at >10% frequency)
# 3. for user-provided reference sequences, uses FIMO (MEME suite, http://meme-suite.org/tools/fimo) and
#    user-provided positional weight matrix file to find matches to TFBS motifs
# 4. for top 5 reads, uses FIMO and user-provided positional weight matrix file to find matches to TFBS motifs
# 5. compares TFBS in reads to TFBS in specific reference sequence, outputting 'new' and 'lost' TFBS relative to the reference sequence.
#
# As part of these operations, a BLAST alignment database is created by MAKEBLASTDB from a user-supplied reference sequence file,
# and a Markov background file is created for FIMO by FASTA-GET-MARKOV from a user-supplied reference sequence file.

# Input notes:
# ==============================================
# You will be prompted for the following user-specific information (11 items--10 required, 1 optional):
#      ** required (4 directory paths):
#         * where should output files go?                     path to output directory for output files
#         * where are input files found?                      path to single directory containing
#                                                               demultiplexed fastq files
#         * where is BLASTN executable found?                 path to BLASTN installation
#         * where is the reference sequence                   path to directory containing six files
#             database used for alignment?                      that compose the reference sequence database
#                                                               reference for BLASTN alignments
#                                                               (.nhr, .nin, .nog, .nsd, .nsi, .nsg)
#
#      ** required:
#         2 directory paths
#           * where should output files go?                     path to output directory for output files
#           * where are input files found?                      path to single directory containing
#                                                                 demultiplexed fastq files
#         4 executable paths
#           * where is BLASTN executable found?                 path to BLASTN installation
#           * where is the FIMO executable found?               path to FIMO installation
#           * where is the MAKEBLASTDB executable found?        path to MAKEBLASTDB installation
#           * where is the FASTA-GET-MARKOV executable found?   path to FASTA-GET-MARKOV installation
#         3 file paths
#           * what are your reference sequence(s),              path to single fasta file, containing reference
#             to which you will (a) align sequenced               sequence(s) for processing by
#             reads, and (b) compare sequenced reads              (a) MAKEBLASTDB, to generate a database
#             for TFBS occurrence?                                    reference for BLASTN, and
#                                                                 (b) FIMO, to establish TFBS occurrence(s)
#                                                                     to be evaluated relative to sequenced
#                                                                     reads
#           * what are the TFBS motif(s)                        path to single text file, containing
#             for which you will search,                          position frequency matrix(ces) for TFs
#             and for which you will draw comparisons             
#             for presence/absence between sequences?
#           * what DNA sequence(s) will you use as a basis      path to single text/fasta file, containing
#             for Markov background estimation, to be used        DNA sequence(s) from which a Markov background
#             by FIMO                                             file will be generated for use by FIMO
#
#      ** optional:
#           * name of transcription factor (TF) of interest,    single string with TF name
#             selected from standardized Entrez gene name
#             options available in the positional weight
#             matrix file you provide

# Output notes:
# ==============================================
# This script produces 6 output files in the user-specified output directory, plus three directories:
# two directories and subsidiary files created by FIMO (fimo_out and fimo_out_ref) and one directory
# and subsidiary files created by MAKEBLASTDB (alignment_database).
# CollatedMotifs.py output files include: 
#	  1. fasta.fa
#	  2. blastn_alignments.txt (output of BLASTN operation on fasta.fa)
#     3. markov_background.txt (output of FASTA-GET-MARKOV operation on user-supplied fasta reference file)
#     4. collated_TFBS.txt (output of script operation on FIMO-generated .tsv files in fimo_out and fimo_out_ref)
#     5. collated_TFBS.xlsx (output of script operation on FIMO-generated data, compiling TFBS lost/gained in
#        sample alleles relative to reference sequence; includes special focus on TFBS for TF of interest,
#        if optional 'TF of interest' is provided during user input)
#     6. script_metrics.txt (summary/analysis of script operation metrics [metadata])
#
#           Directory structure under an output directory specified as 'CollatedMotifs', for example,
#           would contain the following subdirectories and files following CollatedMotifs.py operations:
#
#           /CollatedMotifs 
#                          `-----/alignment_database
#                                        `----------*.nin
#                                        `----------*.nhr
#                                        `----------*.nog
#                                        `----------*.nsd
#                                        `----------*.nsg
#                                        `----------*.nsi
#                          `-----blastn_alignments.txt
#                          `-----collated_TFBS.txt
#                          `-----collated_TFBS.xlsx
#                          `-----fasta.fa
#                          `-----/fimo_out
#                                        `----------cisml.xml
#                                        `----------fimo.gff
#                                        `----------fimo.html
#                                        `----------fimo.tsv
#                                        `----------fimo.xml
#                          `-----/fimo_out_ref
#                                        `----------cisml.xml
#                                        `----------fimo.gff
#                                        `----------fimo.html
#                                        `----------fimo.tsv
#                                        `----------fimo.xml
#                          `-----markov_background.txt
#                          `-----script_metrics.txt
#############################################################################

#############################################################################
# SCRIPT:

# Check for availability of Python dependencies (libraries, modules) in path
missing_dependencies_list = []

try:
    import psutil
except ImportError:
    missing_dependencies_list.append('psutil')
    
try:
    import numpy
except ImportError:
    missing_dependencies_list.append('numpy')

try:
    import scipy
except ImportError:
    missing_dependencies_list.append('scipy')
    
if len(missing_dependencies_list) > 0:
    print('ModuleNotFoundError\n')
    print('Please note, the following required Python module(s) are not found in your Python system path:')
    for i in missing_dependencies_list:
        print('   '+i)
    print('\nPlease exit the script and install these Python dependencies in your system path.')
    print("""\nGuidelines for installation of Python dependencies can be found in the README file for
    CollatedMotifs.py ('System Setup')""")
    print("""    (Creation of a Python virtual environment is recommended)""")

# Import libraries, modules
# Operating system interfaces
import os

# Time access and conversions, Basic data and time types
import time
from datetime import datetime

# System-specific parameters and functions
import sys

# Process and system utilities
import psutil
from psutil import virtual_memory

# Gzip to read GNU zipped files
import gzip

# Low-level networking interface
import socket

# System version information
import platform

# Unix-style pathname pattern expansion
import glob

# NumPy (numeric operations)
import numpy

# SciPy (for percentile) 
from scipy import stats

# Container datatypes (for Counter operation)
from collections import Counter

# Decimal fixed point and floating point arithmetic
from decimal import Decimal

# Regular expression operations
import re

# Object-oriented filesystem paths
from pathlib import Path

# Internationalization services (for use of thousands separator in numbers where appropriate)
import locale
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

# pandas for panel data
import pandas as pd

# start time
initialTime = datetime.now()

# Define 'prompts' function for coached user input
def prompts():
    """Coached prompts to collect user input"""
    # Make variables assigned in prompts() function globally available
    global output_directory
    global fastq_directory
    global fasta_ref
    global blastn_path
    global makeblastdb_path
    global db_prefix
    global fimo_path
    global fimo_motifs_path
    global fasta_get_markov_path
    global markov_background_file
    global TF_of_interest
    # 1-Specify output directory.
    print(r"""
    ---------------------------------------------
    Location of OUTPUT DIRECTORY for output files
    ---------------------------------------------
    
    This script produces 6 output files in the user-specified output directory, plus three directories:
    two directories and subsidiary files created by FIMO (fimo_out and fimo_out_ref) and one directory
    and subsidiary files created by MAKEBLASTDB (alignment_database).
    
    CollatedMotifs.py output files include:
    
        1. fasta.fa

        2. blastn_alignments.txt
            (output of BLASTN operation on fasta.fa)

        3. markov_background.txt
            (output of FASTA-GET-MARKOV operation on user-supplied fasta reference file)

        4. collated_TFBS.txt
            (output of script operation on FIMO-generated .tsv files in fimo_out and fimo_out_ref)
            
        5. collated_TFBS.xlsx
            (output of script interpretation of lost and gained TFBS, detailed for inferred alleles in spreadsheet)
            
        6. script_metrics.txt (summary/analysis of script operation metrics [metadata])
  
            Note: 
            * These files do not exist before the script is run. The files are made by the script.
            * The primary data outputs for TFBS comparisons are found in collated_TFBS.txt
        
    At this prompt, indicate an absolute path to a ** directory ** that will be created by the script as the location
    for output files.  This directory should not exist yet -- it will be created as an output of this script, and will
    be populated with the file outputs of this specific instance of the script operation.

    Use only forward slashes ('/') as directory separators.

    Example: if you'd like to create a directory ('CollatedMotifs') in an existing directory ('Illumina'), accessed
    with absolute path of '/Users/myname/Illumina/CollatedMotifs' (Mac) or 'C:\Users\myname\Illumina\CollatedMotifs'
    (Windows), enter '/Users/myname/Illumina/CollatedMotifs' at the command line prompt. Replace 'myname' with the
    appropriate intervening directory identifiers. Do *not* flank your entry with quotation marks (') at the
    command-line.
    
    Alternatively, simply enter a desired directory name (e.g., 'CollatedMotifs') and run this script from
    within a directory where you'd like to create this new directory."""+'\n')
    output_directory = input(r"""    -----> Output directory name and path:  """)
    # 2-Specify the fastq files to be used for input, by indicating directory location of the file list.
    print(r"""
    ------------------------------------------------------------------------------
    Location of INPUT FILES (single directory containing demutiplexed fastq files)
    ------------------------------------------------------------------------------

    You will now be asked to enter the path to the directory containing the fastq files
    to be processed as CollatedMotifs.py input.  

    Example: if your fastq input files are named file1.fastq, file2.fastq, etc. and are found in a directory
    named 'Sequences' with absolute path of '/Users/myname/Sequences' (Mac) or 'C:\Users\myname\Sequences' (PC),
    enter '/Users/myname/Sequences' at the command line prompt.

    When you're done entering the fastq file location, press 'Enter' again to proceed in the script."""+'\n')
    fastq_directory = input(r"""    -----> Directory name and path:  """)
    # 3-Specify fasta file containing reference sequences as basis for TFBS motif comparisons/contrasts.
    print(r"""
    -----------------------------------------
    Location of FIMO REFERENCE SEQUENCES FILE
    -----------------------------------------

    This script aligns and compares your top sample read sequence(s) to a defined reference sequence,
    as its basis for determining distinct vs. common TFBS motifs.
    Please indicate the absolute path to a fasta file containing reference sequences.

    Example: if you have three sample names relating to three different reference sequences, enter these
    sequences in fasta format, saved in a single text file. Each fasta entry definition line (defline)
    should be named such that the defline name matches a unique descriptor in fastq file names.
    
        >Sample1
        GATCGACTAGAGCGAGCATTCATCATATCACGAGTAGCATCGACGTGCACGATCGATCGTAGCTAGCTAGTCATGCATGCATGCTAGATTCGAGCATGCATGCTAC
        >Sample2
        AGTAGCTGTGATGCTAGTCATCTAGCTAGCAGCGTAGCTAGCGATCGATCTAGAGCCGATCGATCGAGCATCTAGCTATCAGCGGCGGGATCATCTATCTACGGG
        >Sample3
        CGATGCAGCGCGATCGAGCGCGATCGATATTAGCATGCGCAGCTAGCTAGCTGGCGATCGATGCATGCTAGCTGTGTCAGTCGACGATCACACGATCACACTGTGTG

    When you're done entering the list of reference sequences, press 'Enter' again to proceed in the script."""+'\n')
    fasta_ref = input(r"""    -----> Path to fasta file containing reference sequences:  """)
    # 4-Collect path to blastn executable.
    print(r"""
    -----------------------------
    Location of BLASTN EXECUTABLE
    -----------------------------

    This script uses BLASTN (NCBI) to align reads from your fastq files to a reference sequence database.
    Please indicate the absolute path to the BLASTN executable.

    Example: if your BLASTN executable is found at absolute path /Users/myname/blastn, type '/Users/myname/blastn'
    and press Enter."""+'\n')
    blastn_path = input(r"""    -----> Path to BLASTN executable:  """)
    # 5-Collect path to makeblastdb executable.
    print(r"""
    ----------------------------------
	Location of MAKEBLASTDB EXECUTABLE
    ----------------------------------

    Because this script uses BLASTN (NCBI) to align reads from your fastq files to a reference sequence database,
    a compatible reference sequence database is required. This script uses MAKEBLASTDB (NCBI) to generate
    a reference sequence database from the reference sequences in the fasta file you provided earlier.
    
    Please indicate the absolute path to the MAKEBLASTDB executable.

    Example: if your MAKEBLASTDB executable is found at absolute path /Users/myname/makeblastdb,
    type '/Users/myname/makeblastdb' and press Enter."""+'\n')
    makeblastdb_path = input(r"""    -----> Path to MAKEBLASTDB executable:  """)
    # 6-Specify prefix to files in database
    print(r"""
    ---------------------------------------------
    Prefix for files in BLASTN ALIGNMENT DATABASE
    ---------------------------------------------

    Because this script uses BLASTN (NCBI) and an alignment reference database, a common prefix identifier for the six
    database files generated by MAKEBLASTDB is needed.

    Please indicate a prefix to assign to each of the database files.

    Example: if your alignment reference was generated by MAKEBLASTDB from a fasta file called GRCh38.fa,
    the alignment database files will have been assigned the prefix 'GRCh38'; you would type 'GRCh38'
    and press Enter."""+'\n')
    db_prefix = input(r"""    -----> Prefix for alignment reference sequence database files:  """)
    # 7-Specify path to FIMO installation
    print(r"""
    ---------------------------
    Location of FIMO EXECUTABLE
    ---------------------------

    This script uses FIMO from the MEME suite of sequence analysis tools as its basis for determining distinct vs.
    common TFBSs.
    Please indicate the absolute path to the FIMO installation.

    Example: if your FIMO executable is found at absolute path /Users/myname/fimo, type '/Users/myname/fimo'
    and press Enter."""+'\n')
    fimo_path = input(r"""    -----> Path to FIMO executable:  """)
    # 8-Specify path to FIMO motif file.
    print(r"""
    ----------------------------
    Location of FIMO MOTIFS FILE
    ----------------------------

    This script uses FIMO from the meme suite of sequence analysis tools as its basis for determining distinct vs.
    common TFBS motifs.
    Please indicate the absolute path to the FIMO motifs file (containing position frequency matrix/matrices).

    When you're done entering the location of the motifs file, press Enter."""+'\n')
    fimo_motifs_path = input(r"""    -----> Path to FIMO motifs file:  """)
    # 9-Specify path to FIMO fasta-get-markov installation.
    print(r"""
    ---------------------------------------------
    Location of FIMO FASTA-GET-MARKOV EXECUTABLE
    ---------------------------------------------

    This script uses FIMO from the MEME suite of sequence analysis tools as its basis for determining distinct vs.
    common TFBSs.
    Please indicate an absolute path to the location of the FASTA-GET-MARKOV executable.

    When you're done entering the location of the executable, press Enter."""+'\n')
    fasta_get_markov_path = input(r"""    -----> Path to FASTA-GET-MARKOV executable:  """)
    # 10-Specify path to markov background file.
    print(r"""
    ------------------------------------------------------------
    Location of FIMO FASTA-GET-MARKOV BACKGROUND REFERENCE FILE
    ------------------------------------------------------------

    This script uses FIMO from the MEME suite of sequence analysis tools as its basis for determining distinct vs.
    common TFBSs.
    Please indicate an absolute path to the location of the fasta file you will use as your background reference
    (on which FASTA-GET-MARKOV will operate to generate a markov background file).
    
    When you're done entering the location of the reference sequence, press Enter."""+'\n')
    markov_background_file = input(r"""    -----> Path to background reference file:  """)
    # 11-Specify transcription factor (TF) of interest, for which to search for lost TFBS occurrences in alleles.
    print(r"""
    ------------------------------------------------
    TRANSCRIPTION FACTOR (TF) of interest (optional)
    ------------------------------------------------

    This script collates lost and gained TFBS for sample-associated allele(s) relative to a reference sequence;
    if detailed analysis of alleles that have lost TFBS matches for a specific transcription factor (TF) are
    desired,the identity of an individual TF of interest can be provided (optional).

    If you would like the script to further analyze alleles for TFBS matches to a specific TF, please indicate the
    TF here.  Otherwise, press Enter. 
    
    Important: Use only the standardized Entrez gene name for the TF of interest (such as NR3C1), rather than the
    matrix model stable ID (for example, MA0113 for NR3C1) or stable ID with version number (for example, MA0113.3
    for NR3C1).

    Example: if you are interested in losses of TFBS for the TF NR3C1, you would type 'NR3C1'
    and press Enter."""+'\n')
    TF_of_interest = input(r"""    -----> Transcription Factor (TF) of interest:  """)

    
# Define 'convert_bytes' and 'path_size' functions to be used in data collection for script_metrics.txt        
def convert_bytes(num):
    """
    This function converts bytes to convenient order of magnitude prefixes
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0

def path_size(given_path):
    """
    This function returns file or directory size
    """
    if os.path.isfile(given_path):
        file_info = os.stat(given_path)
        return convert_bytes(file_info.st_size)
    elif os.path.isdir(given_path):
        dir_info = os.stat(given_path)
        return convert_bytes(dir_info.st_size)
    
# Define 'merge' function to merge R1 & R2 reads
def merge(s1, s2):
    i = 0
    while not s2.startswith(s1[i:]):
        i += 1
    if i < len(s2):
        return s1[:i] + s2
    else:
        return 'no overlap'
      
# Define 'merge1' function to append two strings that do not overlap
def merge1(s1, s2):
    i = 0
    while not s2.startswith(s1[i:]):
        i += 1
    return s1[:i] + s2

# Define nt complement dictionary      
nt_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}
            
# Welcome/orient to script
print("""
    ==============================================
    CollatedMotifs.py v1.1
    ==============================================
    https://github.com/YamamotoLabUCSF/CollatedMotifs
    Committed 7-01-2021
    ----------------------------------------------
    This script accepts text to standard input, and returns a comparison of matches to
    transcription factor binding site (TFBS) motifs between a reference sequence and sample sequence(s).
    
    Requirements for operation:
    Python3
    BLASTN & MAKEBLASTDB from NCBI
    FIMO & FASTA-GET-MARKOV from the MEME suite of tools
    
    For usage details, please refer to README file at GitHub location and to the following manuscript:
        Ehmsen, Knuesel, Martinez, Asahina, Aridomi, Yamamoto (2021)
    
    Please cite usage as:
        CollatedMotifs.py
        Ehmsen, Knuesel, Martinez, Aridomi, Yamamoto (2021)
    
    ===========================================================================
    Welcome.  You will be prompted for the following user-specific information:
           * absolute paths to the following...
                * directories:
                   * output directory for output files
                   * directory containing fastq files for processing
                * files:
                   * single fasta file, containing reference sequence(s)
                     for processing by MAKEBLASTDB and FIMO
                   * single text file, to be used as FIMO motif file
                   * single text or fasta file, to be used by FASTA-GET-MARKOV
                     as basis for markov background file
                * executables:
                   * blastn
                   * makeblastdb
                   * fimo
                   * fasta-get-markov
            * single string as prefix to assign to 6 files of alignment database
            * single string indicating transcription factor (TF) of interest (optional)
              (use TF name as present in positional weight matrix file)

           You can expect the following directories and files in an output directory you specify:
           
                          `-----/alignment_database
                                        `----------*.nin
                                        `----------*.nhr
                                        `----------*.nog
                                        `----------*.nsd
                                        `----------*.nsg
                                        `----------*.nsi
                          `-----blastn_alignments.txt
                          `-----collated_TFBS.txt
                          `-----collated_TFBS.xlsx
                          `-----fasta.fa
                          `-----/fimo_out
                                        `----------cisml.xml
                                        `----------fimo.gff
                                        `----------fimo.html
                                        `----------fimo.tsv
                                        `----------fimo.xml
                          `-----/fimo_out_ref
                                        `----------cisml.xml
                                        `----------fimo.gff
                                        `----------fimo.html
                                        `----------fimo.tsv
                                        `----------fimo.xml
                          `-----markov_background.txt
                          `-----script_metrics.txt
        
    """)
    
input("Press Enter to continue...")

# Specify whether user input is provided at individual coached prompts or as single-list entry
print(r"""
    ---------------------------------------------------------------------
    User-specified input: choice of coached prompts vs. single list entry
    ---------------------------------------------------------------------

    Values for the user-specified input indicated above can be entered at individually coached command-line prompts
    (default), or as a single list of variables provided in a single command-line entry without coached prompts.

    To proceed with input at individual command-line PROMPTS, type 'Prompt' and press Enter;
    To proceed with input provided as a single LIST in one command-line entry, type 'List' and press Enter:  
    """)
user_input = input(r"""    -----> List or Prompt: """)

if user_input == 'Prompt':
    prompts()
elif user_input == 'List':
    print("""
    ----------------------------------
    User-specified input (list format)
    ----------------------------------
    
    Please paste a single list of input values directly at the command line prompt, specifying the following 10 (or 11) values.
    Press 'Enter' twice to complete.
    
    1-Location of OUTPUT DIRECTORY for output files
    2-Location of INPUT FILES (directory containing fastq files)
    3-Location of REFERENCE FASTA FILE
    4-Location of BLASTN EXECUTABLE
    5-Location of MAKEBLASTDB EXECUTABLE
    6-Prefix common to BLASTN sequence database files
    7-Location of FIMO EXECUTABLE
    8-Location of POSITION FREQUENCY MATRIX FILE
    9-Location of FASTA-GET-MARKOV EXECUTABLE
    10-Location of MARKOV BACKGROUND FILE
    11-Identity of TRANSCRIPTION FACTOR (TF) of interest (optional)
    
    """)
    input_list = []
    stopword = ""
    while True:
        input_str = input()
        if input_str.strip() == stopword:
            break
        else:
            input_list.append(input_str)
    output_directory = input_list[0].strip()
    fastq_directory = input_list[1].strip()
    fasta_ref = input_list[2].strip()
    blastn_path = input_list[3].strip()
    makeblastdb_path = input_list[4].strip()
    db_prefix = input_list[5].strip()
    fimo_path = input_list[6].strip()
    fimo_motifs_path = input_list[7].strip()
    fasta_get_markov_path = input_list[8].strip()
    markov_background_file = input_list[9].strip()
    if len(input_list) == 11:
        TF_of_interest = input_list[10].strip()
    else:
        TF_of_interest = ''
        
# Wait to create the directories and files until after input has been reviewed and accepted.
# Convert fastq_directory input to operating system-appropriate filepath.
output_directory = Path(str(output_directory))
# Convert fastq_directory input to operating system-appropriate filepath.
fastq_directory = Path(str(fastq_directory))
# Convert fasta_ref input to operating system-appropriate filepath.
fasta_ref = Path(str(fasta_ref))
# Convert blastn_path input to operating system-appropriate filepath.
blastn_path = Path(str(blastn_path))
# Convert makeblastdb_path input to operating system-appropriate filepath.
makeblastdb_path = Path(str(makeblastdb_path))
# Convert fimo_path input to operating system-appropriate filepath.
fimo_path = Path(str(fimo_path))
# Convert fimo_motifs_path input to operating system-appropriate filepath.
fimo_motifs_path = Path(str(fimo_motifs_path))
# Convert fasta_get_markov_path input to operating system-appropriate filepath.
fasta_get_markov_path = Path(str(fasta_get_markov_path))
# Convert markov_background_file input to operating system-appropriate filepath.
markov_background_file = Path(str(markov_background_file))

# Collect fastq files from directory
myFastqFilenames = [file for file in glob.glob(str(fastq_directory)+'/*') if Path(file).suffix in [".gz",".fastq"]]

#Sort fastq file names
myFastqFilenames = sorted(myFastqFilenames)
            
# Double-check whether entries look good:
print("""
---------------------------------------------------------------
Preparation for output:
Please double-check that your inputs were recorded as expected.
---------------------------------------------------------------""")

print("""
Your OUTPUT DIRECTORY was recorded as:
""")
print(str(output_directory))

print("""
Your directory containing fastq INPUT FILES was recorded as:
""")
print(str(fastq_directory))

print("""
Please wait a moment while the following file properties are retrieved and/or
calculated across the fastq files to be processed:

* Illumina sequencing run ID(s)
* Total number of fastq files
* Total number of sequencing reads
* Size distribution of fastq files

""")

# Collect Illumina run IDs from fastq files, consolidate to unique run IDs
runIDlist = []
for sourcefile in myFastqFilenames:
    if Path(sourcefile).suffix == ".gz":
        with gzip.open(sourcefile, "rt") as f:
            runID = ":".join(f.readline().split(":",-2)[:2])
            if not runID in runIDlist:
                runIDlist.append(runID) 
    elif Path(sourcefile).suffix == ".fastq":
        with open(sourcefile, "r") as f:
            runID = ":".join(f.readline().split(":",-2)[:2])
            if not runID in runIDlist:
                runIDlist.append(runID)

# Collect total read counts for fastq files
readcount = []
for sourcefile in myFastqFilenames:
    if Path(sourcefile).suffix == ".gz":
        with gzip.open(sourcefile, "rt") as f:    
            readcount.append(int(len((f).readlines())/4))
    elif Path(sourcefile).suffix == ".fastq":
        with open(sourcefile, "r") as f:
            readcount.append(int(len((f).readlines())/4))
        
# Collect file sizes for fastq files
filesize = []
for sourcefile in myFastqFilenames:
    if Path(sourcefile).suffix == ".gz":
        with gzip.open(sourcefile, "rt") as f:
            filesize.append(round((os.path.getsize(sourcefile)/1048576),5))
    elif Path(sourcefile).suffix == ".fastq":
        filesize.append(round((os.path.getsize(sourcefile)/1048576),5))

# fastq_overview prepares summation of fastq file names, their sizes, and read counts, to be reported in script_metrics.txt    
fastq_overview = list(zip(myFastqFilenames, filesize, readcount))

print("""The following data were collected:  """)
print("    Illumina sequencing run ID(s): ")
for i in runIDlist:
    print('        '+i)

print("    # of fastq files to process: {0}".format(len(myFastqFilenames)))

print("    size distribution of fastq files to process: \n      total... "+str(round((sum(file for file in filesize))))+' MB \n      range... max: '+str(round((max(file for file in filesize)),2))+' MB; min: '+str(round((min(file for file in filesize)),5))+' MB; median: '+str(round((numpy.median([file for file in filesize])),3))+' MB; mean +/- stdev: '+str(round((numpy.mean([file for file in filesize])),3))+' +/- '+str(round((numpy.std([file for file in filesize])),3))+' MB')

print("    read distribution within fastq files to process: \n      total... "+locale.format_string("%d", sum(readcount), grouping=True)+' reads \n      range... max: '+str((max(file for file in readcount)))+' reads; min: '+str((min(file for file in readcount)))+' reads; median: '+str((numpy.median([file for file in readcount])))+' reads; mean +/- stdev: '+str(round((numpy.mean([file for file in readcount]))))+' +/- '+str(round((numpy.std([file for file in readcount]))))+' reads')

# Activate this code block if wish to print total read # for each individual fastq file
#print("fastq files to process for genotyping:")
#for sourcefile in myFastqFilenames:
#    with open(sourcefile, "r") as f:
#        reads = int(len((f).readlines())/4)
#    print('    '+sourcefile, reads)

print("""
Your FASTA REFERENCE FILE location was recorded as:
""")
print(str(fasta_ref))

print("""
Your BLASTN EXECUTABLE location was recorded as:
""")
print(str(blastn_path))

print("""
Your MAKEBLASTDB EXECUTABLE location was recorded as:
""")
print(makeblastdb_path)

print("""
Your BLASTN DATABASE FILE PREFIX was recorded as:
""")
print(db_prefix)

print("""
Your FIMO EXECUTABLE location was recorded as:
""")
print(fimo_path)

print("""
Your POSITION FREQUENCY FILE location was recorded as:
""")
print(fimo_motifs_path)

# Examine the reference file and indicate the ID, number of motifs, etc. print out list of factors for query

motifcountlist = []
with open(fimo_motifs_path, 'r') as file:
    for line in file:
        if bool(re.search('MOTIF', line)): 
            motifcountlist.append(line.strip())

print("""
# of TFBS motifs examined: """+str(len(motifcountlist)))

motifID = [i.split(' ')[2] for i in motifcountlist]
motifID = sorted(motifID)
chunked_motifID = [motifID[i: i+10] for i in range(0, len(motifID), 10)]

print('Identities of TFBS motifs examined: ')

for row in chunked_motifID:
    itemnumber = (len(row)*'{: ^13} ').rstrip()
    print(itemnumber.format(*row))

print("""
Your FASTA-GET-MARKOV EXECUTABLE location was recorded as:
""")
print(fasta_get_markov_path)

print("""
Your MARKOV BACKGROUND FILE location was recorded as:
""")
print(markov_background_file)

print("""
Your TF of interest was recorded as:
""")
if TF_of_interest != '':
    print(TF_of_interest)
else:
    print('No TF of interest was provided')

check = input("""
Is this list accurately recorded? Type 'Y' or 'N': 
""")

if check == 'Y':
    pass
elif check == 'N':
    checkup = input("""
If you have corrections to make, please quit the active script and start again.
To continue in the script, type 'Continue' and press enter.
To quit the script, type 'Exit' and press enter, or press 'Ctrl+C'.  """)
    if checkup == 'Exit':
        exit(0)
    elif checkup == 'Continue':
        pass

# Start the clock on script operation duration
startTime = datetime.now()
startTimestr = str(startTime).split(' ')[1].split('.')[0]
          
# Proceed to file processing
# Generate the directory and its files (to accept content later in script)
path = str(output_directory)

if not os.path.exists(path):
    os.makedirs(path)

output_path = Path(output_directory)
    
# Create output files
filename_list = ['fasta.fa', 'blastn_alignments.txt', 'collated_TFBS.txt', 'markov_background.txt', 'script_metrics.txt']
      
# Define current date as prefix to all filenames
processdate = datetime.today().strftime("%m%d%Y")
      
for filename in filename_list:
    with open(os.path.join(path, processdate+'_'+filename), 'wb') as file:
        pass
    
# Collect RAM info
mem = virtual_memory()
ramem = mem.total/1073741824
   
# The file script_metrics.txt records script operation metadata (summarizes script input and performance)
# Use print redirection to write to target file, in append mode (begin script_metrics.txt)
with open(fasta_ref, "r") as f:
    ref_seqs = f.readlines()

filename = Path(str(output_path)+'/'+processdate+'_script_metrics.txt')
with open(filename, 'a') as f:
    print("""CollatedMotifs.py: Script Metrics
Date: """ + (datetime.today().strftime("%m/%d/%Y")) +
"""\n\nOperating system information:
    name: """ + socket.gethostname() +
'\n    platform: ' + platform.platform() +
'\n    RAM (GB): ' + str(ramem) +
'\n    physical CPU/effective CPU: ' + str(psutil.cpu_count(logical=False)) +'/'+ str(psutil.cpu_count()) +
'\n    executable: ' + psutil.Process().exe() +
"""\n\nUser-entered variables:
    output_directory: """+ str(output_directory) +
"\n    fastq_directory: "+ str(fastq_directory) +
"\n    fasta_ref: "+ str(fasta_ref) +
"\n    blastn_path: "+ str(blastn_path) +
"\n    makeblastdb_path: "+ str(makeblastdb_path) +
"\n    db_prefix: "+ str(db_prefix) +
"\n    fimo_path: "+ str(fimo_path) +
"\n    fimo_motifs_path: "+ str(fimo_motifs_path) +
"\n    fasta_get_markov_path: "+ str(fasta_get_markov_path) +
"\n    markov_background_file: "+ str(markov_background_file), file = f)
    if TF_of_interest == '':
        print("    TF_of_interest: none specified", file = f)
    else:
        print("    TF_of_interest: "+ TF_of_interest, file = f)    
    print("""\nfastq file information:
    Illumina sequencing run ID(s): """+ str(runIDlist).strip('[]').replace("'","") +
"\n    Number of fastq files processed: "+ str(len(myFastqFilenames)) +
"""\n    Size distribution of fastq files processed: 
        total... """ +str(round((sum(file for file in filesize))))+' MB \n        range... max: '+str(round((max(file for file in filesize)),2))+' MB; min: '+str(round((min(file for file in filesize)),5))+' MB; median: '+str(round((numpy.median([file for file in filesize])),3))+' MB; mean +/- stdev: '+str(round((numpy.mean([file for file in filesize])),3))+' +/- '+str(round((numpy.std([file for file in filesize])),3))+' MB' +
"\n    Read distribution within fastq files to process: \n        total... "+locale.format_string("%d", sum(readcount), grouping=True)+' reads \n        range... max: '+str((max(file for file in readcount)))+' reads; min: '+str((min(file for file in readcount)))+' reads; median: '+str((numpy.median([file for file in readcount])))+' reads; mean +/- stdev: '+str(round((numpy.mean([file for file in readcount]))))+' +/- '+str(round((numpy.std([file for file in readcount]))))+' reads', file = f)
    print("\nfastq files processed (name, size (MB), reads): ", file = f)
    for i in (sorted(fastq_overview)):
        print("    " + str(i).strip("()").replace("'",""), file = f)
    print("\nReference sequences provided in fasta_ref file: ", file = f)
    for i in ref_seqs:
        print("    " + i.strip('\n'), file = f)     
    print("\n# of TFBS motifs examined: "+str(len(motifcountlist))+
"\nIdentities of TFBS motifs examined: ", file = f)
    for row in chunked_motifID:
        itemnumber = (len(row)*'{: ^13} ').rstrip()
        print(itemnumber.format(*row), file = f)        
f.close()

print("""
Script is now using MAKEBLASTDB (NCBI) to prepare a reference sequence database for BLASTN alignments.
The output of this operation is a set of 6 database files in alignments_directory.""")

# Start the clock on makeblastdb and fastgetmarkov operations
startTime_makeblastdb_fastagetmarkov_operations = datetime.now()
      
# Construct alignment database, in alignment_database directory
# Reference sequence input
mydb_input = Path(fasta_ref)

# Alignment database directory
mydb_output = Path(str(output_directory)+'/alignment_database')

os.makedirs(mydb_output)

# 'Make blastn database' command (usage: makeblastdb -in mydb.fsa -parse_seqids -dbtype nucl -out path)
# 'Make blastn database' command (usage: makeblastdb -in mydb.fsa -parse_seqids -dbtype nucl -title string -out path)
cmd_makeblastndb = str(makeblastdb_path)+' -in '+str(mydb_input)+' -parse_seqids -dbtype nucl -out '+str(mydb_output)+'/'+db_prefix

os.system(cmd_makeblastndb)

print("""
Script is now using FASTA-GET-MARKOV (MEME) to prepare a background markov file for FIMO statistical operations.
The output of this operation is a single file, markov_background.txt, supplied to FIMO with sample and reference fasta
files.""")

# Construct background markov file to be used by FIMO
# Markov input file (markov background file)
markovbackground_input = Path(markov_background_file)

# Markov background output file  
markovbackground_output = Path(str(output_directory)+'/'+processdate+'_markov_background.txt')
      
# 'Make markov background file' command (usage: fasta-get-markov [options] [<sequence file> [<background file>]])
cmd_fastagetmarkov = str(fasta_get_markov_path)+' -dna '+str(markovbackground_input)+' '+str(markovbackground_output)

os.system(cmd_fastagetmarkov)

# Log makeblastdb and fastgetmarkov operations time duration
makeblastdb_fastagetmarkov_operationsDuration = str(datetime.now()- startTime_makeblastdb_fastagetmarkov_operations).split(':')[0]+' hr|'+str(datetime.now() - startTime_makeblastdb_fastagetmarkov_operations).split(':')[1]+' min|'+str(datetime.now() - startTime_makeblastdb_fastagetmarkov_operations).split(':')[2].split('.')[0]+' sec|'+str(datetime.now()- startTime_makeblastdb_fastagetmarkov_operations).split(':')[2].split('.')[1]+' microsec'

      
# Prepare fasta file for reads to be queried for matches to TFBS motifs
print("""
Script is now processing the top 5 unique read sequences for each fastq file (read identity and count [% of reads in sample]).
The output of this step is a fasta file (.fa) that will be created in the OUTPUT directory you indicated.

""")

# Start the clock on read count duration      
startTime_readcount = datetime.now()

# Populate fasta files for fasta.fa, in preparation for fimo analysis      
query_input = Path(str(output_directory)+'/'+processdate+'_fasta.fa')

# define Nextera adaptor sequence, in preparation to trim read 3' ends if necessary
adaptor_str = 'CTGTCTCTTATACACATCT'
adaptor_str_rev = 'AGATGTGTATAAGAGACAG'

# Merge R1 and R2 reads, if present, as single sequence
R1_file_list = [sourcefile for sourcefile in myFastqFilenames if bool(re.split('_',os.path.basename(sourcefile))[3] == 'R1')]
R2_file_list = [sourcefile for sourcefile in myFastqFilenames if bool(re.split('_',os.path.basename(sourcefile))[3] == 'R2')]  

# R1, R2 cluster mapping
processed_files_list = []
R1_R2_map_list = []
for sourcefile in myFastqFilenames:
    if sourcefile in processed_files_list:
        pass
    else:
        testname = ''.join(re.split('_',os.path.basename(sourcefile))[0:3])
        for sourcefile1 in R1_file_list:
            if testname == ''.join(re.split('_',os.path.basename(sourcefile1))[0:3]):
                R1 = sourcefile
                if sourcefile not in processed_files_list:
                    processed_files_list.append(sourcefile)
        for sourcefile2 in R2_file_list:
            if testname == ''.join(re.split('_',os.path.basename(sourcefile2))[0:3]):
                R2 = sourcefile2
                if sourcefile2 not in processed_files_list:
                    processed_files_list.append(sourcefile2)
        R1_R2_map_list.append((R1, R2))
        
# Make fasta file of read entries, direct top read count output and annotation to fasta.fa   
for file_pair in R1_R2_map_list:
    R1_file = file_pair[0]
    R2_file = file_pair[1]
    fastaname = re.split('_', os.path.basename(R1_file))
    cluster_sequence_R1_dict = {}
    cluster_sequence_R2_dict = {}
    cluster_sequence_R2_revcomp_dict = {}
    cluster_merged_R1_R2revcomp_dict = {}
    cluster_merged_R1_R2revcomp_dict2 = {}
    merged_read_list = []
    counter=()
    if Path(R1_file).suffix == ".gz":
        with gzip.open(R1_file, "rt") as f:
            lines_R1 = f.readlines()
    elif Path(R1_file).suffix == ".fastq":
        with open(R1_file, 'r') as f:
            lines_R1 = f.readlines()    
    for x in range(0,len(lines_R1),4):
        # trim adaptor sequence and up to 3' end of read from R1 sequence, if adaptor sequence found
        cluster_sequence_R1_dict[lines_R1[x].split(':')[5]+':'+lines_R1[x].split(':')[6].split(' ')[0]] = lines_R1[x+1].strip('\n')[:lines_R1[x+1].strip('\n').index(adaptor_str)] if adaptor_str in lines_R1[x+1].strip('\n') else lines_R1[x+1].strip('\n') 
    #cluster_IDs_list_R1 = [x.split(':')[5]+':'+x.split(':')[6].split(' ')[0] for x in lines_R1[0::4]]
    if Path(R2_file).suffix == ".gz":
        with gzip.open(R2_file, "rt") as f:
            lines_R2 = f.readlines()
    elif Path(R2_file).suffix == ".fastq":
        with open(R2_file, 'r') as f:
            lines_R2 = f.readlines()
    for x in range(0,len(lines_R2),4):
        # trim adaptor sequence and up to 3' end of read from R2 sequence, if adaptor sequence found
        cluster_sequence_R2_dict[lines_R2[x].split(':')[5]+':'+lines_R2[x].split(':')[6].split(' ')[0]] = lines_R2[x+1].strip('\n')[:lines_R2[x+1].strip('\n').index(adaptor_str)] if adaptor_str in lines_R2[x+1].strip('\n') else lines_R2[x+1].strip('\n') 
    #cluster_IDs_list_R2 = [x.split(':')[5]+':'+x.split(':')[6].split(' ')[0] for x in lines_R2[0::4]]
    for cluster in cluster_sequence_R2_dict:
        cluster_sequence_R2_revcomp_dict[cluster] = ''.join(reversed(''.join(nt_dict.get(nt) for nt in cluster_sequence_R2_dict.get(cluster))))
    for cluster in cluster_sequence_R1_dict:
        if cluster in cluster_sequence_R2_revcomp_dict:
            if merge(cluster_sequence_R1_dict.get(cluster), cluster_sequence_R2_revcomp_dict.get(cluster)) != 'no overlap':
                cluster_merged_R1_R2revcomp_dict[cluster] = merge(cluster_sequence_R1_dict.get(cluster), cluster_sequence_R2_revcomp_dict.get(cluster))
            else:
                cluster_merged_R1_R2revcomp_dict2[cluster] = merge1(cluster_sequence_R1_dict.get(cluster), cluster_sequence_R2_revcomp_dict.get(cluster))
    for cluster in cluster_merged_R1_R2revcomp_dict:
        merged_read_list.append(cluster_merged_R1_R2revcomp_dict.get(cluster))
    # create dictionary (counter) relating unique read sequence to its # of occurrences
    counter=Counter(merged_read_list)
    modified_read_list_top5 = []
    for index, i in enumerate(counter.most_common(5)):
        filtered1 = sum([x for x in counter.values() if x/(sum(counter.values())) > 0.01])
        filtered10 = sum([x for x in counter.values() if x/(sum(counter.values())) > 0.1])
        raw_freq = round((100*i[1]/sum(counter.values())),2)
        modified_read_list_top5.append([i[0], '['+str(i[1])+'/'+str(sum(counter.values()))+']', 'rank'+str(index+1), raw_freq, int(stats.percentileofscore([i for i in counter.values()], i[1], 'rank')), round((100*i[1]/sum([i[1] for i in counter.most_common(5)])),2), round((100*i[1]/filtered1),2) if filtered1 > 0 and raw_freq >= 1 else 'None', round((100*i[1]/filtered10),2) if filtered10 > 0 and raw_freq >= 10 else 'None'])
    with open(str(query_input), 'a+') as file:
        for i in modified_read_list_top5:
              file.write('>'+fastaname[0]+'_'+'R1+R2'+'_'+str(i[1])+'_'+i[2]+'_%totalreads:'+str(i[3])+'_percentile:'+str(i[4])+'_%top5reads:'+str(i[5])+'_%readsfilteredfor1%:'+str(i[6])+'_%readsfilteredfor10%:'+str(i[7])+'\n'+i[0]+'\n')
                
# Log read count time duration      
readcountDuration = str(datetime.now()- startTime_readcount).split(':')[0]+' hr|'+str(datetime.now() - startTime_readcount).split(':')[1]+' min|'+str(datetime.now() - startTime_readcount).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_readcount).split(':')[2].split('.')[1]+' microsec'

# Prepare alignments
print("""
Script is now aligning the top 5 unique read sequences for each fastq file to the reference sequence (BLASTN); the output of this step is a blastn text file that will be created in the OUTPUT directory you indicated.

""")

# Start the clock on blastn alignments duration
startTime_alignments = datetime.now()

# Process alignments relative to reference sequence database, using blastn
# Reference database
db_input = mydb_output / db_prefix

# Alignment output
query_output = str(output_directory)+'/'+processdate+'_blastn_alignments.txt'
                   
# Alignment command
cmd_align = str(blastn_path)+' -strand plus -query '+str(query_input)+' -db '+str(db_input)+' -out '+str(query_output)+' -gapopen 1 -gapextend 1 -outfmt "5"'

os.system(cmd_align)

# Log alignment time duration
alignmentsDuration = alignmentsDuration = str(datetime.now()- startTime_alignments).split(':')[0]+' hr|'+str(datetime.now()- startTime_alignments).split(':')[1]+' min|'+str(datetime.now()- startTime_alignments).split(':')[2].split('.')[0]+' sec|'+str(datetime.now()- startTime_alignments).split(':')[2].split('.')[1]+' microsec'

# Start the clock on allele definitions duration
startTime_alleles = datetime.now()

# Infer genotypes from blastn_alignments.txt file
# Import blastn alignments output as a list of strings (each string corresponds to a query alignment)
alignments_list = []
with open(str(query_output), 'r') as file:
    reader = file.read()
    for i,part in enumerate(reader.split('<Iteration_iter-num>')):
        alignments_list.append(part)
# Remove blastn header line from alignments_list
alignments_list = alignments_list[1:]

# Convert alignments_list to list of lists (i.e., each query alignment string is encapsulateed into its own sublist within alignments_list2)
alignments_list2 = [alignments_list[i:i+1] for i in range(0, len(alignments_list))]

# Identify & subset queries for which no alignments were found in reference database ('no hits found')
no_hits_list = []
for i in alignments_list2:
    if re.search('No hits found', str(i)):
        no_hits_list.append(str(i).split('<Iteration_query-def>')[1].split('</Iteration_query-def>')[0])

# Record sample names having reads with no alignment hits      
no_hits_samplename_list = []
for i in no_hits_list:
    samplename = i.split('_')[0]
    if samplename not in no_hits_samplename_list:
        no_hits_samplename_list.append(samplename)

# Within each sublist of alignments_list2, split each line into an individual string, remove beginning and trailing whitespace, and recapture specified subset of alignment information in alignments_list3
alignments_list3 = []
for i in alignments_list2:
    if str(i).split('<Iteration_query-def>')[1].split('</Iteration_query-def>')[0] not in no_hits_list:
        alignments_list3.append([y.strip() for x in i for y in x.split('\n') if y.strip().startswith(('<Iteration_query-ID>', '<Iteration_query-def>', '<Hit_num>', '<Hit_id>', '<Hit_def>', '<Hsp_hit-from>', '<Hsp_hit-to>', '<Hsp_qseq>', '<Hsp_hseq>', '<Hsp_midline>'))])
    
# Identify & subset reads with >1 alignment to sequences in reference database
# Some reads with >1 alignment will be recovered to 'reconstitute' hypothesized allele (if BLASTN has split the read into multiple 'hits' or 'high-scoring pairs' (hsp's) within the span of the user-provided reference sequence)

# There are in principle at least 3 ways a read could potentially align to >1 position in reference database (1 & 2a,b below):
# (1) same sequence span aligns to >1 different locus (disparate <Hit_id>'s)
# [unlike in Genotypes.py, this scenario is not anticipated in CollatedMotifs.py, unless significant sequence overlap occurs in user-provided fasta file containing reference sequence(s)]
# (2) one sequence span may be split into two (ore more) different alignment matches, because of intervening gap(s) or insertion(s) that exceed ~60 bp (an apparent BLASTN gap limit)
#    (a) if the two (or more) 'split matches' align to the same <Hit_id>, but to different coordinates of that <Hit_id>, they will be presented by BLASTN as belonging to the same <Hit_num>, but to different <Hsp_num> (Hsp=high scoring pair)
#    (b) if the two (or more) 'split matches' span different <Hit_id>'s (essentially different 'chunks' of sequence with unique names, as organized within the alignment database), they will be presented by BLASTN as belonging to different <Hit_num>
# These observations suggest that it is important to distinguish a read with alignment to >1 sequence as either one with poor resolution among >1 reference sequences (if >1 reference sequence is provided), vs. one that harbors sizeable deletions or insertions relative to the reference sequence
# CollatedMotifs.py assumes continuity of 2 or more hsp's if they are assigned to the same user-provided reference sequence,
# and therefore attempts to reconstitute hypothesized alleles that span multiple non-overlapping hsp's (but does not attempt to reconstitute across multiple hits or for ambiguous reconstructions from overlapping hsp's)

# Organize reads with multiple hit IDs
# These reads are deprecated (not further analyzed)
multiple_alignments_hits_list = []
for i in alignments_list3:
    if len(re.findall('<Hit_num>', str(i))) > 1:
        multiple_alignments_hits_list.append(i)

# Prepare dictionary linking sample names to their reads having >1 alignment to sequences in reference database    
multiple_alignments_samplename_list = []
for i in multiple_alignments_hits_list:
    multiple_alignments_samplename_list.append(i[1].split('>')[1].split('_')[0])
    
multiple_alignments_dict = {}
for i in multiple_alignments_samplename_list:
    multiple_alignments_dict ["{0}".format(i)] = tuple(x for x in multiple_alignments_hits_list if bool(re.search(i, x[1])))

# Organize reads with single hit, but multiple associated hsp IDs
# These reads will be processed separately to 'reconstitute' potential alleles, with high-scoring alignment pairs matched to the reference sequence, but split into separate matches due to intervening non-aligning span between the alignment matches
multiple_alignments_hsp_list = []
for i in alignments_list3:
    if len(re.findall('<Hit_num>', str(i))) > 1:
        pass
    elif len(re.findall('<Hsp_hit-from>', str(i))) > 1:
        multiple_alignments_hsp_list.append(i)
        
# filter for multiple hsp's that can be reasonably reconstructed (i.e., >1 hsp's that do not overlap. Overlapping hsp's cannot readily be reconstructed and alleles with overlapping hsp's will be deprecated from analysis)
alleles_with_multiple_hsps_that_can_be_reconstructed_list = []
for i in alignments_list3:
    count = 0
    for x in i:
        if re.search('<Hsp_hit-from>', x):
            count = count+1
    if count > 1:
        overlapping_hsp = False
        index_split_list = []
        hsp_list = []
        hsp_to_from_list = []
        for index, x in enumerate(i):
            if re.search('<Hsp_hit-from>', x):
                index_split_list.append(index)
        for index, y in enumerate(index_split_list):
            hsp_list.append(i[y:y+(index+1)*5])
        hsp_to_from_list_chunked = []
        for index, w in enumerate(range(0,len(hsp_list))):
            hsp_to_from_list_chunked.append(hsp_list[index])
        for index, hsp in enumerate(hsp_list):
            range_to_check = range(int(hsp[0].split('>')[1].split('<')[0]),int(hsp[1].split('>')[1].split('<')[0]))
            # ranges of other hsp's in hsp_list
            other_hsps_list = hsp_list.copy()
            del other_hsps_list[index]
            for other_hsp in other_hsps_list:
                other_hsp_range = range(int(other_hsp[0].split('>')[1].split('<')[0]),int(other_hsp[1].split('>')[1].split('<')[0]))
                if set(range_to_check).intersection(other_hsp_range):
                    overlapping_hsp = True
                else:
                    pass
        if overlapping_hsp is False:
            alleles_with_multiple_hsps_that_can_be_reconstructed_list.append(i)
            
# identify alleles with multiple hsp's that were not able to be reconstructed; note these later in script_metrics.txt
alleles_with_multiple_hsps_that_cannot_be_reconstructed_list = []
for i in multiple_alignments_hsp_list:
    if i not in alleles_with_multiple_hsps_that_can_be_reconstructed_list:
        alleles_with_multiple_hsps_that_cannot_be_reconstructed_list.append(i)

# attempt to reconstruct long-ranging alignment for alleles with >1 BLASTN hsp that can be reasonably reconstituted
reconstructed_alleles_with_multiple_hsps_list = []

with open(mydb_input) as file:
    ref_candidates = file.readlines()
ref_candidates2 = [i.strip('>\n') for i in ref_candidates]
iterator = iter(ref_candidates2)
ref_candidates3 = list(zip(iterator,iterator))

for allele in alleles_with_multiple_hsps_that_can_be_reconstructed_list:
    # get allele as it appears intact in fasta file
    with open(query_input) as file:
        for index, line in enumerate(file):
            if line.strip('\n') == '>'+allele[1].split('>')[1].split('<')[0]:
                allele_fasta = next(file).strip()
    allele_reconstruction_list = [i for i in allele[0:5]]
    allele_reconstruction_temp_list = []
    allele_ref = allele[3].split('>')[1].split('<')[0]
    for ref in ref_candidates3:
        if ref[0] == allele_ref:
            allele_ref_sequence = ref[1]
        else:
            pass
    # collect info re: hsp match spans and alignments to reference
    for hsp in range(0,int(((len(allele)-5)/5))):
        hsp_from = int(allele[5+hsp*5].split('>')[1].split('<')[0])
        hsp_to = int(allele[6+hsp*5].split('>')[1].split('<')[0])
        hsp_qseq = allele[7+hsp*5].split('>')[1].split('<')[0]
        hsp_hseq = allele[8+hsp*5].split('>')[1].split('<')[0]
        hsp_midline = allele[9+hsp*5].split('>')[1].split('<')[0]
        
        if hsp_from == 1:
            allele_reconstruction_temp_list.append((hsp_from, str(hsp_from)+':'+str(hsp_to), allele_ref_sequence[0:hsp_to], hsp_qseq, hsp_hseq, hsp_midline))
        else:
            allele_reconstruction_temp_list.append((hsp_from, str(hsp_from)+':'+str(hsp_to), allele_ref_sequence[hsp_from-1:hsp_to], hsp_qseq, hsp_hseq, hsp_midline))

    allele_span_list = []
    reference_span_list = []
    alignment_midline_list = []
    # prepare to account for bp spans from allele read as represented in fasta file, which are represented in hsps
    allele_fasta_span_bp_accounting_set = set()
    allele_fasta_span_bp = range(1,len(allele_fasta)+1)
    allele_fasta_spans_in_hsps = set()
    # prepare to account for bp spans from reference as represented in reference file, which are represented in hsps
    reference_span_bp_accounting_set = set()
    reference_span_bp = range(1,len(allele_ref_sequence)+1)
    reference_spans_in_hsps = set()
    # assess hsp match positions relative to ref span, and re-order if needed based on relative order of start and stop positions of hsp alignment
    for subregion in sorted(allele_reconstruction_temp_list):
        reference_spans_in_hsps.update(set(range(int(subregion[1].split(':')[0]),int(subregion[1].split(':')[1]))))
        match = re.search(subregion[3].replace('-',''), allele_fasta)
        allele_fasta_spans_in_hsps.update(range(match.span()[0],match.span()[1]))
        
    for index, subregion in enumerate(sorted(allele_reconstruction_temp_list)):
        # subregion[0] is start position of hsp alignment match in reference; subregion[1] is string form of hsp coordinate span (start:end) relative to reference
        # subregion[2] is direct sequence span from reference sequence; subregion[3] is query from allele; subregion[4] is hit from reference; subregion[5] is midline
        if index == 0:
            if subregion[0] == 1:
                allele_span_list.append(subregion[3])
                reference_span_list.append(subregion[4])
                alignment_midline_list.append(subregion[5])
                reference_span_bp_accounting_set.update(range(1,int(subregion[1].split(':')[1])))
                # check for coverage in allele as represented in fasta file
                match = re.search(subregion[3].replace('-',''), allele_fasta)
                allele_fasta_span_bp_accounting_set.update(range(match.span()[0],match.span()[1]))
                
            else:
                allele_span_list.append(subregion[3])
                reference_span_list.append(subregion[4])
                alignment_midline_list.append(subregion[5])
                reference_span_bp_accounting_set.update(range(int(subregion[1].split(':')[0]),int(subregion[1].split(':')[1])))
                # check for coverage in allele as represented in fasta file
                match = re.search(subregion[3].replace('-',''), allele_fasta)
                allele_fasta_span_bp_accounting_set.update(range(match.span()[0],match.span()[1]))
                #allele_span_list.append('-'*subregion[0])
                #reference_span_list.append(allele_ref_sequence[0:subregion[0]])
                #alignment_midline_list.append(' '*subregion[0])
                #reference_span_bp_accounting_set.update(range(1,int(subregion[1].split(':')[1]))) 
        elif len(sorted(allele_reconstruction_temp_list)) > index > 0:
            test_span = range(int(sorted(allele_reconstruction_temp_list)[index-1][1].split(':')[1]), int(subregion[1].split(':')[0]))
            
            if reference_span_bp_accounting_set.intersection(test_span):
                allele_span_list.append(subregion[3])
                reference_span_list.append(subregion[4])
                alignment_midline_list.append(subregion[5])
                reference_span_bp_accounting_set.update(range(int(subregion[1].split(':')[0]),int(subregion[1].split(':')[1])))
                match = re.search(subregion[3].replace('-',''), allele_fasta)
                allele_fasta_span_bp_accounting_set.update(range(match.span()[0],match.span()[1]))
            else:
                match = re.search(subregion[3].replace('-',''), allele_fasta)
                allele_fasta_span_bp_accounting_set.update(range(match.span()[0],match.span()[1]))
                bases_in_fasta_allele_not_accounted_for_in_alignment = sorted(set(range(sorted(allele_fasta_span_bp_accounting_set)[0],sorted(allele_fasta_span_bp_accounting_set)[-1]))-allele_fasta_span_bp_accounting_set)
                if len(bases_in_fasta_allele_not_accounted_for_in_alignment) > 0:
                    bases_to_add = allele_fasta[bases_in_fasta_allele_not_accounted_for_in_alignment[0]:bases_in_fasta_allele_not_accounted_for_in_alignment[-1]+1]
                    allele_fasta_span_bp_accounting_set.update(range(bases_in_fasta_allele_not_accounted_for_in_alignment[0],bases_in_fasta_allele_not_accounted_for_in_alignment[1]))
                else:
                    bases_to_add = ''
                allele_span_list.append(bases_to_add)
                allele_span_list.append('-'*(int(subregion[0])-int(sorted(allele_reconstruction_temp_list)[index-1][1].split(':')[1])-1-len(bases_to_add)))
                reference_span_list.append(allele_ref_sequence[int(sorted(allele_reconstruction_temp_list)[index-1][1].split(':')[1]):int(subregion[0])-1])
                alignment_midline_list.append(' '*(int(subregion[0])-int(sorted(allele_reconstruction_temp_list)[index-1][1].split(':')[1])-1))
                reference_span_bp_accounting_set.update(range(int(subregion[1].split(':')[0]),int(subregion[1].split(':')[1])))
                allele_span_list.append(subregion[3])
                reference_span_list.append(subregion[4])
                alignment_midline_list.append(subregion[5])
                reference_span_bp_accounting_set.update(range(int(sorted(allele_reconstruction_temp_list)[index-1][1].split(':')[1])-1,int(subregion[1].split(':')[0])))

    
    # missing region of allele sequence as it appears in fasta
    bases_in_fasta_allele_not_accounted_for_in_alignment = sorted(set(range(sorted(allele_fasta_span_bp_accounting_set)[0],sorted(allele_fasta_span_bp_accounting_set)[-1]))-allele_fasta_span_bp_accounting_set)
    if len(bases_in_fasta_allele_not_accounted_for_in_alignment) > 0:   
        bases_to_add = allele_fasta[bases_in_fasta_allele_not_accounted_for_in_alignment[0]:bases_in_fasta_allele_not_accounted_for_in_alignment[-1]+1]
    else:
        bases_to_add = ''
        
    reconstructed_hsp_from = str(int(sorted(list(reference_span_bp_accounting_set))[0])+1)
    reconstructed_hsp_to = str(int(sorted(list(reference_span_bp_accounting_set))[-1])+1)
    reconstructed_hsp_qseq =  ''.join(allele_span_list).strip('-')
    reconstructed_hsp_hseq = ''.join(reference_span_list).strip('-')
    reconstructed_hsp_midline = ''.join(alignment_midline_list)
    
    allele_reconstruction_list.append('<Hsp_hit-from>'+str(reconstructed_hsp_from)+'</Hsp_hit-from>')
    allele_reconstruction_list.append('<Hsp_hit-to>'+str(reconstructed_hsp_to)+'</Hsp_hit-to>')
    allele_reconstruction_list.append('<Hsp_qseq>'+reconstructed_hsp_qseq+'</Hsp_qseq>')
    allele_reconstruction_list.append('<Hsp_hseq>'+reconstructed_hsp_hseq+'</Hsp_hseq>')
    allele_reconstruction_list.append('<Hsp_midline>'+reconstructed_hsp_midline+'</Hsp_midline>')
    
    reconstructed_alleles_with_multiple_hsps_list.append(allele_reconstruction_list)  
    
    
# Prepare alignment_list4 for reads with exclusively 1 alignment hit in reference database, or hsp's that can be reasonably reconstructed
alignments_list4 = []
for i in alignments_list3:
    if i not in multiple_alignments_hsp_list:
        alignments_list4.append(i)
for i in reconstructed_alleles_with_multiple_hsps_list:
        alignments_list4.append(i)

# Use print redirection to write to target file, in append mode (append to script_metrics.txt)

filename = Path(str(output_path)+'/'+processdate+'_script_metrics.txt')
with open(filename, 'a') as f:
    print("\nRecord of ranked alleles deprecated from analysis output:", file = f)
    print("\n    No hits identified by BLASTN in alignment database: ", file = f)
    if len(no_hits_list) == 0:
        print("        None", file = f)    
    else:
        for i in no_hits_list:
            print("        "+i, file = f)
    print("\n    Multiple hits identified by BLASTN in alignment database: ", file = f)
    if len(multiple_alignments_hits_list) == 0:
        print("        None", file = f) 
    else:
        for i in multiple_alignments_hits_list:
            print("        "+i[1].split('>')[1].split('<')[0], file = f)  
    print("\n    >1 high-scoring pair (hsp) identified by BLASTN, but hsp's could not be reconstructed into a hypothesized allele: ", file = f)
    if len(alleles_with_multiple_hsps_that_cannot_be_reconstructed_list) == 0:
        print("        None", file = f) 
    else:
        for i in alleles_with_multiple_hsps_that_cannot_be_reconstructed_list:
            print("        "+i[1].split('>')[1].split('<')[0], file = f)
    print("\nRecord of ranked alleles reconstructed from >1 high-scoring pair (hsp):", file = f)
    print("\n    >1 high-scoring pair (hsp) identified by BLASTN, and hsp's were reconstructed into a hypothesized allele: ", file = f)
    if len(alleles_with_multiple_hsps_that_can_be_reconstructed_list) == 0:
        print("        None", file = f) 
    else:
        for i in alleles_with_multiple_hsps_that_can_be_reconstructed_list:
            print("        "+i[1].split('>')[1].split('<')[0], file = f)   
    print("\n", file = f)
f.close()

# Among lists containing alignment data in alignments_list4, determine which queries (reads) correspond to the same sample; where querydef = i[1].split(">")[1].split("_[")[0], reads belonging to the same sample share identical querydef
# Fasta deflines encode frequency metrics for reads, based on defline format:
# sampleID_[reads/total reads]_percentile_% read abundance_% top 10 reads_% reads filtered for 1%_% reads filtered for 10%
querydef_list = []
for i in alignments_list3:
    querydef = i[1].split(">")[1].split("_")[0]
    querydef_list.append(querydef)
    
querydef_uniq_list = []
for i in querydef_list:
    if i in querydef_uniq_list:
        pass
    else:
        querydef_uniq_list.append(i)

# Prepare dictionary relating sample IDs to their associated reads ('alleles')      
alignmentoutput_dict = {}
for i in querydef_uniq_list:
    alignmentoutput_dict["{0}".format(i)] = tuple(x for x in alignments_list4 if bool(re.search(i, x[1])))
    
# Identify & subset sample ID's that do not have output alleles (empty tuple values in dictionary)
empty_sampleIDs_list = []
for i in alignmentoutput_dict:
    if bool(alignmentoutput_dict.get(i) == ()):
        empty_sampleIDs_list.append(i)

# Make a copy of alignmentoutput_dict, removing dictionary keys with empty tuple values
alignmentoutput_dict2 = { k : v for k,v in alignmentoutput_dict.items() if v}
# Alignmentoutput_dict2 is the key dictionary for alignment information
      
# Log allele definitions time duration
allele_definitionsDuration = allele_definitionsDuration = str(datetime.now() - startTime_alleles).split(':')[0]+' hr|'+str(datetime.now() - startTime_alleles).split(':')[1]+' min|'+str(datetime.now() - startTime_alleles).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_alleles).split(':')[2].split('.')[1]+' microsec'

# FIMO operations: identify TFBS
print("""
Script is now using FIMO to identify TFBSs in reference(s) and allele(s).

""")

# Start the clock on FIMO operations duration
startTime_fimo = datetime.now()

# Identify matches to TFBS motifs in each reference sequence.
# Notes about FIMO: "FIMO takes as input one or more fixed-length motifs, represented as position-specific frequency matrices. These motifs can be generated from the MEME motif discovery algorithm, extracted from an existing motif database or created by hand using a simple text format. The program computes a log-likelihood ratio score (often referred to incorrectly as a ‘log-odds score’) for each motif with respect to each sequence position and converts these scores to P-values using dynamic programming (Staden, 1994), assuming a zero-order null model in which sequences are generated at random with user-specified per-letter background frequencies. Finally, FIMO employs a bootstrap method (Storey, 2002) to estimate false discovery rates (FDRs). Because the FDR is not monotonic relative to the P-value, FIMO instead reports for each P-value a corresponding q-value, which is defined as the minimal FDR threshold at which the P-value is deemed significant (Storey, 2003)".
# Usage: fimo [options] <motif file> <sequence file>
# Install MEME suite (updated July 31 2018), download Motif Databases (updated Dec 7 2017)
# http://meme-suite.org/doc/install.html?man_type=web#quick

# Reference sequence input
ref_input = Path(fasta_ref)

# Reference sequence(s): FIMO file output directory
ref_TFBS_output = output_path / 'fimo_out_ref'

# alleles: FIMO file output directory      
allele_TFBS_output = output_path / 'fimo_out'

# Reference sequence(s): FIMO command (usage: fimo --bfile <background file> <motif file> <sequence file>)       
cmd_TFBS = str(fimo_path)+' --bfile '+str(markovbackground_output)+' --o '+str(ref_TFBS_output)+' '+str(fimo_motifs_path)+' '+str(ref_input)

os.system(cmd_TFBS)

# Alleles: FIMO command (usage: fimo --bfile <background file> <motif file> <sequence file>) 
cmd_TFBS = str(fimo_path)+' --bfile '+str(markovbackground_output)+' --o '+str(allele_TFBS_output)+' '+str(fimo_motifs_path)+' '+str(query_input)

os.system(cmd_TFBS)

# Log FIMO operations time duration
fimoDuration = fimoDuration = str(datetime.now() - startTime_fimo).split(':')[0]+' hr|'+str(datetime.now() - startTime_fimo).split(':')[1]+' min|'+str(datetime.now() - startTime_fimo).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_fimo).split(':')[2].split('.')[1]+' microsec' 

# Start the clock on TFBS collation operations duration
startTime_TFBScollation = datetime.now()
      
# Now, have TFBS output files for both reference sequences and alleles, in two separate files, in two separate directories (fimo_out & fimo_out_ref)  
# Prepare dictionary of TFBSs ID'ed, for each reference sample 
with open(str(ref_TFBS_output)+'/fimo.tsv', 'r') as file:
    ref_lines = file.readlines()  

# Remove FIMO header lines, etc.
ref_lines = ref_lines[1:]
for line in ref_lines.copy():
    if len(line.split('\t')) < 10:
        ref_lines.remove(line)

# Convert to set for faster processing
ref_lines = set(ref_lines)

with open(str(ref_input)) as file:
    fasta_lines = file.readlines()
fasta_names = fasta_lines[0::2]
fasta_names = [i.strip('\n').strip('>') for i in fasta_names]
ref_set = set(fasta_names)

dict_ref_TFBS = {}
for ref in ref_set:
    dict_ref_TFBS[ref] = []

# Take 3rd field of all lines as search for presence of key in dictionary, and add line as string in value list of key (allele)
for line in ref_lines:
    if line.split('\t')[2].strip() in dict_ref_TFBS:
        dict_ref_TFBS[line.split('\t')[2]].append(line.strip())

# Prepare allele dictionary; first, populate with 'all_sites'.  Then, run comparison to 'dict_TFBS_ref' to define sites that are lost vs. gained relative to reference sequence.
# Prepare dictionary of TFBSs ID'ed, for each sample allele
    
dict_allele_TFBS = {}
for allele in alignmentoutput_dict2:
    dict_allele_TFBS[allele] = {}

for allele in alignmentoutput_dict2:
    for x in range(0, len(alignmentoutput_dict2.get(allele))):
        dict_allele_TFBS[allele].update({alignmentoutput_dict2.get(allele)[x][1].split(">")[1].split("<")[0]:[]})

with open(str(allele_TFBS_output)+'/fimo.tsv', 'r') as file:
    allele_lines = file.readlines()

# Remove FIMO header lines, etc.
allele_lines = allele_lines[1:]
for line in allele_lines.copy():
    if len(line.split('\t')) < 10:
        allele_lines.remove(line)

# Convert to set for faster processing
allele_lines = set(allele_lines)
      
# Populate each allele with its 'all_sites' information
for line in allele_lines:
    dict_allele_TFBS_sample_key = line.split('\t')[2].strip().split('_')[0]
    dict_allele_TFBS_allele_key = line.split('\t')[2].strip()
    if dict_allele_TFBS_sample_key in dict_allele_TFBS:
        if dict_allele_TFBS.get(dict_allele_TFBS_sample_key).get(dict_allele_TFBS_allele_key) is not None:
            dict_allele_TFBS.get(dict_allele_TFBS_sample_key).get(dict_allele_TFBS_allele_key).append(line.strip())
        
# Prepare synopsis dictionary with alleles as keys, and list of 3 sublists: gained, lost, all_sites      
# Run comparison of 'all_sites' information relative to reference allele information.
dict_allele_TFBS_synopsis = {}
for allele in alignmentoutput_dict2:
    dict_allele_TFBS_synopsis[allele] = {}
    
for allele in alignmentoutput_dict2:
    for x in range(0, len(alignmentoutput_dict2.get(allele))):
        dict_allele_TFBS_synopsis[allele].update({alignmentoutput_dict2.get(allele)[x][1].split(">")[1].split("<")[0]:{'gained':[],'lost':[],'all_sites':[], 'TFs':{}, 'allele_sequence':[
            alignmentoutput_dict2.get(allele)[x][7]]+[alignmentoutput_dict2.get(allele)[x][8]]+[alignmentoutput_dict2.get(allele)[x][9]]}})
        
for sample in dict_allele_TFBS_synopsis:
    for allele in dict_allele_TFBS_synopsis.get(sample):
        for motif in dict_allele_TFBS.get(sample).get(allele):
            dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites').append(motif.split('\t')[1]+' ('+motif.split('\t')[0]+'),'+motif.split('\t')[5]+','+motif.split('\t')[9]+','+motif.split('\t')[7]+','+motif.split('\t')[3]+','+motif.split('\t')[4])
            if motif.split('\t')[0]+' ('+motif.split('\t')[1]+')' not in dict_allele_TFBS_synopsis.get(sample).get(allele).get('TFs'):
                count = 1
                dict_allele_TFBS_synopsis.get(sample).get(allele).get('TFs').update({motif.split('\t')[0]+' ('+motif.split('\t')[1]+')':count})
            else:
                count = dict_allele_TFBS_synopsis.get(sample).get(allele).get('TFs').get(motif.split('\t')[0]+' ('+motif.split('\t')[1]+')')+1
                dict_allele_TFBS_synopsis.get(sample).get(allele).get('TFs').update({motif.split('\t')[0]+' ('+motif.split('\t')[1]+')':count})

# Run comparisons:
# Make ref_TFBS_synopsis dictionary
dict_ref_TFBS_synopsis = {}
for ref in dict_ref_TFBS:
    dict_ref_TFBS_synopsis[ref] = {'all_sites':[], 'TFs':{}}
    
# Summarize TF counts in reference sequences      
for ref in dict_ref_TFBS_synopsis:    
    for motif in dict_ref_TFBS.get(ref):
        if motif.split('\t')[0]+' ('+motif.split('\t')[1]+')' not in dict_ref_TFBS_synopsis.get(ref).get('TFs'):
            count = 1
            dict_ref_TFBS_synopsis.get(ref).get('TFs').update({motif.split('\t')[1]+' ('+motif.split('\t')[0]+')':count})
        else:
            count = dict_ref_TFBS_synopsis.get(ref).get('TFs').get(motif.split('\t')[0]+' ('+motif.split('\t')[1]+')')+1
            dict_ref_TFBS_synopsis.get(ref).get('TFs').update({motif.split('\t')[1]+' ('+motif.split('\t')[0]+')':count})
            
# Catalog TFBSs in reference sequences, in format akin to 'all_sites' format in dict_allele_TFBS_synopsis          
for ref in dict_ref_TFBS_synopsis:
    for motif in dict_ref_TFBS.get(ref):
        dict_ref_TFBS_synopsis.get(ref).get('all_sites').append(motif.split('\t')[1]+' ('+motif.split('\t')[0]+'),'+motif.split('\t')[5]+','+motif.split('\t')[9]+','+motif.split('\t')[7]+','+motif.split('\t')[3]+','+motif.split('\t')[4])
        
# Run comparisons, populating into dict_allele_TFBS_synopsis
ref_options = [ref for ref in dict_ref_TFBS]
for sample in dict_allele_TFBS_synopsis:
    # define reference sequence appropriate to sample
    for ref in ref_options:
        if re.search(ref, sample):
            sample_ref = ref
    for allele in dict_allele_TFBS_synopsis.get(sample):
        # check only for motifs that overlap aligned region between allele and ref
        for x in alignmentoutput_dict2.get(sample):
            if x[1].split('>')[1].split('<')[0] == allele:
                to_from_range = range(int(x[5].split('>')[1].split('<')[0]),int(x[6].split('>')[1].split('<')[0]))
        dict_allele_TFBS_synopsis.get(sample).get(allele)['ref_coordinates_span'] = to_from_range
        ref_spans_represented_in_allele_hsps_temp = []
        # get hsp aligment spans relative to reference coordinates
        for index, x in enumerate(alignments_list3):
            if x[1].split('>')[1].split('<')[0] == allele:
                for i in range(5, len(alignments_list3[index]), 5):
                    ref_spans_represented_in_allele_hsps_temp.append(range(int(alignments_list3[index][i].split('>')[1].split('<')[0])+1, int(alignments_list3[index][i+1].split('>')[1].split('<')[0])))
                    ref_spans_represented_in_allele_hsps = sorted(list(i) for i in ref_spans_represented_in_allele_hsps_temp)
        dict_allele_TFBS_synopsis.get(sample).get(allele)['hsp_alignment_spans'] = ref_spans_represented_in_allele_hsps   
        ref_TFBS_set_to_include_in_evaluation = []
        for ref_motif in dict_ref_TFBS_synopsis.get(sample_ref).get('all_sites'):
            ref_motif_range = range(int(ref_motif.split(',')[4]), int(ref_motif.split(',')[5]))
            if set(to_from_range).intersection(ref_motif_range):
                ref_TFBS_set_to_include_in_evaluation.append(ref_motif)  
        for motif in dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites'):
            # limit ref range
            if motif.split(',')[:4] in [i.split(',')[:4] for i in ref_TFBS_set_to_include_in_evaluation]:
                pass
            else:
                dict_allele_TFBS_synopsis.get(sample).get(allele).get('gained').append(motif)
        for motif in ref_TFBS_set_to_include_in_evaluation:
            if motif.split(',')[:4] in [i.split(',')[:4] for i in dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites')]:
                pass
            else:
                dict_allele_TFBS_synopsis.get(sample).get(allele).get('lost').append(motif)

# Add allele ranks to allele names
dict_allele_TFBS_synopsis_allele_ranks = {}

for sample in dict_allele_TFBS_synopsis:
    index_frequency_list = []
    for index, allele in enumerate(dict_allele_TFBS_synopsis.get(sample)):
        index_frequency_list.append((float(allele.split('_')[6].split(':')[1]), allele.split('_')[6], allele, index))
        index_frequency_list_sorted = sorted(index_frequency_list, reverse=True)
    dict_allele_TFBS_synopsis_allele_ranks[sample] = index_frequency_list_sorted
    
for sample in dict_allele_TFBS_synopsis:
    for index, allele in enumerate(dict_allele_TFBS_synopsis.get(sample)):
        for index, ranked_allele in enumerate(dict_allele_TFBS_synopsis_allele_ranks.get(sample)):
            if allele == ranked_allele[2]:
                allele_rank = index+1
        dict_allele_TFBS_synopsis.get(sample).get(allele).update({'allele_rank': allele_rank}) 

# TFBS interpretations
print("""
Script is now interpreting potential lost-regained TFBS pairs & lost-gained TFBS pairs.
If TF of interest was provided, script will also further assess alleles with loss of TFBS for TF of interest.

""")
        
# Take stock of gained and lost TFBS that positionally overlap in reference vs. allele.
# Begin interpretive assessment of whether TFBS 'gained' for given TFs may in fact be 'regained' TFBS, in alleles where TFBS for the same TF has been lost
# In other words, small local base changes that disrupt a TFBS for a given TF may nevertheless supply a distinct TFBS for the same TF,
# amounting to, in principle, a 'reconstitution' or preservation of potential TFBS for TF
# First, convert lost/gained TFBS for each allele (per sample) to dataframe, allele_TFBS_synopsis_df

sample_list = [] 
allele_count_list = []
allele_list = []
allele_sequence_list = []
reference_sequence_list = []
alignment_midline_list = []
TF_list = []
strand_list = []
gained_TFBS_sequence_list = []
lost_TFBS_sequence_list = []
p_val_list = []
lostvsgained_list = []
allele_start_coordinate_list = []
allele_stop_coordinate_list = []
ref_start_coordinate_list = []
ref_stop_coordinate_list = []

for sample in dict_allele_TFBS_synopsis:
    allele_count = 0
    for allele in dict_allele_TFBS_synopsis.get(sample):
        allele_count = allele_count+1
        for TF_class in dict_allele_TFBS_synopsis.get(sample).get(allele):
            if TF_class == 'lost':
                # if allele alignment is to a subset of a longer user-provided reference span, coordinates
                # must be converted to allele span coordinates (because coordinates in the 'lost' category are derived
                # from reference sequence coordinates in dict_allele_TFBS_synopsis
                if len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('lost')) == 0:
                    sample_list.append(sample)
                    allele_count_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_rank'))
                    allele_list.append(allele)
                    allele_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].strip())
                    reference_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[1].split('>')[1].split('<')[0].strip())
                    alignment_midline_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[2].split('>')[1].split('<')[0].strip())   
                    TF_list.append('n/a')
                    strand_list.append('n/a')
                    gained_TFBS_sequence_list.append('n/a')
                    lost_TFBS_sequence_list.append('n/a')
                    p_val_list.append('n/a')
                    allele_start_coordinate_list.append('n/a')
                    allele_stop_coordinate_list.append('n/a')
                    ref_start_coordinate_list.append('n/a')
                    ref_stop_coordinate_list.append('n/a')
                    lostvsgained_list.append("No TFBS lost")
                else:
                    # retrieve allele's alignment span to reference sequence
                    ref_span = [dict_allele_TFBS_synopsis.get(sample).get(allele).get('ref_coordinates_span')[0],dict_allele_TFBS_synopsis.get(sample).get(allele).get('ref_coordinates_span')[-1]]
                    # retrieve allele's intact span
                    allele_length = len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].replace('-',''))
                    # retrieve allele's alignment span(s) relative to reference sequence
                    dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')
                    # retrieve intervening span between allele alignment spans relative to reference sequence
                    if len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')) < 2:
                        span = int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].count('-'))
                    else:
                        span = int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[1][0])-int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[0][-1])  
                    for TF in dict_allele_TFBS_synopsis.get(sample).get(allele).get('lost'):
                        sample_list.append(sample)
                        allele_count_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_rank'))
                        allele_list.append(allele)
                        allele_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].strip())
                        reference_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[1].split('>')[1].split('<')[0].strip())
                        alignment_midline_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[2].split('>')[1].split('<')[0].strip())   
                        TF_list.append(TF.split(',')[0])
                        strand_list.append(TF.split(',')[1])
                        gained_TFBS_sequence_list.append('n/a')
                        lost_TFBS_sequence_list.append(TF.split(',')[2])                            
                        p_val_list.append(TF.split(',')[3])
                        allele_start_coordinate_list.append('n/a')
                        allele_stop_coordinate_list.append('n/a')
                        # make adjustments in recorded 'lost' TFBS coordinates, to reflect coordinates as defined in allele span rather than coordinates as defined in reference span
                        if set(range(int(TF.split(',')[4]), int(TF.split(',')[5]))).intersection(range(int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[0][0]),
                                                                                      int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[0][-1]+1))):
                            ref_start_coordinate_list.append(int(TF.split(',')[4])-ref_span[0])
                            ref_stop_coordinate_list.append(int(TF.split(',')[5])-ref_span[0])
                        elif set(range(int(TF.split(',')[4]), int(TF.split(',')[5]))).intersection(range(int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[0][1]+1),
                                                                                      int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[1][0]))):
                            ref_start_coordinate_list.append(int(TF.split(',')[4])-ref_span[0])
                            ref_stop_coordinate_list.append(int(TF.split(',')[5])-ref_span[0])
                        elif set(range(int(TF.split(',')[4]), int(TF.split(',')[5]))).intersection(range(int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[1][0]),
                                                                                      int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[1][-1]+1))):
                            ref_start_coordinate_list.append(int(TF.split(',')[4])-ref_span[0])
                            ref_stop_coordinate_list.append(int(TF.split(',')[5])-ref_span[0])                       
                        lostvsgained_list.append('lost')
            elif TF_class == 'gained':
                if len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('gained')) == 0:
                    sample_list.append(sample)
                    allele_count_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_rank'))
                    allele_list.append(allele)
                    allele_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].strip())
                    reference_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[1].split('>')[1].split('<')[0].strip())
                    alignment_midline_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[2].split('>')[1].split('<')[0].strip())
                    TF_list.append('n/a')
                    strand_list.append('n/a')
                    gained_TFBS_sequence_list.append('n/a')
                    lost_TFBS_sequence_list.append('n/a')
                    p_val_list.append('n/a')
                    allele_start_coordinate_list.append('n/a')
                    allele_stop_coordinate_list.append('n/a')
                    ref_start_coordinate_list.append('n/a')
                    ref_stop_coordinate_list.append('n/a')
                    lostvsgained_list.append("No TFBS gained")
                else:
                    # retrieve allele's alignment span to reference sequence
                    ref_span = [dict_allele_TFBS_synopsis.get(sample).get(allele).get('ref_coordinates_span')[0],dict_allele_TFBS_synopsis.get(sample).get(allele).get('ref_coordinates_span')[-1]]
                    # retrieve allele's intact span
                    allele_length = len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].replace('-',''))
                    # retrieve allele's alignment span(s) relative to reference sequence
                    dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')
                    # retrieve intervening span between allele alignment spans relative to reference sequence
                    if len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')) < 2:
                        span = int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].count('-'))
                    else:
                        span = int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[1][0])-int(dict_allele_TFBS_synopsis.get(sample).get(allele).get('hsp_alignment_spans')[0][-1])   
                    for TF in dict_allele_TFBS_synopsis.get(sample).get(allele).get('gained'):
                        sample_list.append(sample)
                        allele_count_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_rank'))
                        allele_list.append(allele)
                        allele_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].strip())
                        reference_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[1].split('>')[1].split('<')[0].strip())
                        alignment_midline_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[2].split('>')[1].split('<')[0].strip())
                        TF_list.append(TF.split(',')[0])
                        strand_list.append(TF.split(',')[1])
                        gained_TFBS_sequence_list.append(TF.split(',')[2])
                        lost_TFBS_sequence_list.append('n/a')                            
                        p_val_list.append(TF.split(',')[3])
                        allele_start_coordinate_list.append(TF.split(',')[4])
                        allele_stop_coordinate_list.append(TF.split(',')[5]) 
                        ref_start_coordinate_list.append('n/a')
                        ref_stop_coordinate_list.append('n/a')
                        lostvsgained_list.append('gained')           
            
# Prepare dataframe synopsis of samples and alleles, with individual rows detailing TFBS lost or gained 
allele_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_count_list, "allele ID":allele_list, "alignment query\n(allele sequence)":allele_sequence_list,
                                   "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list, "TF":TF_list,
                                   "strand":strand_list, 
                                   "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":lost_TFBS_sequence_list,
                                   "Lost TFBS start coordinate (in reference)":ref_start_coordinate_list,
                                   "Lost TFBS end coordinate (in reference)":ref_stop_coordinate_list,
                                   "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":gained_TFBS_sequence_list, 
                                   "Gained TFBS start coordinate (in allele)":allele_start_coordinate_list,
                                   "Gained TFBS end coordinate (in allele)":allele_stop_coordinate_list, 
                                   "p-val":p_val_list, "lost or gained in allele (relative to ref)?":lostvsgained_list}
allele_TFBS_synopsis_df = pd.DataFrame(allele_TFBS_synopsis_df_columns)

allele_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","lost or gained in allele (relative to ref)?"],ascending=[True, True, True, True, False])

# Add read counts and calculated frequencies to allele_TFBS_synopsis_df
read_count_list = [i.split('_')[2].strip('[]').split('/')[0] for i in allele_TFBS_synopsis_df['allele ID'].to_list()]
total_reads_list = [i.split('_')[2].strip('[]').split('/')[1] for i in allele_TFBS_synopsis_df['allele ID'].to_list()]
pct_total_reads_list = [i.split('_')[4].split(':')[1] for i in allele_TFBS_synopsis_df['allele ID'].to_list()]
pct_reads_filtered_for_1pct_list = [float(i.split('_')[7].split(':')[1]) if i.split('_')[7].split(':')[1] != 'None' else 0 for i in allele_TFBS_synopsis_df['allele ID'].to_list()]
pct_reads_filtered_for_10pct_list = [float(i.split('_')[8].split(':')[1]) if i.split('_')[8].split(':')[1] != 'None' else 0 for i in allele_TFBS_synopsis_df['allele ID'].to_list()]

# Add column with allele comment (comment if appropriate)
# Note, pre-processing reads with a read cleaning utility such as cutadapt, trimmomatic, or fastp may remove such reads/
# inferred alleles in advance, obviating need for this read length flag
allele_TFBS_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in allele_TFBS_synopsis_df['alignment query\n(allele sequence)'].to_list()]]

allele_TFBS_synopsis_df.insert(loc=3, column='reads', value=read_count_list)
allele_TFBS_synopsis_df.insert(loc=4, column='total reads', value=total_reads_list)
allele_TFBS_synopsis_df.insert(loc=5, column='% total reads', value=pct_total_reads_list)
allele_TFBS_synopsis_df.insert(loc=6, column='% reads filtered for reads <1%', value=pct_reads_filtered_for_1pct_list)
allele_TFBS_synopsis_df.insert(loc=7, column='% reads filtered for reads <10%', value=pct_reads_filtered_for_10pct_list)

# May-June 2021, revisited July
# assess whether allele_TFBS_synopsis_df lost/gained TFBS may be TFBS 'replacements'/cognates
lost_TFBS_list = []
gained_TFBS_list = []

allele_TFBS_synopsis_df_coordinates_updated = pd.DataFrame(columns=allele_TFBS_synopsis_df.columns)
span_between_aligning_blocks_allele_list = []

# coordinates for 'lost' TFBS have already been adjusted to reflect allele alignment span relative to user-provided reference sequence span
# coordinates for 'gained' TFBS (novel to an allele relative to reference sequence) need to be corrected for positions in allele that are beyond a deletion/insertion span
# (and would therefore not enable comparison to cognate coordinate position in reference sequence unless coordinates are adjusted for missing span)
for index, row in allele_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","lost or gained in allele (relative to ref)?"],ascending=[True, True, True, True, False]).iterrows():
    if row['lost or gained in allele (relative to ref)?'] == "lost":
        lost_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Lost TFBS start coordinate (in reference)'])+','+str(row['Lost TFBS end coordinate (in reference)'])+
                              ','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
        allele_TFBS_synopsis_df_coordinates_updated.loc[index] = row
        if re.search('-', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0]): 
            # position of longest non-corresponding span in allele, relative to reference (in alignment) (characteristic if deletion allele)
            span_between_aligning_blocks_allele_temp = re.search(max(re.findall(r'-+', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0])),dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0]).span()
            span_between_aligning_blocks_allele = tuple(value+1 for value in span_between_aligning_blocks_allele_temp)
            calculated_span_between_aligning_blocks_allele = span_between_aligning_blocks_allele[1]-span_between_aligning_blocks_allele[0]
            # position of longest non-corresponding span in reference, relative to allele (in alignment) (characteristic of insertion allele)
            # span_between_aligning_blocks_reference = re.search(max(re.findall(r'-+', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[1])),dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0]).span()
            span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
        else:
            span_between_aligning_blocks_allele_list.append('n/a') 
    elif row['lost or gained in allele (relative to ref)?'] == "gained":
        # retrieve allele's alignment span to reference sequence
        ref_span = [dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('ref_coordinates_span')[0],dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('ref_coordinates_span')[-1]]
        # retrieve allele's intact span
        allele_length = len(dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0].replace('-',''))
        # retrieve allele's alignment span(s) relative to reference sequence
        dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('hsp_alignment_spans')
        # retrieve intervening span between allele alignment spans relative to reference sequence
        # make adjustments in 'gained' TFBS coordinates, to reflect coordinates as defined in reference span rather than coordinates as defined in allele span
        # scenario where there was not >1 hsp detected by BLASTN (alignment is unsplit by BLASTN)
        if len(dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('hsp_alignment_spans')) < 2:
            # print(row['allele rank'], dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('hsp_alignment_spans'))
            # adjust allele coordinates relative to reference coordinates
            hsp_spans_relative_to_reference_seq = dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('hsp_alignment_spans')    
            basal_number = int(hsp_spans_relative_to_reference_seq[0][0])
            hsp_spans_relative_to_reference_seq_adjusted = []
            for x in hsp_spans_relative_to_reference_seq:
                hsp_spans_relative_to_reference_seq_adjusted.append([int(y)-basal_number for y in x])  
            if re.search('-', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0]): 
            # position of longest non-corresponding span in allele, relative to reference (in alignment) (characteristic if deletion allele)
                span_between_aligning_blocks_allele_temp = re.search(max(re.findall(r'-+', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0])),dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0]).span()
                span_between_aligning_blocks_allele = tuple(value+1 for value in span_between_aligning_blocks_allele_temp)
                calculated_span_between_aligning_blocks_allele = span_between_aligning_blocks_allele[1]-span_between_aligning_blocks_allele[0]
            # position of longest non-corresponding span in reference, relative to allele (in alignment) (characteristic of insertion allele)
                # span_between_aligning_blocks_reference = re.search(max(re.findall(r'-+', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[1])),dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0]).span()
                # scenario if no coordinate adjustment is needed (TFBS end coordinate occurs before largest alignment gap:
                if int(row['Gained TFBS end coordinate (in allele)']) in range(int(hsp_spans_relative_to_reference_seq_adjusted[0][0]),
                                                                                      int(span_between_aligning_blocks_allele[0])):
                    gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Gained TFBS start coordinate (in allele)'])+','+
                              str(row['Gained TFBS end coordinate (in allele)'])+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
                    allele_TFBS_synopsis_df_coordinates_updated.loc[index] = row
                    span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
                # scenario if coordinate adjustment is needed (TFBS end coordinate occurs between start of alignment gap and alignment end):
                elif int(row['Gained TFBS end coordinate (in allele)']) in range(int(span_between_aligning_blocks_allele[0]), int(hsp_spans_relative_to_reference_seq_adjusted[0][-1])):
                    gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(int(row['Gained TFBS start coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele)+','+
                              str(int(row['Gained TFBS end coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele)+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val']) 
                    # add row with updated coordinates
                    allele_TFBS_synopsis_df_coordinates_updated.loc[index] = {'sample':row['sample'],'allele rank':row['allele rank'], 
                                                                              'allele ID':row['allele ID'], 'reads':row['reads'], 'total reads':row['total reads'],
                                                                              '% total reads':row['% total reads'], '% reads filtered for reads <1%':row['% reads filtered for reads <1%'],
                                                                              '% reads filtered for reads <10%':row['% reads filtered for reads <10%'], 'alignment query\n(allele sequence)':row['alignment query\n(allele sequence)'],
                                                                              'alignment midline':row['alignment midline'], 'alignment hit\n(reference)':row['alignment hit\n(reference)'], 'TF':row['TF'],
                                                                              'strand':row['strand'], "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"],
                                                                              'Lost TFBS start coordinate (in reference)':row['Lost TFBS start coordinate (in reference)'], 'Lost TFBS end coordinate (in reference)':row['Lost TFBS end coordinate (in reference)'],
                                                                              "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"],
                                                                              'Gained TFBS start coordinate (in allele)':int(row['Gained TFBS start coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele,
                                                                              'Gained TFBS end coordinate (in allele)':int(row['Gained TFBS end coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele,
                                                                              'p-val':row['p-val'], 'lost or gained in allele (relative to ref)?':row['lost or gained in allele (relative to ref)?'], 'comment':row['comment']}
                    span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
            else:
                span_between_aligning_blocks_allele = 'n/a'
                gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(int(row['Gained TFBS start coordinate (in allele)']))+','+
                              str(int(row['Gained TFBS end coordinate (in allele)']))+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
                allele_TFBS_synopsis_df_coordinates_updated.loc[index] = row
                span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
        # scenario where there was >1 hsp detected by BLASTN (aligning segments were split by BLASTN and required reconstruction)
        else:
            span_between_aligning_blocks_allele_temp = re.search(max(re.findall(r'-+', dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0])),dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('allele_sequence')[0].split('>')[1].split('<')[0]).span()
            span_between_aligning_blocks_allele = tuple(value+1 for value in span_between_aligning_blocks_allele_temp)
            calculated_span_between_aligning_blocks_allele = span_between_aligning_blocks_allele[1]-span_between_aligning_blocks_allele[0]
            hsp_spans_relative_to_reference_seq = dict_allele_TFBS_synopsis.get(row['sample']).get(row['allele ID']).get('hsp_alignment_spans')    
            basal_number = int(hsp_spans_relative_to_reference_seq[0][0])
            hsp_spans_relative_to_reference_seq_adjusted = []
            for x in hsp_spans_relative_to_reference_seq:
                hsp_spans_relative_to_reference_seq_adjusted.append([int(y)-basal_number for y in x])  
            # condition for coordinates within first alignment block/hsp (no coordinate adjustment needed)
            if int(row['Gained TFBS end coordinate (in allele)']) in range(int(hsp_spans_relative_to_reference_seq_adjusted[0][0]),
                                                                                      int(span_between_aligning_blocks_allele[0])):
                gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Gained TFBS start coordinate (in allele)'])+','+
                              str(row['Gained TFBS end coordinate (in allele)'])+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
                allele_TFBS_synopsis_df_coordinates_updated.loc[index] = row
                span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
            # condition for coordinates within gap between alignments blocks/hsp's
            elif int(row['Gained TFBS end coordinate (in allele)']) in range(int(hsp_spans_relative_to_reference_seq_adjusted[0][-1])+1,
                                                                                      int(hsp_spans_relative_to_reference_seq_adjusted[1][0])):
                gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Gained TFBS start coordinate (in allele)'])+','+
                              str(row['Gained TFBS end coordinate (in allele)'])+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
                allele_TFBS_synopsis_df_coordinates_updated.loc[index] = row
                span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
            # condition for coordinates beyond gap between alignments blocks/hsp's (coordinate adjustment needed)
            elif int(row['Gained TFBS end coordinate (in allele)']) in range(int(hsp_spans_relative_to_reference_seq_adjusted[1][0]),
                                                                                      int(hsp_spans_relative_to_reference_seq_adjusted[1][-1])):
                gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(int(row['Gained TFBS start coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele)+','+
                              str(int(row['Gained TFBS end coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele)+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
                # add row with updated coordinates
                allele_TFBS_synopsis_df_coordinates_updated.loc[index] = {'sample':row['sample'],'allele rank':row['allele rank'], 
                                                                              'allele ID':row['allele ID'], 'reads':row['reads'], 'total reads':row['total reads'],
                                                                              '% total reads':row['% total reads'], '% reads filtered for reads <1%':row['% reads filtered for reads <1%'],
                                                                              '% reads filtered for reads <10%':row['% reads filtered for reads <10%'], 'alignment query\n(allele sequence)':row['alignment query\n(allele sequence)'],
                                                                              'alignment midline':row['alignment midline'], 'alignment hit\n(reference)':row['alignment hit\n(reference)'], 'TF':row['TF'],
                                                                              'strand':row['strand'], "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"],
                                                                              'Lost TFBS start coordinate (in reference)':row['Lost TFBS start coordinate (in reference)'], 'Lost TFBS end coordinate (in reference)':row['Lost TFBS end coordinate (in reference)'],
                                                                              "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"],
                                                                              'Gained TFBS start coordinate (in allele)':int(row['Gained TFBS start coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele,
                                                                              'Gained TFBS end coordinate (in allele)':int(row['Gained TFBS end coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele,
                                                                              'p-val':row['p-val'], 'lost or gained in allele (relative to ref)?':row['lost or gained in allele (relative to ref)?'], 'comment':row['comment']}
                span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)
            else:
                gained_TFBS_list.append(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(int(row['Gained TFBS start coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele)+','+
                              str(int(row['Gained TFBS end coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele)+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
                # add row with updated coordinates
                allele_TFBS_synopsis_df_coordinates_updated.loc[index] = {'sample':row['sample'],'allele rank':row['allele rank'], 
                                                                              'allele ID':row['allele ID'], 'reads':row['reads'], 'total reads':row['total reads'],
                                                                              '% total reads':row['% total reads'], '% reads filtered for reads <1%':row['% reads filtered for reads <1%'],
                                                                              '% reads filtered for reads <10%':row['% reads filtered for reads <10%'], 'alignment query\n(allele sequence)':row['alignment query\n(allele sequence)'],
                                                                              'alignment midline':row['alignment midline'], 'alignment hit\n(reference)':row['alignment hit\n(reference)'], 'TF':row['TF'],
                                                                              'strand':row['strand'], "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"],
                                                                              'Lost TFBS start coordinate (in reference)':row['Lost TFBS start coordinate (in reference)'], 'Lost TFBS end coordinate (in reference)':row['Lost TFBS end coordinate (in reference)'],
                                                                              "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"],
                                                                              'Gained TFBS start coordinate (in allele)':int(row['Gained TFBS start coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele,
                                                                              'Gained TFBS end coordinate (in allele)':int(row['Gained TFBS end coordinate (in allele)'])+calculated_span_between_aligning_blocks_allele,
                                                                              'p-val':row['p-val'], 'lost or gained in allele (relative to ref)?':row['lost or gained in allele (relative to ref)?'], 'comment':row['comment']}
                span_between_aligning_blocks_allele_list.append(span_between_aligning_blocks_allele)

allele_TFBS_synopsis_df_coordinates_updated['span between alignment blocks'] = span_between_aligning_blocks_allele_list
                

potential_matched_TFBS_pairs_list = []

for i in lost_TFBS_list:
    for x in gained_TFBS_list:
        if x.split(',')[4] == 'n/a' or x.split(',')[5] == 'n/a':
            pass
        else:
            if i.split(',')[:4] == x.split(',')[:4]:
                lost_range = range(int(i.split(',')[4]), int(i.split(',')[5])+1)
                gained_range = range(int(x.split(',')[4]), int(x.split(',')[5])+1)
                if len(set(lost_range) & set(gained_range)) > 0:
                    potential_matched_TFBS_pairs_list.append((i,x))
                
unpaired_TFBS_gains_list = list(set(gained_TFBS_list) - 
                                set([i[0] for i in potential_matched_TFBS_pairs_list]+[i[1] for i in potential_matched_TFBS_pairs_list]))
unpaired_TFBS_losses_list = list(set(lost_TFBS_list) - 
                                set([i[0] for i in potential_matched_TFBS_pairs_list]+[i[1] for i in potential_matched_TFBS_pairs_list]))

# July 2021
# Reconstitute data for categories ('no TFBS predicted as lost or gained in allele', 'predicted TFBS loss (TFBS lost in allele)',
# 'predicted TFBS gain (novel to allele)')
sample_list = []
allele_rank_list = []
allele_list = []
read_count_list = []
total_reads_count_list = []
pct_total_reads_list = []
pct_reads_filtered_for_1pct_list = []
pct_reads_filtered_for_10pct_list = []
allele_sequence_list = []
reference_sequence_list = []
alignment_midline_list = []
TF_list = []
strand_list = []
lost_TFBS_sequence_list = []
gained_TFBS_sequence_list = []
allele_start_coordinate_list = []
allele_stop_coordinate_list = []
ref_start_coordinate_list = []
ref_stop_coordinate_list = []
lost_TFBS_pval_list = []
gained_TFBS_pval_list = []
predicted_lost_gained_pair_list = []

for index, row in allele_TFBS_synopsis_df_coordinates_updated.iterrows():
    if row['lost or gained in allele (relative to ref)?'] == 'lost':
        lost_test_phrase = ''.join(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Lost TFBS start coordinate (in reference)'])+','+str(row['Lost TFBS end coordinate (in reference)'])+
                              ','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
        if lost_test_phrase in set([i[0] for i in potential_matched_TFBS_pairs_list]+[i[1] for i in potential_matched_TFBS_pairs_list]):
            for match_pair in potential_matched_TFBS_pairs_list:
                if lost_test_phrase in match_pair:
                    sample_list.append(row['sample'])
                    allele_rank_list.append(row['allele rank'])
                    allele_list.append(row['allele ID'])
                    read_count_list.append(row['reads'])
                    total_reads_count_list.append(row['total reads'])
                    pct_total_reads_list.append(row['% total reads'])
                    pct_reads_filtered_for_1pct_list.append(row['% reads filtered for reads <1%'])
                    pct_reads_filtered_for_10pct_list.append(row['% reads filtered for reads <10%'])
                    allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                    reference_sequence_list.append(row['alignment hit\n(reference)'])
                    alignment_midline_list.append(row['alignment midline'])
                    TF_list.append(row['TF'])
                    strand_list.append(row['strand'])
                    lost_TFBS_sequence_list.append(match_pair[0].split(',')[7])
                    ref_start_coordinate_list.append(match_pair[0].split(',')[4])
                    ref_stop_coordinate_list.append(match_pair[0].split(',')[5])
                    gained_TFBS_sequence_list.append(match_pair[1].split(',')[7])
                    allele_start_coordinate_list.append(match_pair[1].split(',')[4])
                    allele_stop_coordinate_list.append(match_pair[1].split(',')[5])
                    lost_TFBS_pval_list.append(match_pair[0].split(',')[8])
                    gained_TFBS_pval_list.append(match_pair[1].split(',')[8])
                    predicted_lost_gained_pair_list.append('predicted lost-regained TFBS pair')
                else:
                    pass
        elif lost_test_phrase in unpaired_TFBS_losses_list:
            sample_list.append(row['sample'])
            allele_rank_list.append(row['allele rank'])
            allele_list.append(row['allele ID'])
            read_count_list.append(row['reads'])
            total_reads_count_list.append(row['total reads'])
            pct_total_reads_list.append(row['% total reads'])
            pct_reads_filtered_for_1pct_list.append(row['% reads filtered for reads <1%'])
            pct_reads_filtered_for_10pct_list.append(row['% reads filtered for reads <10%'])
            allele_sequence_list.append(row['alignment query\n(allele sequence)'])
            reference_sequence_list.append(row['alignment hit\n(reference)'])
            alignment_midline_list.append(row['alignment midline'])
            TF_list.append(row['TF'])
            strand_list.append(row['strand'])
            lost_TFBS_sequence_list.append(lost_test_phrase.split(',')[7])
            ref_start_coordinate_list.append(lost_test_phrase.split(',')[4])
            ref_stop_coordinate_list.append(lost_test_phrase.split(',')[5])
            gained_TFBS_sequence_list.append('n/a')
            allele_start_coordinate_list.append('n/a')
            allele_stop_coordinate_list.append('n/a')
            lost_TFBS_pval_list.append(lost_test_phrase.split(',')[8])
            gained_TFBS_pval_list.append('n/a (>1e-4 threshold)')
            predicted_lost_gained_pair_list.append('predicted TFBS loss (TFBS lost in allele)')
    elif row['lost or gained in allele (relative to ref)?'] == 'gained':    
        gained_test_phrase = ''.join(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Gained TFBS start coordinate (in allele)'])+','+
                              str(row['Gained TFBS end coordinate (in allele)'])+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
        if gained_test_phrase in unpaired_TFBS_gains_list:
            sample_list.append(row['sample'])
            allele_rank_list.append(row['allele rank'])
            allele_list.append(row['allele ID'])
            read_count_list.append(row['reads'])
            total_reads_count_list.append(row['total reads'])
            pct_total_reads_list.append(row['% total reads'])
            pct_reads_filtered_for_1pct_list.append(row['% reads filtered for reads <1%'])
            pct_reads_filtered_for_10pct_list.append(row['% reads filtered for reads <10%'])
            allele_sequence_list.append(row['alignment query\n(allele sequence)'])
            reference_sequence_list.append(row['alignment hit\n(reference)'])
            alignment_midline_list.append(row['alignment midline'])
            TF_list.append(row['TF'])
            strand_list.append(row['strand'])
            lost_TFBS_sequence_list.append('n/a')
            ref_start_coordinate_list.append('n/a')
            ref_stop_coordinate_list.append('n/a')
            gained_TFBS_sequence_list.append(gained_test_phrase.split(',')[7])
            allele_start_coordinate_list.append(gained_test_phrase.split(',')[4])
            allele_stop_coordinate_list.append(gained_test_phrase.split(',')[5])
            lost_TFBS_pval_list.append('n/a (>1e-4 threshold)')
            gained_TFBS_pval_list.append(gained_test_phrase.split(',')[8])
            predicted_lost_gained_pair_list.append('predicted TFBS gain (novel to allele)')
    elif row['lost or gained in allele (relative to ref)?'] == 'No TFBS gained':
        lost_test_phrase = ''.join(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Lost TFBS start coordinate (in reference)'])+','+str(row['Lost TFBS end coordinate (in reference)'])+
                              ','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val'])
        if lost_test_phrase in unpaired_TFBS_losses_list:
            pass
        else:
            if row['allele ID'] in allele_list:
                pass
            else:
                sample_list.append(row['sample'])
                allele_rank_list.append(row['allele rank'])
                allele_list.append(row['allele ID'])
                read_count_list.append(row['reads'])
                total_reads_count_list.append(row['total reads'])
                pct_total_reads_list.append(row['% total reads'])
                pct_reads_filtered_for_1pct_list.append(row['% reads filtered for reads <1%'])
                pct_reads_filtered_for_10pct_list.append(row['% reads filtered for reads <10%'])
                allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                reference_sequence_list.append(row['alignment hit\n(reference)'])
                alignment_midline_list.append(row['alignment midline'])
                TF_list.append(row['TF'])
                strand_list.append(row['strand'])
                lost_TFBS_sequence_list.append('n/a')
                ref_start_coordinate_list.append('n/a')
                ref_stop_coordinate_list.append('n/a')
                gained_TFBS_sequence_list.append('n/a')
                allele_start_coordinate_list.append('n/a')
                allele_stop_coordinate_list.append('n/a')
                lost_TFBS_pval_list.append('n/a')
                gained_TFBS_pval_list.append('n/a')
                predicted_lost_gained_pair_list.append('no TFBS predicted as lost or gained in allele')
    elif row['lost or gained in allele (relative to ref)?'] == 'No TFBS lost':
        gained_test_phrase = ''.join(row['sample']+','+str(row['allele rank'])+','+row['TF']+','+row['strand']+','+
                              str(row['Gained TFBS start coordinate (in allele)'])+','+
                              str(row['Gained TFBS end coordinate (in allele)'])+','+row['lost or gained in allele (relative to ref)?']+','
                              +row["Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'"]+
                             ','+row['p-val']) 
        if gained_test_phrase in unpaired_TFBS_gains_list:
            if row['allele ID'] in allele_list:
                pass
            else:
                sample_list.append(row['sample'])
                allele_rank_list.append(row['allele rank'])
                allele_list.append(row['allele ID'])
                read_count_list.append(row['reads'])
                total_reads_count_list.append(row['total reads'])
                pct_total_reads_list.append(row['% total reads'])
                pct_reads_filtered_for_1pct_list.append(row['% reads filtered for reads <1%'])
                pct_reads_filtered_for_10pct_list.append(row['% reads filtered for reads <10%'])
                allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                reference_sequence_list.append(row['alignment hit\n(reference)'])
                alignment_midline_list.append(row['alignment midline'])
                TF_list.append(row['TF'])
                strand_list.append(row['strand'])
                lost_TFBS_sequence_list.append('n/a')
                ref_start_coordinate_list.append('n/a')
                ref_stop_coordinate_list.append('n/a')
                gained_TFBS_sequence_list.append('n/a')
                allele_start_coordinate_list.append('n/a')
                allele_stop_coordinate_list.append('n/a')
                lost_TFBS_pval_list.append('n/a')
                gained_TFBS_pval_list.append('n/a')
                predicted_lost_gained_pair_list.append('no TFBS predicted as lost or gained in allele')
                         
# Prepare dataframe synopsis of samples and alleles, with individual rows mapping potential lost-regained TFBS pairs for in-common TFs
interpreted_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list, 
                                        'read count':read_count_list,
                                        'total reads':total_reads_count_list,
                                        '% total reads':pct_total_reads_list,
                                        '% reads filtered for reads <1%':pct_reads_filtered_for_1pct_list,
                                        '% reads filtered for reads <10%':pct_reads_filtered_for_10pct_list,
                                        "alignment query\n(allele sequence)":allele_sequence_list,
                                        "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                        "TF":TF_list,
                                        "strand":strand_list, 
                                        "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":lost_TFBS_sequence_list,
                                        "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":gained_TFBS_sequence_list, 
                                        "Lost TFBS coordinate start (in reference)":ref_start_coordinate_list,
                                        "Lost TFBS coordinate end (in reference)":ref_stop_coordinate_list,
                                        "Gained TFBS coordinate start (in allele)":allele_start_coordinate_list,
                                        "Gained TFBS coordinate end (in allele)":allele_stop_coordinate_list,
                                        "Lost TFBS p-val (in reference)":lost_TFBS_pval_list,
                                        "Gained TBFS p-val (in allele)":gained_TFBS_pval_list,
                                        "interpretation":predicted_lost_gained_pair_list}

interpreted_TFBS_synopsis_df = pd.DataFrame(interpreted_TFBS_synopsis_df_columns)

# Add column with allele comment (comment if appropriate)
interpreted_TFBS_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in interpreted_TFBS_synopsis_df['alignment query\n(allele sequence)'].to_list()]]

interpreted_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","interpretation"],ascending=[True, True, True, True, False])

# clean up interpreted_TFBS_synopsis_df; some sample alleles have been assigned rows that include the
# label 'no TFBS predicted as lost or gained in allele', because of order of operations above,
# but in fact have TFBS(s) predicted as lost or gained; find and remove these rows from interpreted_TFBS_synopsis_df
# (in new dataframe called interpreted_TFBS_syopsis_df_updated)
allele_rank_count_list = []
sample_set = set(interpreted_TFBS_synopsis_df['sample'].to_list())

for sample in sample_set:
    allele_rank_list = []
    for index, row in interpreted_TFBS_synopsis_df.iterrows():
        if row['sample'] == sample:
            if row['allele rank'] not in allele_rank_list:
                allele_rank_list.append(row['allele rank'])
    allele_rank_count_list.append((sample, sorted(allele_rank_list)))

suspects = []
for sample in allele_rank_count_list:
    for allele_rank in sample[1]:
        for index1, row in interpreted_TFBS_synopsis_df.iterrows():
            if row['sample'] == sample[0] and row['allele rank'] == allele_rank and row['interpretation'] == 'no TFBS predicted as lost or gained in allele':
                test_singularity = True
                for index, row in interpreted_TFBS_synopsis_df.iterrows():
                    if row['sample'] == sample[0] and row['allele rank'] == allele_rank and row['interpretation'] != 'no TFBS predicted as lost or gained in allele':
                        test_singularity = False
                if test_singularity == True:
                    pass
                elif test_singularity == False:
                    suspects.append((sample[0], allele_rank, index1))
                    
index_drop_list = [i[2] for i in suspects]
interpreted_TFBS_synopsis_df_updated = interpreted_TFBS_synopsis_df.drop(index_drop_list)

# Now, query output for specific TF of interest, based on the following properties:
# (1) lost without corresponding re-gain, (2) lost without corresponding re-gain and without positionally coinciding TFBS for new TF
# TF of interest was provided by user at script outset (encoded by variable samples_list = set(interpreted_TFBS_synopsis_df['sample'].to_list())

samples_list = set(interpreted_TFBS_synopsis_df['sample'].to_list())
                   
if TF_of_interest == '':
    pass
else:
    interpretation_dict = {}
    for sample in samples_list:
        sample_interpretation_dict = {}
        for index, row in interpreted_TFBS_synopsis_df_updated.sort_values(by=['sample','allele rank','TF',"strand","interpretation"],ascending=[True, True, True, True, False]).iterrows():
            if row['sample'] == sample:
                if row['interpretation'] == "no TFBS predicted as lost or gained in allele":
                    sample_interpretation_dict[row['allele rank']] = [row.to_list()]  
                elif re.search(r'\b'+TF_of_interest+r'\b', row['TF']) and row['interpretation'] == 'predicted TFBS loss (TFBS lost in allele)':
                    if row['allele rank'] not in sample_interpretation_dict:
                        sample_interpretation_dict[row['allele rank']] = [row.to_list()]
                    elif row['allele rank'] in sample_interpretation_dict:
                        sample_interpretation_dict[row['allele rank']].append(row.to_list())
                elif re.search(r'\b'+TF_of_interest+r'\b', row['TF']) and row['interpretation'] == 'predicted TFBS gain (novel to allele)':
                    if row['allele rank'] not in sample_interpretation_dict:
                        sample_interpretation_dict[row['allele rank']] = [row.to_list()]
                    elif row['allele rank'] in sample_interpretation_dict:
                        sample_interpretation_dict[row['allele rank']].append(row.to_list())
                elif re.search(r'\b'+TF_of_interest+r'\b', row['TF']) and row['interpretation'] == 'predicted lost-regained TFBS pair':
                    if row['allele rank'] not in sample_interpretation_dict:
                        sample_interpretation_dict[row['allele rank']] = [row.to_list()]
                    elif row['allele rank'] in sample_interpretation_dict:
                        sample_interpretation_dict[row['allele rank']].append(row.to_list())         
        interpretation_dict[sample] = sample_interpretation_dict

# July 2021
# Now assess potential samples of particular interest based on 2 criteria above
# Iterate through alleles to bin alleles for samples among the indicated lists below
# (1) lost without corresponding gain, (2) lost without corresponding gain and without positionally coinciding TFBS for new TF
# for (2), check for 'predicted TFBS gain (novel to allele)' that coincides with position of TFBS loss for TF of interest

if TF_of_interest == '':
    pass
else:               
    predicted_loss_of_target_TFBS_list = []
    predicted_exclusive_loss_of_target_TFBS_list = []
    predicted_loss_with_regain_of_different_TFBS_for_same_TF_list = []
    predicted_loss_with_gain_of_different_TFBS_list = []

    for sample in interpretation_dict:
        for allele in interpretation_dict.get(sample):
            for instance in interpretation_dict.get(sample).get(allele):
                # loss with corresponding "re-gain" of 'replacement' TFBS for TF of interest
                if instance[21] == 'predicted lost-regained TFBS pair':
                    predicted_loss_with_regain_of_different_TFBS_for_same_TF_list.append(instance)
                # loss without corresponding gain of 'replacement' TFBS for TF of interest
                elif instance[21] == 'predicted TFBS loss (TFBS lost in allele)':
                    predicted_loss_of_target_TFBS_list.append((sample,allele,instance))
                    exclusive_loss_check = True
                # further filter, for loss without corresponding gain of 'replacement' TFBS for TF of interest & no predicted positionally coinciding novel TFBS
                    # consider two coordinate ranges to check (regarding losses):
                    # first, check for allele's span between alignment blocks
                    span_between_alignment_blocks = allele_TFBS_synopsis_df_coordinates_updated.loc[
                        (allele_TFBS_synopsis_df_coordinates_updated['sample'] == sample) & 
                        (allele_TFBS_synopsis_df_coordinates_updated['allele rank'] == allele) &
                        (allele_TFBS_synopsis_df_coordinates_updated['strand'] == instance[12]) &
                        (allele_TFBS_synopsis_df_coordinates_updated['TF'].str.match(r'\b'+TF_of_interest+r'\b')) &
                        (allele_TFBS_synopsis_df_coordinates_updated['Lost TFBS start coordinate (in reference)'] == int(instance[15])) &
                        (allele_TFBS_synopsis_df_coordinates_updated['Lost TFBS end coordinate (in reference)'] == int(instance[16]))]['span between alignment blocks'].values[0]
                    if span_between_alignment_blocks != 'n/a':
                        for index, row in interpreted_TFBS_synopsis_df_updated.iterrows():
                            if row['sample'] == sample and row['allele rank'] == allele and row['interpretation'] == 'predicted TFBS gain (novel to allele)':
                                if set(range(int(row['Gained TFBS coordinate start (in allele)']),int(row['Gained TFBS coordinate end (in allele)']))).intersection(range(span_between_alignment_blocks[0],span_between_alignment_blocks[1])):
                                    predicted_loss_with_gain_of_different_TFBS_list.append(((sample,allele,instance, row.to_list())))
                                    exclusive_loss_check = False
                    else:
                        coordinate_range_to_check = range(int(interpreted_TFBS_synopsis_df_updated.loc[(interpreted_TFBS_synopsis_df_updated['sample'] == sample) &
                                                 (interpreted_TFBS_synopsis_df_updated['allele rank'] == allele) &
                                                  (interpreted_TFBS_synopsis_df_updated['strand'] == instance[12]) &
                                                 (interpreted_TFBS_synopsis_df_updated['interpretation'] == 'predicted TFBS loss (TFBS lost in allele)') &
                                                (interpreted_TFBS_synopsis_df_updated['TF'].str.match(r'\b'+TF_of_interest+r'\b'))]['Lost TFBS coordinate start (in reference)'].values[0]),
                                                  int(interpreted_TFBS_synopsis_df_updated.loc[(interpreted_TFBS_synopsis_df_updated['sample'] == sample) &
                                                 (interpreted_TFBS_synopsis_df_updated['allele rank'] == allele) &
                                                  (interpreted_TFBS_synopsis_df_updated['strand'] == instance[12]) &
                                                 (interpreted_TFBS_synopsis_df_updated['interpretation'] == 'predicted TFBS loss (TFBS lost in allele)') &
                                                (interpreted_TFBS_synopsis_df_updated['TF'].str.match(r'\b'+TF_of_interest+r'\b'))]['Lost TFBS coordinate end (in reference)'].values[0]))
                        for index, row in interpreted_TFBS_synopsis_df_updated.iterrows():
                            if row['sample'] == sample and row['allele rank'] == allele and row['interpretation'] == 'predicted TFBS gain (novel to allele)':
                                if set(coordinate_range_to_check).intersection(range(int(row['Gained TFBS coordinate start (in allele)']), 
                                                                             int(row['Gained TFBS coordinate end (in allele)']))):
                                    predicted_loss_with_gain_of_different_TFBS_list.append(((sample,allele,instance, row.to_list())))
                                    exclusive_loss_check = False
                    if exclusive_loss_check == True:
                        predicted_exclusive_loss_of_target_TFBS_list.append((sample,allele,instance))
                                                  
# Prepare output to indicate **loss of TFBS for TF of interest**, without regain of different TFBS for same TF
if TF_of_interest == '':
    pass
else:
    sample_list = []
    allele_rank_list = []
    allele_list = []
    read_count_list = []
    total_reads_count_list = []
    pct_total_reads_list = []
    pct_reads_filtered_for_1pct_list = []
    pct_reads_filtered_for_10pct_list = []
    allele_sequence_list = []
    reference_sequence_list = []
    alignment_midline_list = []
    TF_lost_list = []
    TF_lost_strand_list = []
    lost_TFBS_sequence_list = []
    lost_TFBS_start_coordinate_list = []
    lost_TFBS_end_coordinate_list = []
    lost_TFBS_pval_list = []
    for pair in predicted_loss_of_target_TFBS_list:
        sample_list.append(pair[0])
        allele_rank_list.append(pair[1])
        allele_list.append(pair[2][2])
        read_count_list.append(pair[2][2].split('_')[2].strip('[]').split('/')[0])
        total_reads_count_list.append(pair[2][2].split('_')[2].strip('[]').split('/')[1])
        pct_total_reads_list.append(pair[2][5])
        pct_reads_filtered_for_1pct_list.append(pair[2][6])
        pct_reads_filtered_for_10pct_list.append(pair[2][7])
        allele_sequence_list.append(pair[2][8])
        reference_sequence_list.append(pair[2][10])
        alignment_midline_list.append(pair[2][9])
        TF_lost_list.append(pair[2][11])
        TF_lost_strand_list.append(pair[2][12])
        lost_TFBS_sequence_list.append(pair[2][13])
        lost_TFBS_start_coordinate_list.append(pair[2][15])
        lost_TFBS_end_coordinate_list.append(pair[2][16])
        lost_TFBS_pval_list.append(pair[2][19])
        
# create dataframe for predicted_loss_of_TFBS_synopsis
if TF_of_interest == '':
    pass
else:
    predicted_loss_of_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list, 
                                        'read count':read_count_list,
                                        'total reads':total_reads_count_list,
                                        '% total reads':pct_total_reads_list,
                                        '% reads filtered for reads <1%':pct_reads_filtered_for_1pct_list,
                                        '% reads filtered for reads <10%':pct_reads_filtered_for_10pct_list,
                                        "alignment query\n(allele sequence)":allele_sequence_list,
                                        "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                        "TF lost":TF_lost_list,                  
                                        "TF lost strand":TF_lost_strand_list,                            
                                        "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":lost_TFBS_sequence_list,
                                        "Lost TFBS coordinate start (in reference)":lost_TFBS_start_coordinate_list,
                                        "Lost TFBS coordinate end (in reference)":lost_TFBS_end_coordinate_list,
                                        "Lost TFBS p-val (in reference)":lost_TFBS_pval_list}
    predicted_loss_of_TFBS_synopsis_df = pd.DataFrame(predicted_loss_of_TFBS_synopsis_df_columns)
    # Add column with allele comment (comment if appropriate)
    predicted_loss_of_TFBS_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in predicted_loss_of_TFBS_synopsis_df['alignment query\n(allele sequence)'].to_list()]]
    predicted_loss_of_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True])
    
# Prepare output to indicate **loss of TFBS for TF of interest with regain of different TFBS for same TF**
if TF_of_interest == '':
    pass
else:
    sample_list = []
    allele_rank_list = []
    allele_list = []
    read_count_list = []
    total_reads_count_list = []
    pct_total_reads_list = []
    pct_reads_filtered_for_1pct_list = []
    pct_reads_filtered_for_10pct_list = []
    allele_sequence_list = []
    reference_sequence_list = []
    alignment_midline_list = []
    TF_lost_list = []
    TF_lost_strand_list = []
    lost_TFBS_sequence_list = []
    lost_TFBS_start_coordinate_list = []
    lost_TFBS_end_coordinate_list = []
    lost_TFBS_pval_list = []
    TF_gained_list = []
    TF_gained_strand_list = []
    gained_TFBS_sequence_list = []
    gained_TFBS_start_coordinate_list = []
    gained_TFBS_end_coordinate_list = []
    gained_TFBS_pval_list = []
    for allele in predicted_loss_with_regain_of_different_TFBS_for_same_TF_list:
        sample_list.append(allele[0])
        allele_rank_list.append(allele[1])
        allele_list.append(allele[2])
        read_count_list.append(allele[2].split('_')[2].strip('[]').split('/')[0])
        total_reads_count_list.append(allele[2].split('_')[2].strip('[]').split('/')[1])
        pct_total_reads_list.append(allele[5])
        pct_reads_filtered_for_1pct_list.append(allele[6])
        pct_reads_filtered_for_10pct_list.append(allele[7])
        allele_sequence_list.append(allele[8])
        reference_sequence_list.append(allele[10])
        alignment_midline_list.append(allele[9])
        TF_lost_list.append(allele[11])
        TF_lost_strand_list.append(allele[12])
        lost_TFBS_sequence_list.append(allele[13])
        lost_TFBS_start_coordinate_list.append(allele[15])
        lost_TFBS_end_coordinate_list.append(allele[16])
        lost_TFBS_pval_list.append(allele[19])
        TF_gained_list.append(allele[11])
        TF_gained_strand_list.append(allele[12])
        gained_TFBS_sequence_list.append(allele[14])
        gained_TFBS_start_coordinate_list.append(allele[17])
        gained_TFBS_end_coordinate_list.append(allele[18])
        gained_TFBS_pval_list.append(allele[20])
        
# create dataframe for predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis
if TF_of_interest == '':
    pass
else:
    predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list, 
                                        'read count':read_count_list,
                                        'total reads':total_reads_count_list,
                                        '% total reads':pct_total_reads_list,
                                        '% reads filtered for reads <1%':pct_reads_filtered_for_1pct_list,
                                        '% reads filtered for reads <10%':pct_reads_filtered_for_10pct_list,
                                        "alignment query\n(allele sequence)":allele_sequence_list,
                                        "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                        "TF lost":TF_lost_list,
                                        "TF gained":TF_gained_list,                    
                                        "TF lost strand":TF_lost_strand_list,
                                        "TF gained strand":TF_gained_strand_list,                        
                                        "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":lost_TFBS_sequence_list,
                                        "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":gained_TFBS_sequence_list, 
                                        "Lost TFBS coordinate start (in reference)":lost_TFBS_start_coordinate_list,
                                        "Lost TFBS coordinate end (in reference)":lost_TFBS_end_coordinate_list,
                                        "Gained TFBS coordinate start (in allele)":gained_TFBS_start_coordinate_list,
                                        "Gained TFBS coordinate end (in allele)":gained_TFBS_end_coordinate_list,
                                        "Lost TFBS p-val (in reference)":lost_TFBS_pval_list,
                                        "Gained TBFS p-val (in allele)":gained_TFBS_pval_list}
    predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df = pd.DataFrame(predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df_columns)
    # Add column with allele comment (comment if appropriate)
    predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df['alignment query\n(allele sequence)'].to_list()]]
    predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True])

# Prepare output to indicate **loss of TFBS for TF of interest with gain of TFBS for novel TF**
if TF_of_interest == '':
    pass
else:
    sample_list = []
    allele_rank_list = []
    allele_list = []
    read_count_list = []
    total_reads_count_list = []
    pct_total_reads_list = []
    pct_reads_filtered_for_1pct_list = []
    pct_reads_filtered_for_10pct_list = []
    allele_sequence_list = []
    reference_sequence_list = []
    alignment_midline_list = []
    TF_lost_list = []
    TF_lost_strand_list = []
    lost_TFBS_sequence_list = []
    lost_TFBS_start_coordinate_list = []
    lost_TFBS_end_coordinate_list = []
    lost_TFBS_pval_list = []
    TF_gained_list = []
    TF_gained_strand_list = []
    gained_TFBS_sequence_list = []
    gained_TFBS_start_coordinate_list = []
    gained_TFBS_end_coordinate_list = []
    gained_TFBS_pval_list = []
    for pair in predicted_loss_with_gain_of_different_TFBS_list:
        sample_list.append(pair[0])
        allele_rank_list.append(pair[1])
        allele_list.append(pair[2][2])
        read_count_list.append(pair[2][2].split('_')[2].strip('[]').split('/')[0])
        total_reads_count_list.append(pair[2][2].split('_')[2].strip('[]').split('/')[1])
        pct_total_reads_list.append(pair[2][5])
        pct_reads_filtered_for_1pct_list.append(pair[2][6])
        pct_reads_filtered_for_10pct_list.append(pair[2][7])
        allele_sequence_list.append(pair[2][8])
        reference_sequence_list.append(pair[2][10])
        alignment_midline_list.append(pair[2][9])
        TF_lost_list.append(pair[2][11])
        TF_lost_strand_list.append(pair[2][12])
        lost_TFBS_sequence_list.append(pair[2][13])
        lost_TFBS_start_coordinate_list.append(pair[2][15])
        lost_TFBS_end_coordinate_list.append(pair[2][16])
        lost_TFBS_pval_list.append(pair[2][19])
        TF_gained_list.append(pair[3][11])
        TF_gained_strand_list.append(pair[3][12])
        gained_TFBS_sequence_list.append(pair[3][14])
        gained_TFBS_start_coordinate_list.append(pair[3][17])
        gained_TFBS_end_coordinate_list.append(pair[3][18])
        gained_TFBS_pval_list.append(pair[3][20])
        
# create dataframe for predicted_loss_with_gain_of_different_TFBS_synopsis
if TF_of_interest == '':
    pass
else:
    predicted_loss_with_gain_of_different_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list, 
                                        'read count':read_count_list,
                                        'total reads':total_reads_count_list,
                                        '% total reads':pct_total_reads_list,
                                        '% reads filtered for reads <1%':pct_reads_filtered_for_1pct_list,
                                        '% reads filtered for reads <10%':pct_reads_filtered_for_10pct_list,
                                        "alignment query\n(allele sequence)":allele_sequence_list,
                                        "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                        "TF lost":TF_lost_list,
                                        "TF gained":TF_gained_list,                    
                                        "TF lost strand":TF_lost_strand_list,
                                        "TF gained strand":TF_gained_strand_list,                        
                                        "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":lost_TFBS_sequence_list,
                                        "Gained TFBS sequence (not in reference at this position, novel to allele)\n*Note: this TFBS sequence is in the allele, 5'-3' on strand indicated in 'strand'":gained_TFBS_sequence_list, 
                                        "Lost TFBS coordinate start (in reference)":lost_TFBS_start_coordinate_list,
                                        "Lost TFBS coordinate end (in reference)":lost_TFBS_end_coordinate_list,
                                        "Gained TFBS coordinate start (in allele)":gained_TFBS_start_coordinate_list,
                                        "Gained TFBS coordinate end (in allele)":gained_TFBS_end_coordinate_list,
                                        "Lost TFBS p-val (in reference)":lost_TFBS_pval_list,
                                        "Gained TBFS p-val (in allele)":gained_TFBS_pval_list}
    predicted_loss_with_gain_of_different_TFBS_synopsis_df = pd.DataFrame(predicted_loss_with_gain_of_different_TFBS_synopsis_df_columns)
    # Add column with allele comment (comment if appropriate)
    predicted_loss_with_gain_of_different_TFBS_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in predicted_loss_with_gain_of_different_TFBS_synopsis_df['alignment query\n(allele sequence)'].to_list()]]
    predicted_loss_with_gain_of_different_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True])

# Prepare output to indicate **loss of TFBS for TF of interest with no predicted gain of TFBS for novel TF** (at pval threshold)
if TF_of_interest == '':
    pass
else:
    sample_list = []
    allele_rank_list = []
    allele_list = []
    read_count_list = []
    total_reads_count_list = []
    pct_total_reads_list = []
    pct_reads_filtered_for_1pct_list = []
    pct_reads_filtered_for_10pct_list = []
    allele_sequence_list = []
    reference_sequence_list = []
    alignment_midline_list = []
    TF_lost_list = []
    TF_lost_strand_list = []
    lost_TFBS_sequence_list = []
    lost_TFBS_start_coordinate_list = []
    lost_TFBS_end_coordinate_list = []
    lost_TFBS_pval_list = []
    TF_gained_list = []
    TF_gained_strand_list = []
    gained_TFBS_sequence_list = []
    gained_TFBS_start_coordinate_list = []
    gained_TFBS_end_coordinate_list = []
    gained_TFBS_pval_list = []
    for pair in predicted_exclusive_loss_of_target_TFBS_list:
        sample_list.append(pair[0])
        allele_rank_list.append(pair[1])
        allele_list.append(pair[2][2])
        read_count_list.append(pair[2][2].split('_')[2].strip('[]').split('/')[0])
        total_reads_count_list.append(pair[2][2].split('_')[2].strip('[]').split('/')[1])
        pct_total_reads_list.append(pair[2][5])
        pct_reads_filtered_for_1pct_list.append(pair[2][6])
        pct_reads_filtered_for_10pct_list.append(pair[2][7])
        allele_sequence_list.append(pair[2][8])
        reference_sequence_list.append(pair[2][10])
        alignment_midline_list.append(pair[2][9])
        TF_lost_list.append(pair[2][11])
        TF_lost_strand_list.append(pair[2][12])
        lost_TFBS_sequence_list.append(pair[2][13])
        lost_TFBS_start_coordinate_list.append(pair[2][15])
        lost_TFBS_end_coordinate_list.append(pair[2][16])
        lost_TFBS_pval_list.append(pair[2][19])
        
# create dataframe for predicted_exclusive_loss_of_TFBS_synopsis
if TF_of_interest == '':
    pass
else:
    predicted_exclusive_loss_of_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list, 
                                        'read count':read_count_list,
                                        'total reads':total_reads_count_list,
                                        '% total reads':pct_total_reads_list,
                                        '% reads filtered for reads <1%':pct_reads_filtered_for_1pct_list,
                                        '% reads filtered for reads <10%':pct_reads_filtered_for_10pct_list,
                                        "alignment query\n(allele sequence)":allele_sequence_list,
                                        "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                        "TF lost":TF_lost_list,                  
                                        "TF lost strand":TF_lost_strand_list,                            
                                        "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":lost_TFBS_sequence_list,
                                        "Lost TFBS coordinate start (in reference)":lost_TFBS_start_coordinate_list,
                                        "Lost TFBS coordinate end (in reference)":lost_TFBS_end_coordinate_list,
                                        "Lost TFBS p-val (in reference)":lost_TFBS_pval_list}
    predicted_exclusive_loss_of_TFBS_synopsis_df = pd.DataFrame(predicted_exclusive_loss_of_TFBS_synopsis_df_columns)
    # Add column with allele comment (comment if appropriate)
    predicted_exclusive_loss_of_TFBS_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in predicted_exclusive_loss_of_TFBS_synopsis_df['alignment query\n(allele sequence)'].to_list()]]
    predicted_exclusive_loss_of_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True])

# July 2021
# Take stock of samples with alleles having lost TFBS for TF of interest without regain/gain of TFBS for distinct TF
# Re-populate remaining ranked alleles for these samples, to facilitate genotype inference
if TF_of_interest == '':
    pass
else:
    sample_list = []
    allele_rank_list = []
    allele_list = []
    read_count_list = []
    total_reads_list = []
    total_reads_pct_list = []
    total_reads_1pct_list = []
    total_reads_10pct_list = []
    allele_sequence_list = []
    alignment_midline_list = []
    reference_sequence_list = []
    TF_list = []
    TF_exlusively_lost_list = []
    TF_lost_with_regain_of_TFBS_for_same_TF_list = []
    TF_lost_with_gain_of_distinct_TF_list = []
    reference_TFBS_unchanged_list = []
    strand_list = []
    TFBS_sequence_list = []
    allele_start_coordinate_list = []
    allele_stop_coordinate_list = []
    p_val_list = []
    comment_list = []

    for sample in set(predicted_loss_of_TFBS_synopsis_df['sample'].to_list()).union(
        set(predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df['sample'].to_list())).union(
        set(predicted_loss_with_gain_of_different_TFBS_synopsis_df['sample'].to_list())):
        allele_lost_TFBS_list = []
        allele_rank_lost_TFBS_list = []
        allele_rank_other_list = []
        
        for index, row in predicted_exclusive_loss_of_TFBS_synopsis_df.iterrows():
            if row['sample'] == sample:
                sample_list.append(row['sample'])
                allele_rank_list.append(row['allele rank'])
                allele_rank_lost_TFBS_list.append(row['allele rank'])
                allele_list.append(row['allele ID'])
                read_count_list.append(row['read count'])
                total_reads_list.append(row['total reads'])
                total_reads_pct_list.append(row['% total reads'])
                total_reads_1pct_list.append(row['% reads filtered for reads <1%']) 
                total_reads_10pct_list.append(row['% reads filtered for reads <10%'])
                allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                alignment_midline_list.append(row['alignment midline'])
                reference_sequence_list.append(row['alignment hit\n(reference)'])
                TF_list.append(row['TF lost'])
                TF_exlusively_lost_list.append('x')
                TF_lost_with_regain_of_TFBS_for_same_TF_list.append('')
                TF_lost_with_gain_of_distinct_TF_list.append('')
                reference_TFBS_unchanged_list.append('')
                strand_list.append(row['TF lost strand']) 
                TFBS_sequence_list.append(row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"])
                allele_start_coordinate_list.append(row['Lost TFBS coordinate start (in reference)'])
                allele_stop_coordinate_list.append(row['Lost TFBS coordinate end (in reference)'])
                p_val_list.append(row['Lost TFBS p-val (in reference)'])
                comment_list.append(row['comment']) 
                allele_rank_lost_TFBS_list.append(row['allele rank'])
                allele_lost_TFBS_list.append(row['allele ID'])
                
        for index, row in predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df.iterrows():
             if row['sample'] == sample:
                allele_rank_other_list.append(row['allele rank'])
                allele_list.append(row['allele ID'])
                allele_rank_list.append(row['allele rank'])
                sample_list.append(row['sample'])
                read_count_list.append(row['read count'])
                total_reads_list.append(row['total reads'])
                total_reads_pct_list.append(row['% total reads'])
                total_reads_1pct_list.append(row['% reads filtered for reads <1%']) 
                total_reads_10pct_list.append(row['% reads filtered for reads <10%'])
                allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                alignment_midline_list.append(row['alignment midline'])
                reference_sequence_list.append(row['alignment hit\n(reference)'])
                TF_list.append(row['TF lost'])
                TF_exlusively_lost_list.append('')
                TF_lost_with_regain_of_TFBS_for_same_TF_list.append('x')
                TF_lost_with_gain_of_distinct_TF_list.append('')
                reference_TFBS_unchanged_list.append('')
                strand_list.append(row['TF lost strand']) 
                TFBS_sequence_list.append(row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"])
                allele_start_coordinate_list.append(row['Lost TFBS coordinate start (in reference)'])
                allele_stop_coordinate_list.append(row['Lost TFBS coordinate end (in reference)'])
                p_val_list.append(row['Lost TFBS p-val (in reference)'])
                comment_list.append(row['comment'])            
                
        for index, row in predicted_loss_with_gain_of_different_TFBS_synopsis_df.iterrows():
            if row['sample'] == sample:
                allele_rank_other_list.append(row['allele rank'])
                allele_list.append(row['allele ID'])
                allele_rank_list.append(row['allele rank'])
                sample_list.append(row['sample'])
                read_count_list.append(row['read count'])
                total_reads_list.append(row['total reads'])
                total_reads_pct_list.append(row['% total reads'])
                total_reads_1pct_list.append(row['% reads filtered for reads <1%']) 
                total_reads_10pct_list.append(row['% reads filtered for reads <10%'])
                allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                alignment_midline_list.append(row['alignment midline'])
                reference_sequence_list.append(row['alignment hit\n(reference)'])
                TF_list.append(row['TF lost'])
                TF_exlusively_lost_list.append('')
                TF_lost_with_regain_of_TFBS_for_same_TF_list.append('')
                TF_lost_with_gain_of_distinct_TF_list.append('x')
                reference_TFBS_unchanged_list.append('')
                strand_list.append(row['TF lost strand']) 
                TFBS_sequence_list.append(row["Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'"])
                allele_start_coordinate_list.append(row['Lost TFBS coordinate start (in reference)'])
                allele_stop_coordinate_list.append(row['Lost TFBS coordinate end (in reference)'])
                p_val_list.append(row['Lost TFBS p-val (in reference)'])
                comment_list.append(row['comment'])  
                        
        for index, row in allele_TFBS_synopsis_df.iterrows():
            if row['sample'] == sample:
                if row['allele rank'] not in allele_rank_lost_TFBS_list:
                    if row['allele rank'] not in allele_rank_other_list:
                        allele_rank_other_list.append(row['allele rank'])
                        allele_list.append(row['allele ID'])
                        allele_rank_list.append(row['allele rank'])
                        sample_list.append(row['sample'])
                        read_count_list.append(row['reads'])
                        total_reads_list.append(row['total reads'])
                        total_reads_pct_list.append(row['% total reads'])
                        total_reads_1pct_list.append(row['% reads filtered for reads <1%']) 
                        total_reads_10pct_list.append(row['% reads filtered for reads <10%'])
                        allele_sequence_list.append(row['alignment query\n(allele sequence)'])
                        alignment_midline_list.append(row['alignment midline'])
                        reference_sequence_list.append(row['alignment hit\n(reference)'])
                        TF_list.append('n/a')
                        TF_exlusively_lost_list.append('')
                        TF_lost_with_regain_of_TFBS_for_same_TF_list.append('')
                        TF_lost_with_gain_of_distinct_TF_list.append('')
                        reference_TFBS_unchanged_list.append('x')
                        strand_list.append('n/a') 
                        TFBS_sequence_list.append('n/a')
                        allele_start_coordinate_list.append('n/a')
                        allele_stop_coordinate_list.append('n/a')
                        p_val_list.append('n/a')
                        comment_list.append(row['comment'])
    
    samples_predicted_to_have_lost_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list,
                                "read count": read_count_list, "total reads": total_reads_list,
                                "% total reads": total_reads_pct_list, "% reads filtered for reads <1%": total_reads_1pct_list,
                                "% reads filtered for reads <10%": total_reads_10pct_list, "alignment query\n(allele sequence)":allele_sequence_list,
                                "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                "TF lost (lost TFBS; no predicted regain of related for same TF, or gain of novel TFBS for distinct TF)":TF_list,
                                "TFBS for TF exclusively lost": TF_exlusively_lost_list,
                                "TFBS for TF lost with regain of different TFBS for same TF": TF_lost_with_regain_of_TFBS_for_same_TF_list,
                                "TFBS for TF lost with gain of TFBS for different TF": TF_lost_with_gain_of_distinct_TF_list,
                                "TFBS for TF unchanged relative to reference": reference_TFBS_unchanged_list,
                                "TF lost strand":strand_list, 
                                "Lost TFBS sequence (in reference at this position, lost in allele)\n*Note: this TFBS sequence is in the reference, 5'-3' on strand indicated in 'strand'":TFBS_sequence_list,
                                "Lost TFBS coordinate start (in reference)":allele_start_coordinate_list,
                                "Lost TFBS coordinate end (in reference)":allele_stop_coordinate_list,
                                "Lost TFBS p-val (in reference)":p_val_list, "comment":comment_list}

    samples_predicted_to_have_lost_TFBS_synopsis_df = pd.DataFrame(samples_predicted_to_have_lost_TFBS_synopsis_df_columns)
    
    samples_predicted_to_have_lost_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost (lost TFBS; no predicted regain of related for same TF, or gain of novel TFBS for distinct TF)',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True])
    
if TF_of_interest == '':
    pass
else:
    samples_predicted_to_have_lost_TFBS_synopsis_df.drop_duplicates(inplace=True)
    samples_predicted_to_have_lost_TFBS_synopsis_df = samples_predicted_to_have_lost_TFBS_synopsis_df.reset_index(drop=True)

if TF_of_interest == '':
    pass
else:
    genotype_interpretation_dict = {}
    for sample in set(samples_predicted_to_have_lost_TFBS_synopsis_df['sample'].to_list()):
        allele_ranks_read_pct_list = []
        for index, row in samples_predicted_to_have_lost_TFBS_synopsis_df.sort_values(by=['sample','allele rank']).iterrows():
            if row['sample'] == sample:
                allele_ranks_read_pct_list.append((row['allele rank'], row['% reads filtered for reads <10%'], row['TFBS for TF exclusively lost'], row['TFBS for TF lost with regain of different TFBS for same TF'],
                                               row['TFBS for TF lost with gain of TFBS for different TF'], row['TFBS for TF unchanged relative to reference']))
        for index, i in enumerate(sorted(set(allele_ranks_read_pct_list))):
            if sample not in genotype_interpretation_dict:
                if int(i[0]) == 1 and int(i[1]) > 90 and i[2] == 'x':
                    genotype_interpretation_dict[sample] = 'predicted homozygous loss in high-ranking allele (no regain or gain)'
                elif int(i[0]) == 1 and int(i[1]) > 90 and i[3] == 'x':
                    genotype_interpretation_dict[sample] = 'predicted homozygous loss in high-ranking allele, with loss having regained a TFBS for TF'
                elif int(i[0]) == 1 and int(i[1]) > 90 and i[4] == 'x':
                    genotype_interpretation_dict[sample] = 'predicted homozygous loss in high-ranking allele, with loss having gained a novel TFBS for a distinct TF'
                elif int(i[0]) == 1 and int(i[1]) > 90 and i[5] == 'x':
                    genotype_interpretation_dict[sample] = 'predicted loss in inconsequential (low-ranking) allele rank(s)'
 
                elif int(i[0]) == 1 and 35 < int(i[1]) < 90:
                    if i[2] == 'x':
                        if int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][2] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles (no regain or gain)'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][3] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has regained a TFBS for TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][4] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has gained a novel TFBS for a distinct TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][5] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted heterozygous loss among high-ranking alleles'
                        else:
                            genotype_interpretation_dict[sample] = 'predicted loss among high-ranking allele'
                    elif i[3] == 'x':
                        if int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][2] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has regained a TFBS for TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][3] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but both of the 2 losses have regained a TFBS for TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][4] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has gained a novel TFBS for a distinct TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][5] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted heterozygous loss among high-ranking alleles, with single loss having regained a TFBS for TF'
                        else:
                            genotype_interpretation_dict[sample] = 'predicted loss among high-ranking allele, with loss having regained a TFBS for TF'
                    elif i[4] == 'x':
                        if int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][2] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has gained a novel TFBS for a distinct TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][3] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has regained a TFBS for TF and 1 has gained a novel TFBS for a distinct TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][4] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted biallelic loss among high-ranking alleles, but both of the 2 losses have gained a novel TFBS for a distinct TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][5] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted heterozygous loss among high-ranking alleles, with single loss having gained a novel TFBS for a distinct TF'
                        else:
                            genotype_interpretation_dict[sample] = 'predicted loss among high-ranking allele, with loss having gained a novel TFBS for a distinct TF'
                    elif i[5] == 'x':
                        if int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][2] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted heterozygous loss among high-ranking alleles'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][3] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted heterozygous loss among high-ranking alleles, with single loss having regained a TFBS for TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][4] == 'x':
                            genotype_interpretation_dict[sample] = 'predicted heterozygous loss among high-ranking alleles, with single loss having gained a novel TFBS for a distinct TF'
                        elif int(sorted(set(allele_ranks_read_pct_list))[index+1][0]) == 2 and 30 < int(sorted(set(allele_ranks_read_pct_list))[index+1][1]) < 90 and sorted(set(allele_ranks_read_pct_list))[index+1][5] == 'x':
                             genotype_interpretation_dict[sample] = 'predicted loss in inconsequential (low-ranking) allele rank(s)'
                        else:
                             genotype_interpretation_dict[sample] = 'predicted loss in inconsequential (low-ranking) allele rank(s)'  
                else:
                    genotype_interpretation_dict[sample] = 'shrug?'
                
if TF_of_interest == '':
    pass
else:
    genotype_inference_list = []
    for index, row in samples_predicted_to_have_lost_TFBS_synopsis_df.iterrows():
        genotype_inference_list.append(genotype_interpretation_dict.get(row['sample']))
    
    samples_predicted_to_have_lost_TFBS_synopsis_df['genotype inference'] = genotype_inference_list

from pandas.api.types import CategoricalDtype
cat_genotype_order = CategoricalDtype(
    ['predicted loss among high-ranking allele',

'predicted biallelic loss among high-ranking alleles (no regain or gain)',

 'predicted loss among high-ranking allele, with loss having regained a TFBS for TF',

 'predicted loss among high-ranking allele, with loss having gained a novel TFBS for a distinct TF',

 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has regained a TFBS for TF',

 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has gained a novel TFBS for a distinct TF',
     
 'predicted homozygous loss in high-ranking allele, with loss having regained a TFBS for TF',
     
 'predicted homozygous loss in high-ranking allele, with loss having gained a novel TFBS for a distinct TF',
     
 'predicted biallelic loss among high-ranking alleles, but both of the 2 losses have regained a TFBS for TF',

 'predicted biallelic loss among high-ranking alleles, but 1 of the 2 losses has regained a TFBS for TF and 1 has gained a novel TFBS for a distinct TF',

 'predicted biallelic loss among high-ranking alleles, but both of the 2 losses have gained a novel TFBS for a distinct TF',

 'predicted heterozygous loss among high-ranking alleles',

 'predicted heterozygous loss among high-ranking alleles, with single loss having regained a TFBS for TF',

 'predicted heterozygous loss among high-ranking alleles, with single loss having gained a novel TFBS for a distinct TF',

 'predicted loss in inconsequential (low-ranking) allele rank(s)'], 
    ordered=True
)

if TF_of_interest == '':
    pass
else:
    samples_predicted_to_have_lost_TFBS_synopsis_df['genotype inference'] = samples_predicted_to_have_lost_TFBS_synopsis_df['genotype inference'].astype(cat_genotype_order)

# Finally, prepare output that summarizes all TFBS detected for given samples
sample_list = []
allele_rank_list = []
allele_list = []
allele_sequence_list = []
reference_sequence_list = []
alignment_midline_list = []
TF_list = []
strand_list = []
TFBS_sequence_list = []
p_val_list = []
lostvsgained_list = []
allele_start_coordinate_list = []
allele_stop_coordinate_list = []
ref_start_coordinate_list = []
ref_stop_coordinate_list = []
for sample in dict_allele_TFBS_synopsis:
    allele_count = 0
    for allele in dict_allele_TFBS_synopsis.get(sample):
        allele_count = allele_count+1
        for TFBS in dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites'):
            if len(dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites')) == 0:
                sample_list.append(sample)
                allele_rank_list.append(allele_count)
                allele_list.append(allele)
                allele_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].strip())
                reference_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[1].split('>')[1].split('<')[0].strip())
                alignment_midline_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[2].split('>')[1].split('<')[0].strip())   
                TF_list.append(TFBS.split(',')[0])
                strand_list.append(TFBS.split(',')[1])                           
                p_val_list.append(TFBS.split(',')[3])
                allele_start_coordinate_list.append(TFBS.split(',')[4])
                allele_stop_coordinate_list.append(TFBS.split(',')[5])
                TFBS_sequence_list.append(TFBS.split(',')[2])
            else:
                sample_list.append(sample)
                allele_rank_list.append(allele_count)
                allele_list.append(allele)
                allele_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0].strip())
                reference_sequence_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[1].split('>')[1].split('<')[0].strip())
                alignment_midline_list.append(dict_allele_TFBS_synopsis.get(sample).get(allele).get('allele_sequence')[2].split('>')[1].split('<')[0].strip())   
                TF_list.append(TFBS.split(',')[0])
                strand_list.append(TFBS.split(',')[1])
                p_val_list.append(TFBS.split(',')[3])
                allele_start_coordinate_list.append(TFBS.split(',')[4])
                allele_stop_coordinate_list.append(TFBS.split(',')[5])
                TFBS_sequence_list.append(TFBS.split(',')[2])
                
all_TFBS_synopsis_df_columns = {"sample":sample_list, "allele rank":allele_rank_list, "allele ID":allele_list, 
                                "alignment query\n(allele sequence)":allele_sequence_list,
                                "alignment midline":alignment_midline_list, "alignment hit\n(reference)":reference_sequence_list,
                                "TF":TF_list,
                                "strand":strand_list, 
                                "TFBS sequence":TFBS_sequence_list,
                                "TFBS coordinate start (in allele)":allele_start_coordinate_list,
                                "TFBS coordinate end (in allele)":allele_stop_coordinate_list,
                                "TFBS p-val":p_val_list}

all_TFBS_synopsis_df = pd.DataFrame(all_TFBS_synopsis_df_columns)

# Add read count data
read_count_list = [i.split('_')[2].strip('[]').split('/')[0] for i in all_TFBS_synopsis_df['allele ID'].to_list()]
total_reads_list = [i.split('_')[2].strip('[]').split('/')[1] for i in all_TFBS_synopsis_df['allele ID'].to_list()]
pct_total_reads_list = [i.split('_')[4].split(':')[1] for i in all_TFBS_synopsis_df['allele ID'].to_list()]
pct_reads_filtered_for_1pct_list = [float(i.split('_')[7].split(':')[1]) if i.split('_')[7].split(':')[1] != 'None' else 0 for i in all_TFBS_synopsis_df['allele ID'].to_list()]
pct_reads_filtered_for_10pct_list = [float(i.split('_')[8].split(':')[1]) if i.split('_')[8].split(':')[1] != 'None' else 0 for i in all_TFBS_synopsis_df['allele ID'].to_list()]

all_TFBS_synopsis_df.insert(loc=3, column='reads', value=read_count_list)
all_TFBS_synopsis_df.insert(loc=4, column='total reads', value=total_reads_list)
all_TFBS_synopsis_df.insert(loc=5, column='% total reads', value=pct_total_reads_list)
all_TFBS_synopsis_df.insert(loc=6, column='% reads filtered for reads <1%', value=pct_reads_filtered_for_1pct_list)
all_TFBS_synopsis_df.insert(loc=7, column='% reads filtered for reads <10%', value=pct_reads_filtered_for_10pct_list)

# Add column with allele comment (comment if appropriate)
all_TFBS_synopsis_df['comment'] = ['note: inferred allele length <=50 bp; read may be primer dimer; consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs' if i <=50 else '' for i in [len(x) for x in all_TFBS_synopsis_df['alignment query\n(allele sequence)'].to_list()]]

all_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","TFBS coordinate start (in allele)"],ascending=[True, True, True, True, True])

# Print dataframes to output file (Excel)
collatedTFBS_csv_output = Path(str(output_path)+ '/'+processdate+'_collated_TFBS.xlsx')

collatedTFBS_csv_output = Path(str(output_path)+ '/'+processdate+'_collated_TFBS.xlsx')

with pd.ExcelWriter(collatedTFBS_csv_output) as writer:  
    allele_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","lost or gained in allele (relative to ref)?"],ascending=[True, True, True, True, False]).to_excel(writer, sheet_name='1 TFBS, predicted lost, gained', index=False)
    interpreted_TFBS_synopsis_df_updated.sort_values(by=['sample','allele rank','TF',"strand","interpretation"],ascending=[True, True, True, True, False]).to_excel(writer, sheet_name='2 TFBS, lost-regained pairs', index=False)
    if TF_of_interest == '':
        all_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","TFBS coordinate start (in allele)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='3 All TFBS in alleles', index=False)
    # TF_of_interest length is assessed to account for Excel maximum of 31 characters in tab names
    elif len(TF_of_interest) <= 6:
        predicted_loss_of_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='3 '+TF_of_interest+', lost (all)', index=False)
        predicted_exclusive_loss_of_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='4 '+TF_of_interest+', lost (-gain,-regain)', index=False)
        predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='5 '+TF_of_interest+', lost (+regain)', index=False)
        predicted_loss_with_gain_of_different_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='6 '+TF_of_interest+', lost (+gain)', index=False)
        samples_predicted_to_have_lost_TFBS_synopsis_df.sort_values(by=['genotype inference', 'sample','allele rank',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='7 '+TF_of_interest+', curated samples', index=False)
        all_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","TFBS coordinate start (in allele)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='8 All TFBS in alleles', index=False)
    else:
        adjusted_TF_of_interest = TF_of_interest[:7]
        predicted_loss_of_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='3 '+adjusted_TF_of_interest+' lost (all)', index=False)
        predicted_exclusive_loss_of_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='4 '+adjusted_TF_of_interest+' lost (-gain,-regain)', index=False)
        predicted_loss_with_regain_of_new_TFBS_for_same_TF_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='5 '+adjusted_TF_of_interest+' lost (+regain)', index=False)
        predicted_loss_with_gain_of_different_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF lost',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='6 '+adjusted_TF_of_interest+' lost (+gain)', index=False)
        samples_predicted_to_have_lost_TFBS_synopsis_df.sort_values(by=['genotype inference', 'sample','allele rank',"TF lost strand","Lost TFBS coordinate start (in reference)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='7 '+adjusted_TF_of_interest+' curated samples', index=False)
        all_TFBS_synopsis_df.sort_values(by=['sample','allele rank','TF',"strand","TFBS coordinate start (in allele)"],ascending=[True, True, True, True, True]).to_excel(writer, sheet_name='8 All TFBS in alleles', index=False)

# Relate allele definition (alignments in alignmentoutput_dict2) to TFBS collation for each allele of each sample (with focus on lost and gained TFBS for each allele, relative to reference)
collatedTFBS_output = Path(str(output_path)+ '/'+processdate+'_collated_TFBS.txt')
with open(str(collatedTFBS_output), 'a+') as f:
    print('CollatedMotifs.py: Summary of matches to TFBS motifs detected in sample sequence(s) relative to reference\nDate: ' + (datetime.today().strftime("%m/%d/%Y")) + '\n\n', file = f)
    for i in sorted(dict_allele_TFBS_synopsis):
        print((len(i)*'=')+'\n'+i+'\n'+(len(i)*'='), file = f)
        for allele in sorted(dict_allele_TFBS_synopsis.get(i), key=lambda x: x.split('_')[3]):
            for x in range(0, len(alignmentoutput_dict2.get(i))):
                if alignmentoutput_dict2.get(i)[x][1].split('>')[1].split('<')[0] == allele:
                    test = alignmentoutput_dict2.get(i)[x]
            sum_gained_motifs = []
            sum_lost_motifs = []
            sum_motifs = []
            sum_TFs = []
            TFs_gt1 = []
            lost_motifs_plus_strand = []
            lost_motifs_minus_strand = []
            total_lost_motifs = []
            total_lost_motifs_list = []
            gained_motifs_plus_strand = []
            gained_motifs_minus_strand = []
            total_gained_motifs = []
            sum_motifs = str(len(dict_allele_TFBS_synopsis.get(i).get(allele).get('all_sites')))
            #print(allele+' '+sum_motifs)
            sum_TFs = str(len(dict_allele_TFBS_synopsis.get(i).get(allele).get('TFs')))
            TFs_gt1 = str(len([TF for TF in dict_allele_TFBS_synopsis.get(i).get(allele).get('TFs') if dict_allele_TFBS_synopsis.get(i).get(allele).get('TFs').get(TF) > 1]))
            sum_lost_motifs = str(len(dict_allele_TFBS_synopsis.get(i).get(allele).get('lost')))
            sum_gained_motifs = str(len(dict_allele_TFBS_synopsis.get(i).get(allele).get('gained')))
            # lost
            lost_motifs_plus_strand = [motif.split(' ')[0] for motif in dict_allele_TFBS_synopsis.get(i).get(allele).get('lost') if motif.split(',')[1] == '+']
            lost_motifs_minus_strand = [motif.split(' ')[0] for motif in dict_allele_TFBS_synopsis.get(i).get(allele).get('lost') if motif.split(',')[1] == '-']
            total_lost_motifs = lost_motifs_plus_strand+lost_motifs_minus_strand
            total_lost_motifs_dict = dict(Counter(total_lost_motifs))
            total_lost_motifs_list = [i+':'+str(total_lost_motifs_dict.get(i)) for i in total_lost_motifs_dict]
            total_lost_motifs_list = sorted(total_lost_motifs_list)
            # gained
            gained_motifs_plus_strand = [motif.split(' ')[0] for motif in dict_allele_TFBS_synopsis.get(i).get(allele).get('gained') if motif.split(',')[1] == '+']
            gained_motifs_minus_strand = [motif.split(' ')[0] for motif in dict_allele_TFBS_synopsis.get(i).get(allele).get('gained') if motif.split(',')[1] == '-']
            total_gained_motifs = gained_motifs_plus_strand+gained_motifs_minus_strand
            total_gained_motifs_dict = dict(Counter(total_gained_motifs))
            total_gained_motifs_list = [i+':'+str(total_gained_motifs_dict.get(i)) for i in total_gained_motifs_dict]
            total_gained_motifs_list = sorted(total_gained_motifs_list)
            print(3*' '+'Allele: '+allele.replace('_',' | ')+'\n   Motifs: total distinct sites |'+sum_motifs+'|, total unique TFs |'+sum_TFs+'| (motifs for '+TFs_gt1+' TFs occur >1x)', file = f)
            print(' Synopsis: relative to reference sequence--# lost sites |'+sum_lost_motifs+'|, # new sites |'+sum_gained_motifs+'|', file = f)
            print('  Details: lost |'+str(total_lost_motifs_list).strip('[]').replace("'","")+'|', file = f)
            print('            new |'+str(total_gained_motifs_list).strip('[]').replace("'","")+'|', file = f)
            if len(dict_allele_TFBS_synopsis.get(i).get(allele).get('allele_sequence')[0].split('>')[1].split('<')[0]) <= 50:
                print(5*' '+'Note: inferred allele length <=50 bp; read may be primer dimer;\n'+11*' '+'consult fasta file for this inferred allele, and/or consider pre-processing fastq file (filter reads) prior to running CollatedMotifs', file = f)
            # prepare complete visual mapping of new motifs above allele sequence
            if int(sum_gained_motifs) > 0: 
                print('\n'+11*' '+'NEW motifs:', file = f)
                motif_plus_tracker = []
                motif_minus_tracker = []
                new_motif_plus_list = []
                new_motif_minus_list = []
                for motif in dict_allele_TFBS_synopsis.get(i).get(allele).get('gained'):
                    if motif.split(',')[1] == '+':
                        new_motif_plus_list.append(motif)
                    elif motif.split(',')[1] == '-':
                        new_motif_minus_list.append(motif)
                for new_motif_plus in new_motif_plus_list:
                    if len(motif_plus_tracker) == 0:
                        print(11*' '+'plus(+) strand:', file = f)
                        motif_plus_tracker.append('check')
                    else:
                        pass
                    if re.search(new_motif_plus.split(',')[2],test[7].split('>')[1].split('<')[0]):
                        match = re.search(new_motif_plus.split(',')[2],test[7].split('>')[1].split('<')[0])
                        distance = match.span()[0]
                        sequence = match.group()
                        print((11+distance)*' '+sequence+' |-- '+new_motif_plus.split(',')[0]+' (pval '+new_motif_plus.split(',')[3]+')', file = f)
                    elif re.search(new_motif_plus.split(',')[2],test[7].split('>')[1].split('<')[0].replace('-','')):
                        match = re.search(new_motif_plus.split(',')[2],test[7].split('>')[1].split('<')[0].replace('-',''))
                        distance = match.span()[0]
                        sequence = match.group()
                        print((11+distance)*' '+sequence+' |-- '+new_motif_plus.split(',')[0]+' (pval '+new_motif_plus.split(',')[3]+')'+' [note, approx. position]', file = f)
                for new_motif_minus in new_motif_minus_list:
                    if len(motif_minus_tracker) == 0:
                        print(11*' '+'minus(-) strand:', file = f)
                        motif_minus_tracker.append('check')
                    else:
                        pass
                    seq_revcomp = ''.join(reversed(''.join(nt_dict.get(nt) for nt in test[7].split('>')[1].split('<')[0])))
                    if re.search(new_motif_minus.split(',')[2],seq_revcomp):
                        match = re.search(new_motif_minus.split(',')[2],seq_revcomp)
                        distance = len(seq_revcomp)-match.span()[0]-len(new_motif_minus.split(',')[2])
                        sequence = ''.join(reversed(match.group()))
                        print((11*' '+distance*' '+sequence+' |-- '+new_motif_minus.split(',')[0]+' (pval '+new_motif_minus.split(',')[3]+')'), file = f)
                    elif re.search(new_motif_minus.split(',')[2],seq_revcomp.replace('-','')):
                        match = re.search(new_motif_minus.split(',')[2],seq_revcomp.replace('-',''))
                        distance = len(seq_revcomp)-match.span()[0]-len(new_motif_minus.split(',')[2])
                        sequence = ''.join(reversed(match.group()))
                        print((11*' '+distance*' '+sequence+' |-- '+new_motif_minus.split(',')[0]+' (pval '+new_motif_minus.split(',')[3]+')'+' [note, approx. position]'), file = f)  
            else:
                pass
            print('\n'+4*' '+'query  '+test[7].split('>')[1].split('<')[0]+'\n'+11*' '+test[9].split('>')[1].split('<')[0]+'\n'+' reference '+test[8].split('>')[1].split('<')[0]+'\n', file = f)
            if int(sum_lost_motifs) > 0:
                print(11*' '+'LOST motifs:', file = f)
                motif_plus_tracker = []
                motif_minus_tracker = []
                lost_motif_plus_list = []
                lost_motif_minus_list = []
                for motif in dict_allele_TFBS_synopsis.get(i).get(allele).get('lost'):
                    if motif.split(',')[1] == '+':
                        lost_motif_plus_list.append(motif)
                    elif motif.split(',')[1] == '-':
                        lost_motif_minus_list.append(motif)
                for lost_motif_plus in lost_motif_plus_list:
                    if len(motif_plus_tracker) == 0:
                        print(11*' '+'plus(+) strand:', file = f)
                        motif_plus_tracker.append('check')
                    else:
                        pass
                    if re.search(lost_motif_plus.split(',')[2],test[8].split('>')[1].split('<')[0]):
                        match = re.search(lost_motif_plus.split(',')[2],test[8].split('>')[1].split('<')[0])
                        distance = match.span()[0]
                        sequence = match.group()
                        print((11+distance)*' '+sequence+' |-- '+lost_motif_plus.split(',')[0]+' (pval '+lost_motif_plus.split(',')[3]+')', file = f)
                for lost_motif_minus in lost_motif_minus_list:
                    if len(motif_minus_tracker) == 0:
                        print(11*' '+'minus(-) strand:', file = f)
                        motif_minus_tracker.append('check')
                    else:
                        pass
                    seq_revcomp = ''.join(reversed(''.join(nt_dict.get(nt) for nt in test[8].split('>')[1].split('<')[0])))
                    if re.search(lost_motif_minus.split(',')[2],seq_revcomp):
                        match = re.search(lost_motif_minus.split(',')[2],seq_revcomp)
                        distance = len(seq_revcomp)-match.span()[0]-len(lost_motif_minus.split(',')[2])
                        sequence = ''.join(reversed(match.group()))
                        print((11*' '+distance*' '+sequence+' |-- '+lost_motif_minus.split(',')[0]+' (pval '+lost_motif_minus.split(',')[3]+')'), file = f)
                print('', file = f)
            else:
                pass
                print('\n', file = f)
                
# Log TFBS collation operations time duration
TFBScollationDuration = str(datetime.now() - startTime_TFBScollation).split(':')[0]+' hr|'+str(datetime.now() - startTime_TFBScollation).split(':')[1]+' min|'+str(datetime.now() - startTime_TFBScollation).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_TFBScollation).split(':')[2].split('.')[1]+' microsec'
            
# Assess files in output directory
file_set = [file for file in os.listdir(output_directory) if Path(file).suffix in ('.txt','.fa','.xls','.xlsx')] 

# Assign script end time
endTime = datetime.now()
endTimestr = str(endTime).split(' ')[1].split('.')[0]

# Log entire script operations time duration
processingDuration = str(datetime.now() - startTime).split(':')[0]+' hr|'+str(datetime.now() - startTime).split(':')[1]+' min|'+str(datetime.now() - startTime).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime).split(':')[2].split('.')[1]+' microsec'
   
filename = Path(str(output_path)+ '/'+processdate+'_script_metrics.txt')
with open(filename, 'a') as f:
    print("""File output information:
    Output directory: """ + str(output_directory) +
'\n    Total file #: ' + str(len(file_set)) +
'\n    Total file output sizes: '+path_size(str(output_directory)), file = f)
    for file in file_set:
        print('    '+file+': '+path_size(str(output_directory)+'/'+file), file = f)
    print("""\n\nScript operation times:
    start time: """+startTimestr+
    '\n    makeblastdb and fasta-get-markov processing time: '+makeblastdb_fastagetmarkov_operationsDuration+
    '\n    fasta processing time: '+readcountDuration+
    '\n    alignments processing time: '+alignmentsDuration+
    '\n    allele definitions processing time: '+allele_definitionsDuration+
    '\n    TFBS processing time (FIMO): '+fimoDuration+
    '\n    TFBS collation processing time: '+TFBScollationDuration+
    '\n    total processing time: '+processingDuration+
    '\n    end time: '+endTimestr, file = f)
f.close()
          
# End of script operations
print("Script has completed.  Please find output files at "+str(output_directory))
print("\n")

sys.exit(0)

############################################################################# end
    
