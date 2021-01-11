#!/usr/local/bin/anaconda3/bin/python3
# Note: edit shebang line above as appropriate for your system
# AUTH: Kirk Ehmsen
# FILE: CollatedMotifs-with_user-set_pval.py
# DATE: 09-04-2018/8-04-2019
# DESC: This script accepts text to standard input, and returns TFBS motifs for samples
# from a demultiplexed NGS fastq dataset.  Provided with a reference sequence, the script returns TFBS motifs
# collated as 'new' or 'lost' relative to the reference sequence.
# CollatedMotifs-with_user-set_pval.py is a variant of CollatedMotifs.py, with an 11th input option
# to specify the p-value threshold delivered to FIMO (CollatedMotifs.py default is 0.0001, 1e-4)
# USAGE: ./CollatedMotifs-with_user-set_pval.py or python3 CollatedMotifs-with_user-set_pval.py
# REPO: https://github.com/YamamotoLabUCSF/CollatedMotifs

#############################################################################
# Background notes:
# =================================================
# This is CollatedMotifs-with_user-set_pval.py v1.0
# =================================================
# https://github.com/YamamotoLabUCSF/CollatedMotifs
# v1.0/Committed 8-04-2019
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
# You will be prompted for the following user-specific information (11 items):
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
#         1 string variable
#           * what is the user-specified p-value threshold?

# Output notes:
# ==============================================
# This script produces 5 output files in the user-specified output directory, plus three directories:
# two directories and subsidiary files created by FIMO (fimo_out and fimo_out_ref) and one directory
# and subsidiary files created by MAKEBLASTDB (alignment_database).
# CollatedMotifs.py output files include: 
#	  1. fasta.fa
#	  2. blastn_alignments.txt (output of BLASTN operation on fasta.fa)
#     3. markov_background.txt (output of FASTA-GET-MARKOV operation on user-supplied fasta reference file)
#     4. collated_TFBS.txt (output of script operation on FIMO-generated .tsv files in fimo_out and fimo_out_ref)
#     5. script_metrics.txt (summary/analysis of script operation metrics [metadata])
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
    # 1-Specify output directory.
    print(r"""
    ---------------------------------------------
    Location of OUTPUT DIRECTORY for output files
    ---------------------------------------------
    
    This script produces 5 output files in the user-specified output directory, plus three directories:
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
            
        5. script_metrics.txt (summary/analysis of script operation metrics [metadata])
  
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
    #11-Specify user-preferred p-value threshold (FIMO)
    print(r"""
    ----------------------------------------------------------------------------------------
    P-VALUE THRESHOLD to be used by FIMO in identification & reporting of TFBS motif matches
    ----------------------------------------------------------------------------------------

    This script uses FIMO from the MEME suite of sequence analysis tools as its basis for determining distinct vs.
    common TFBSs.  FIMO can be provided with a user-specified p-value threshold to adjust stringency of reported
    TFBS matches to TF position frequency matrix/matrices.

    Please indicate the p-value threshold you'd like FIMO to use in its TFBS reporting. CollatedMotifs.py uses
    a threshold of p-value = 0.001 (1e-4) in its default; in this script ('with pval threshold setting'), specify
    your choice of threshold.

    Input the threshold setting in the scientific notation format based on 1e (10 raised to a specified power), such as,
    '1e-3' for 10^-3 (0.001), or '5e-3' for 5x10^-3 (0.005). Do not include quotes flanking your text, or spaces
    between the characters.
    
    When you're done entering the customized p-value threshold, press Enter."""+'\n')
    pval_threshold = input(r"""    -----> Customized p-value threshold setting:  """)

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
    CollatedMotifs-with_user-set_pval.py v1.0
    ==============================================
    https://github.com/YamamotoLabUCSF/CollatedMotifs
    Committed 8-04-2019
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
    
    Please paste a single list of input values directly at the command line prompt, specifying the following 10 values.
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
    11-P-VALUE THRESHOLD to be used by FIMO in identification & reporting of TFBS motif matches
    
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
    pval_threshold = input_list[10].strip()

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
\nDate: """ + (datetime.today().strftime("%m/%d/%Y")) +
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
"\n    markov_background_file: "+ str(markov_background_file) +
"\n    pval_threshold: "+ str(pval_threshold) +
"""\n\nfastq file information:
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
        cluster_sequence_R1_dict[lines_R1[x].split(':')[5]+':'+lines_R1[x].split(':')[6].split(' ')[0]] = lines_R1[x+1].strip('\n')
    #cluster_IDs_list_R1 = [x.split(':')[5]+':'+x.split(':')[6].split(' ')[0] for x in lines_R1[0::4]]
    if Path(R2_file).suffix == ".gz":
        with gzip.open(R2_file, "rt") as f:
            lines_R2 = f.readlines()
    elif Path(R2_file).suffix == ".fastq":
        with open(R2_file, 'r') as f:
            lines_R2 = f.readlines()
    for x in range(0,len(lines_R2),4): 
        cluster_sequence_R2_dict[lines_R2[x].split(':')[5]+':'+lines_R2[x].split(':')[6].split(' ')[0]] = lines_R2[x+1].strip('\n')
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
    counter=Counter(merged_read_list)
    modified_read_list_top5 = []
    for i in counter.most_common(5):
        filtered1 = sum([x for x in counter.values() if x/(sum(counter.values())) > 0.01])
        filtered10 = sum([x for x in counter.values() if x/(sum(counter.values())) > 0.1])
        raw_freq = round((100*i[1]/sum(counter.values())),2)
        modified_read_list_top5.append([i[0], '['+str(i[1])+'/'+str(sum(counter.values()))+']', raw_freq, int(stats.percentileofscore([i for i in counter.values()], i[1], 'rank')), round((100*i[1]/sum([i[1] for i in counter.most_common(5)])),2), round((100*i[1]/filtered1),2) if filtered1 > 0 and raw_freq >= 1 else 'None', round((100*i[1]/filtered10),2) if filtered10 > 0 and raw_freq >= 10 else 'None'])
    with open(str(query_input), 'a+') as file:
        for i in modified_read_list_top5:
              file.write('>'+fastaname[0]+'_'+'R1+R2'+'_'+str(i[1])+'_%totalreads:'+str(i[2])+'_percentile:'+str(i[3])+'_%top5reads:'+str(i[4])+'_%readsfilteredfor1%:'+str(i[5])+'_%readsfilteredfor10%:'+str(i[6])+'\n'+i[0]+'\n')        

# Log read count time duration      
readcountDuration = readcountDuration = str(datetime.now()- startTime_readcount).split(':')[0]+' hr|'+str(datetime.now() - startTime_readcount).split(':')[1]+' min|'+str(datetime.now() - startTime_readcount).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_readcount).split(':')[2].split('.')[1]+' microsec'

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

# Impute genotypes from blastn_alignments.txt file      
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

# Further subset 'No hits found' queries for R1 vs. R2      
no_hits_R1_read_list = []
no_hits_R2_read_list = []
for i in no_hits_list:
    if i.split('_')[1] == 'R1':
        no_hits_R1_read_list.append(i.split('_')[0]+' '+i.split('_')[2]+' '+i.split('_')[3].split(':')[1]+'%')
    elif i.split('_')[1] == 'R2':
        no_hits_R2_read_list.append(i.split('_')[0]+' '+i.split('_')[2]+' '+i.split('_')[3].split(':')[1]+'%')

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
multiple_alignments_list = []
for i in alignments_list3:
    if len(re.findall('<Hit_num>', str(i))) > 1:
        multiple_alignments_list.append(i)

# Identify read IDs with >1 alignment to sequences in reference database
multiple_alignments_readID_list = []
for i in multiple_alignments_list:
    multiple_alignments_readID_list.append(i[1].split('>')[1].split('<')[0])

# Record sample names having reads with >1 alignment to sequences in reference database
multiple_alignments_samplename_list = []
for i in multiple_alignments_readID_list:
    samplename = i.split('_')[0]
    if samplename not in multiple_alignments_samplename_list:
        multiple_alignments_samplename_list.append(samplename)

# Prepare dictionary linking sample names to their reads having >1 alignment to sequences in reference database      
multiple_alignments_dict = {}
for i in multiple_alignments_samplename_list:
    multiple_alignments_dict ["{0}".format(i)] = tuple(x for x in multiple_alignments_list if bool(re.search(i, x[1])))

# Prepare alignment_list4 for reads with exclusively 1 alignment hit in reference database
alignments_list4 = []
for i in alignments_list3:
    if i not in multiple_alignments_list:
        alignments_list4.append(i)

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
cmd_TFBS = str(fimo_path)+' --bfile '+str(markovbackground_output)+' --o '+str(ref_TFBS_output)+' --thresh '+pval_threshold+' '+str(fimo_motifs_path)+' '+str(ref_input)

os.system(cmd_TFBS)

# Alleles: FIMO command (usage: fimo --bfile <background file> <motif file> <sequence file>) 
cmd_TFBS = str(fimo_path)+' --bfile '+str(markovbackground_output)+' --o '+str(ref_TFBS_output)+' --thresh '+pval_threshold+' '+str(allele_TFBS_output)+' '+str(fimo_motifs_path)+' '+str(query_input)

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
        dict_allele_TFBS_synopsis[allele].update({alignmentoutput_dict2.get(allele)[x][1].split(">")[1].split("<")[0]:{'gained':[],'lost':[],'all_sites':[], 'TFs':{}}})

#for sample in dict_allele_TFBS_synopsis:
#    for allele in dict_allele_TFBS_synopsis.get(sample):
#        for motif in dict_allele_TFBS.get(sample).get(allele):
#            dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites').append(motif.split('\t')[1]+' #('+motif.split('\t')[0]+'),'+motif.split('\t')[5]+','+motif.split('\t')[9]+','+motif.split('\t')[7])
            
for sample in dict_allele_TFBS_synopsis:
    for allele in dict_allele_TFBS_synopsis.get(sample):
        for motif in dict_allele_TFBS.get(sample).get(allele):
            dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites').append(motif.split('\t')[1]+' ('+motif.split('\t')[0]+'),'+motif.split('\t')[5]+','+motif.split('\t')[9]+','+motif.split('\t')[7])
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
        dict_ref_TFBS_synopsis.get(ref).get('all_sites').append(motif.split('\t')[1]+' ('+motif.split('\t')[0]+'),'+motif.split('\t')[5]+','+motif.split('\t')[9]+','+motif.split('\t')[7])

# Run comparisons, populating into dict_allele_TFBS_synopsis
ref_options = [ref for ref in dict_ref_TFBS]
for sample in dict_allele_TFBS_synopsis:
    # define reference sequence appropriate to sample
    for ref in ref_options:
        if re.search(ref, sample):
            sample_ref = ref
    for allele in dict_allele_TFBS_synopsis.get(sample):
        for motif in dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites'):
            if motif in dict_ref_TFBS_synopsis.get(sample_ref).get('all_sites'):
                pass
            else:
                dict_allele_TFBS_synopsis.get(sample).get(allele).get('gained').append(motif)
        for motif in dict_ref_TFBS_synopsis.get(sample_ref).get('all_sites'):
            if motif in dict_allele_TFBS_synopsis.get(sample).get(allele).get('all_sites'):
                pass
            else:
                dict_allele_TFBS_synopsis.get(sample).get(allele).get('lost').append(motif)
                
# Relate allele definition (alignments in alignmentoutput_dict2) to TFBS collation for each allele of each sample (with focus on lost and gained TFBS for each allele, relative to reference)
collatedTFBS_output = Path(str(output_path)+ '/'+processdate+'_collated_TFBS.txt')
with open(str(collatedTFBS_output), 'a+') as f:
    print('CollatedMotifs.py: Summary of matches to TFBS motifs detected in sample sequence(s) relative to reference\nDate: ' + (datetime.today().strftime("%m/%d/%Y")) + '\n\n', file = f)
    for i in sorted(dict_allele_TFBS_synopsis):
        print((len(i)*'=')+'\n'+i+'\n'+(len(i)*'='), file = f)
        for allele in dict_allele_TFBS_synopsis.get(i):
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
file_set = [file for file in os.listdir(output_directory) if Path(file).suffix in ('.txt','.fa')] 

# Assign script end time
endTime = datetime.now()
endTimestr = str(endTime).split(' ')[1].split('.')[0]

# Log entire script operations time duration
processingDuration = str(datetime.now() - startTime).split(':')[0]+' hr|'+str(datetime.now() - startTime).split(':')[1]+' min|'+str(datetime.now() - startTime).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime).split(':')[2].split('.')[1]+' microsec'
   
filename = Path(str(output_path)+ '/'+processdate+'_script_metrics.txt')
with open(filename, 'a') as f:
    print("""\n\nFile output information:
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
print("\nScript has completed.  Please find output files at "+str(output_directory))

sys.exit(0)    

############################################################################# end
    
