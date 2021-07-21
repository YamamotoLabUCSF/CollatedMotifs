** Test files for CollatedMotifs.py **
07-20-2021
This directory contains files compatible with test runs of CollatedMotifs.py


For the following user-entered variables, available test files can be used for the indicated data inputs (underscored). Other inputs must be user-defined (e.g., output directory, paths to executables).

Notes regarding # of fastq files as options for script input:
--fastq_files_subset directory contains R1 and R2 fastq files for 8 samples
--an alternative fastq source directory (fastq_files) is available in the associated Zenodo repository (https://doi.org/10.5281/zenodo.3406861), and contains R1 and R2 fastq files for 384 samples;
--file output sizes and time to complete test runs are relatively small for fastq_files_subset (available for download here in GitHub), large for fastq_files (available for download from Zenodo)

Notes regarding fasta sequence file as Markov background input:
--FKBP5_RefSeq-NG_012645.fa contains the sequence for FKBP5 (RefSeq gene sequence NG_012645) as a background estimate for fimo
--an alternative fasta file containing the entire hg38 human genome sequence is available in the associated Zenodo repository (https://doi.org/10.5281/zenodo.3406861)

===========================================================================

output_directory: ***user-defined (system-specific)***

fastq_directory: SampleTestFiles/fastq_files_subset
---------------------------------------------------

fasta_ref: SampleTestFiles/FKBP5_GOR+86.848kb_KE4.txt
-----------------------------------------------------

blastn_path: ***user-defined (system-specific)***

makeblastdb_path: ***user-defined (system-specific)***

db_prefix: FKBP5_GOR+86.848kb
-----------------------------

fimo_path: ***user-defined (system-specific)***

fimo_motifs_path: SampleTestFiles/JASPAR_CORE_2016_vertebrates.meme
-------------------------------------------------------------------

fasta_get_markov_path: ***user-defined (system-specific)***

markov_background_file: SampleTestFiles/FKBP5_RefSeq-NG_012645.fa
-----------------------------------------------------------------

TF_of_interest: NR3C1
---------------------

