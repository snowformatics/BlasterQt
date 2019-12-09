# BlasterQt

Common sequence alignment tools such as Basic Local Alignment Search Tool (BLAST) or Bowtie short sequence aligner are used by scientist on a daily basis. Often those programs are with command line interface which makes them very practical for developers and advance users but hard for researchers who are more familiar with a graphical user interface. Here we present a tool called “Blaster”, which offers simple and intuitive graphical interface for performing some of the most common tasks in nucleotide and protein sequence work. The software is provided with the following modules:

·        Sequence similarity search (blastn, blastp, blastx, tblastn, tblastx, Bowtie) in local databases

·        Multiple sequence alignment (CLUSTALW2)

·        Sequence extraction from multiple FASTA files by sequence ID

·        Sequence converter (reverse complement, reverse, translate)

·        Sequence format conversions (FASTA, Multiple FASTA, FASTQ, CSV)

·        Local sequence database management



# Installation

We strongly recommend to install "Blaster" with the Microsoft Windows installer which can be downloaded here http://dx.doi.org/10.5447/IPK/2019/23.

For manual installation, the following dependencies are required:

- Python 2.7 (https://www.python.org/download/releases/2.7/)
- Biopython (pip install biopython)
- PyQt4 (https://sourceforge.net/projects/pyqt/)
- Blast+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.19/)
- Bowtie (https://sourceforge.net/projects/bowtie-bio/)
- clustalw2 (http://www.clustal.org/download/current/)

# Manual

http://labtools.ipk-gatersleben.de/help/BlasterQt/BlasterQt.pdf


# Bug reports and tests

Bug reports can be either sent by Email (lueck@ipk-gatersleben.de) or reported as issue on Github (prefered). If possible and not confidential, upload or link the problematic database or sequence files. 

For testing the general functionality of the software, follow the instructions in testing_guideline.txt file in folder test. All required testing files can be found in the same folder.

# Contributions

All kind of contributions (documentation, testing, adding more features, suggestions...) are highly appreciated. 

# History
Blaster v. 1.3

Removed FASTQ conversation feature

Removed editing database names

Updated documentation

Bug fixes

"." or spaces in database names will be replaced with "_" (see https://github.com/snowformatics/BlasterQt/issues/9)


Blaster v. 1.4

Updated documentation

Bug fixes

Added ".fa" extension 

Updated test material



