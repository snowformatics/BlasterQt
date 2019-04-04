# BlasterQt

Common sequence alignment tools such as Basic Local Alignment Search Tool (BLAST) or Bowtie short sequence aligner are used by scientist on a daily basis. Often those programs are with command line interface which makes them very practical for developers and advance users but hard for researchers who are more familiar with a graphical user interface. Here we present a tool called “Blaster”, which offers simple and intuitive graphical interface for performing some of the most common tasks in nucleotide and protein sequence work. The software is provided with the following modules:

·        Sequence similarity search (blastn, blastp, blastx, tblastn, tblastx, Bowtie) in local databases

·        Multiple sequence alignment (CLUSTALW2)

·        Sequence extraction from multiple FASTA files by sequence ID

·        Sequence converter (reverse complement, reverse, translate)

·        Sequence format conversions (FASTA, Multiple FASTA, FASTQ, CSV)

·        Local sequence database management



# Installation

We strongly recommend to install "Blaster" with the Microsoft Windows installer which can be downloaded at sourceforge  https://sourceforge.net/projects/blasterqt/files/BlasterQtv1.2.4.exe/download.

For manual installation, the following dependencies are required:

- Python 2.7 (https://www.python.org/download/releases/2.7/)
- Biopython (pip install biopython)
- PyQt4 (https://sourceforge.net/projects/pyqt/)
- Blast+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.19/)
- Bowtie (https://sourceforge.net/projects/bowtie-bio/)
- clustalw2 (http://www.clustal.org/download/current/)


# Bug reports and tests

Bug reports can be either sent by Email (lueck@ipk-gatersleben.de) or as issue on Github (prefered). If possible and not confidential, upload or link the problematic database or sequence files. 

For testing the general functionality of the software, follow the instructions in testing_guideline.txt file in folder test. All required testing files can be found in the same folder.
