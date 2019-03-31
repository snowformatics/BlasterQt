---
title: 'Blaster'
tags:
  - sequence
  - blast
  - bowtie
authors:
 - name: Stefanie Lueck
   orcid: 0000-0003-0536-835X
   affiliation: 1
 - name: Dimitar Douchkov
   orcid: 0000-0001-6603-4930
   affiliation: 1
affiliations:
 - name: Leibniz-Institut für Pflanzengenetik und Kulturpflanzenforschung Gatersleben, Stadt Seeland, Sachsen-Anhalt
   index: 1
date: 30 March 2019
bibliography: paper.bib
---

# Summary

Common sequence alignment tools such as Basic Local Alignment Search Tool (BLAST) or Bowtie short sequence aligner are used by scientist on a daily basis. Often those programs are with command line interface which makes them very practical for developers and advance users but hard for researchers who are more familiar with a graphical user interface (GUI). Here we present a tool called “Blaster”, which offers simple and intuitive graphical interface for performing some of the most common tasks in nucleotide and protein sequence work. The software is provided with the following modules:

- Sequence similarity search (blastn, blastp, blastx, tblastn, tblastx, Bowtie) in local databases

- Multiple sequence alignment (CLUSTALW2)

- Sequence extraction from multiple FASTA files by sequence ID

- Sequence converter (reverse complement, reverse, translate)

- Sequence format conversions (FASTA, Multiple FASTA, FASTQ, CSV)

- Local  sequence database management

“Blaster” is available as Microsoft Windows installer and the open source code is accessible on GitHub with MIT License.


# Design and Implementation

The Blaster software was implemented using Python 2.7. and PyQt4 was used for the graphical user interface development. Blast 2.2 (Altschul, S.F. et al. 1990), Bowtie (Langmead B et al.), clustalw2 (Thompson JD1 et al. 1994) are required command line sequence alignment tools which will be automatically installed with the Blaster software.

The executable was created with py2exe-0.6.9. and InnoSetup 5 was used to obtain a Windows installer. Blaster was successfully tested  on Microsoft Windows operation systems XP, 7 and 10.

A help file is available under http://labtools.ipk-gatersleben.de/help/BlasterQt/BlasterQt.html

# Conclusions

We developed a Python based desktop application called Blaster (about 200 downloads per year) with a user-friendly GUI to provide an easy access to the most common sequence alignment tools such as Blast and Bowtie with allowing the use of a custom database. With providing a Windows installer we overcome the hurdle for complicated system setup which requires to install the Python interpreter, Blast, clustalw2 and Bowtie command line tools. The results are provided as text files which can be easily shared and analyzed by standard software packages like Microsoft Excel.

# Acknowledgements

Patrick

-![Fidgit deposited in figshare.](figshare_article.png)

# References
