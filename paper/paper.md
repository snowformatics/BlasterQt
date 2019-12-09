---
title: '“Blaster” – a graphical user interface for common sequence analysis tools'
tags:
  - sequences
  - fasta
  - clustalw
  - blast
  - bowtie
  - bioinformatic
  - GUI
  - Windows
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
 
Common sequence alignment tools such as Basic Local Alignment Search Tool (BLAST) or Bowtie short sequence aligner are used by scientists on a daily basis. Often those programs are with command-line interface which makes them very practical for developers and advance users but hard for researchers who are more familiar with a graphical user interface (GUI). Here we present a tool called “Blaster”, which offers a simple and intuitive GUI for performing some of the most common tasks in nucleotide and protein sequence work. The software is provided with the following modules:
 
- Sequence similarity search (blastn, blastp, blastx, tblastn, tblastx, Bowtie) in local databases
 
- Multiple sequence alignment (CLUSTALW2)
 
- Sequence extraction from multiple FASTA files by sequence ID
 
- Sequence converter (reverse complement, reverse, translate)
 
- Sequence format conversions (FASTA, Multiple FASTA, CSV)
 
- Local  sequence database management
 
In contrast to other similar tools, designed to work with online databases, “Blaster” is primarily aimed to be used with local custom sequence databases. This allows the users to use their specific sequence databases, as well as to protect the privacy of their data.
 
“Blaster” provides easy installation, it is free (MIT License) and with an open source code.
 
# Design and Implementation
 
The “Blaster” software was implemented with Python 2.7 and the GUI was build with PyQt4. Blast 2.2 [@ALTSCHUL1990403], Bowtie [@Langmead2009], clustalw2 [@Rédei2008] are required command line sequence alignment tools which will be automatically installed with the “Blaster” software.
 
The executable file was created with py2exe-0.6.9. The Microsoft Windows installer was created with InnoSetup 5. “Blaster” was successfully tested for functioning on the Microsoft Windows operating systems XP, 7 and 10.
 
A help file is available under http://labtools.ipk-gatersleben.de/help/BlasterQt/BlasterQt.pdf
 
# Conclusions
 
We have developed a Python-based desktop application called “Blaster” with a user-friendly GUI to provide easy access to the most common sequence alignment tools, such as Blast and Bowtie, for sequence analysis in local custom sequence databases. The easy one-click installer for Microsoft Windows maximally simplifies the complex system setup that includes installation of a Python interpreter, Blast, clustalw2, and Bowtie command line tools. The plain text file format of the results allows an easy import to other software packages for further analysis.
 
# Links
Github repository: https://github.com/snowformatics/BlasterQt
 
Windows Installer: http://dx.doi.org/10.5447/IPK/2019/23
  
# References
