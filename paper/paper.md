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

Common sequence alignment tools such as Basic Local Alignment Search Tool (BLAST) or Bowtie short sequence aligner are used by scientist on a daily basis. Often those programs are with command line interface which makes them very practical for developers and advance users but hard for researchers who are more familiar with a graphical user interface. Here we present a tool called “Blaster”, which offers simple and intuitive graphical interface for performing some of the most common tasks in nucleotide and protein sequence work. The software is provided with the following modules:

- Sequence similarity search (blastn, blastp, blastx, tblastn, tblastx, Bowtie) in local databases

- Multiple sequence alignment (CLUSTALW2)

- Sequence extraction from multiple FASTA files by sequence ID

- Sequence converter (reverse complement, reverse, translate)

- Sequence format conversions (FASTA, Multiple FASTA, FASTQ, CSV)

- Local  sequence database management

“Blaster” is available as Microsoft Windows installer and the open source code is accessible on GitHub with xx license.


# Design and Implementation

The Blaster software was implemented using Python 2.7. and PyQt4 was used for the graphical user interface development. Blast 2.2 (1), Bowtie (2), clustalw2 (3) are required command line tools which will be automatically installed with the Blaster software.

The executable was created with py2exe-0.6.9. and InnoSetup 5 and runs on Microsoft Windows operation system XP, 7 and 10.

A help file is available under http://labtools.ipk-gatersleben.de/help/BlasterQt/BlasterQt.html

# Conclusions

We developed a Python based desktop application called Blaster (about 200 download per year) with a user-friendly GUI to provide an easy access to the most common sequence alignment tools such as Blast and Bowtie with allowing the use of a custom database. With providing a Windows installer we overcome the hurdle for complicated system setup which requeieres to isntall a Python intepereter, Blast, clustalw and bowtie command line tools. The results are provided as text files which can be easily viwed and analyzed by standard software packages e.g. Excel. 

# Acknowledgements

Patrick

This is a proof of concept integration between a GitHub [@GitHub] repo and figshare [@figshare] in an effort to get a DOI for a GitHub repository. When a repository is tagged for release on GitHub, Fidgit [@Fidgit] will import the release into figshare thus giving the code bundle a DOI. In a somewhat meta fashion, Fidgit is publishing itself to figshare with DOI 'https://doi.org/10.6084/m9.figshare.828487' [@figshare_archive].

-![Fidgit deposited in figshare.](figshare_article.png)

# References
