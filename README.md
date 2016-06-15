genomics scripts
=========

This repository contains several small scripts I use in various steps of my genomics analyses. Mostly they are highly personalized for my own needs and I do not provide warranty for their functionality. In this README you will find a brief description of each script.


#### assembly_filter.py

Filters assemblies based on GC content, proportion of blast hit with particular taxonomic identity and mean read coverage.
Helps to identfy contigs that are potential contaminants. This is somewhat similar to Blobology approaches.

Used non-standard Libraries:
BioPython, scipy, numpy


#### calc_gc.py

Calculated GC content of DNA Sequences (e.g. Assemblies) using a sliding window approach. It creates output as PDF files with graphical representations of the GC content of the sequence.

Used non-standard Libraries:
BioPython, matplotlib

#### codon_alignment.py

Uses Codon Alignments obtained from set of DNA and Protein Sequences to calculate dN/dS.

Used non-standard Libraries:
BioPython, numpy

#### filter_blast.py

Filters BLAST hits in tabular format.

Used non-standard Libraries:
pandas

#### ortho_select.py

Uses OrthoMCL output to select sequences belonging to orthologous groups from multiple species.

Used non-standard Libraries:
pandas, BioPython

#### filter_transcripts.py

Reduces FASTA files to IDs from another FASTA file.

Used non-standard Libraries:
BioPython

#### select_cazys.py

Reduced a FASTA file to the IDS from dbCAN output according to a specific group of CAZymes

Used non-standard Libraries:
pandas, BioPython


