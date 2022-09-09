genomics scripts
=========

This repository contains several small scripts I use in various steps of my genomics analyses. Mostly they are highly personalized for my own needs and I do not provide warranty for their functionality. In this README you will find a brief description of each script.


#### assembly_filter.py

Filters assemblies based on GC content, proportion of blast hit with particular taxonomic identity and mean read coverage.
Helps to identfy contigs that are potential contaminants. This is somewhat similar to Blobology approaches.

Used non-standard Libraries:
BioPython, scipy, numpy


#### gc_plot.py

Calculated GC content of DNA Sequences (e.g. Assemblies) using a sliding window approach. It creates output as PDF files with graphical representations of the GC content of the sequence. Now also includes gene position information.

Used non-standard Libraries:
BioPython, matplotlib

#### get_longest_isoform.py

Identifies longest splicing isoforms from Trinity asssemblies.

Used non-standard Libraries:
BioPython

#### gff3_to_gtf.py

Convert GFF3 (eg. from MAKER) files to GTF


#### codon_alignment.py

Uses Codon Alignments obtained from set of DNA and Protein Sequences to calculate dN/dS.

Used non-standard Libraries:
BioPython, numpy

#### filter_blast.py

Filters BLAST hits in tabular format.

Used non-standard Libraries:
pandas

#### filter_codeml.py

Filters output from codeml.py according to pvalue and number of sites under selection.


#### ortho_select.py

Uses OrthoMCL or blast all vs. all output to select sequences belonging to orthologous groups from multiple species.

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

#### parse_blastxml.py

This script creates a valid xml blast file (-outfmt 5) from a piped concatenated xml file when blast ist used with parallel.

#### parse_to_tab.py

This script converts a blast xml file (-outfmt 5) to tabular (-outfmt 6)

Used non-standard Libraries:
BioPython

#### blast2go_to_tab.py

Converts the output of a Blast2GO BLAST search (several zipped xml files) to a single blast tabular format (-outfmt 6) file.



#### scrape_cazy.py

Downloads information of characterized CAZymes from cazy.org in tab delimited format for downstream analyses.

### get_lineage.py

Downloads and parses lineage information from NCBI. This was used in early versions of phylociraptor.
