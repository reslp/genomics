#!/bin/bash
# script to automatize the different steps for extracting rRNA from (meta)transcriptomes
# written by Philipp Resl, Nov. 2017
# needed software: Trinotate, RNAmmer, GNU parallel, select_transcripts.py, retrieve_lineage_from_taxid.R
# read more on this here: https://reslp.github.io/blog/Extract-rRNA-from-Transcriptomes/

BASEDIR="/home/philipp/data/peltigera/Peltigera_britannica/chloromorph/trinity" #path to folder with trinity assembly
BASENAME="Pbrit_chloro" #used as prefex for filenames

echo "Run RNAmmer"
/usr/local/src/Trinotate-Trinotate-v3.1.0/util/rnammer_support/RnammerTranscriptome.pl  --transcriptome $BASEDIR/Trinity.fasta --path_to_rnammer /usr/local/src/rnammer-1.2/rnammer
echo "...done"


echo "Extract info from GFF file"
#take gff files and extract ids
cat Trinity.fasta.arc.rnammer.gff | awk 'BEGIN {FS="\t"}{print ">"$1}' | uniq > $BASENAME.arc_ids.txt
cat Trinity.fasta.bac.rnammer.gff | awk 'BEGIN {FS="\t"}{print ">"$1}' | uniq > $BASENAME.bac_ids.txt
cat Trinity.fasta.euk.rnammer.gff | awk 'BEGIN {FS="\t"}{print ">"$1}' | uniq > $BASENAME.euk_ids.txt
echo "...done"

#extract fasta sequences for ids
echo "Extract fasta sequences for ids"
select_transcripts.py $BASEDIR/Trinity.fasta $BASENAME.bac_ids.txt normal > $BASENAME.bac.rRNA_transcripts.fasta
select_transcripts.py $BASEDIR/Trinity.fasta $BASENAME.euk_ids.txt normal > $BASENAME.euk.rRNA_transcripts.fasta
select_transcripts.py $BASEDIR/Trinity.fasta $BASENAME.arc_ids.txt normal > $BASENAME.arc.rRNA_transcripts.fasta
echo "...done"

#blastn with taxonid
echo "Blastn sequences against remote nt datase for taxonids"
cat $BASENAME.bac.rRNA_transcripts.fasta | parallel -N 1 -j 6 --recstart ">" --pipe "blastn -db nt -query - -outfmt '6 qseqid staxids sseqid pident qlen length mismatch gapope evalue bitscore' -max_target_seqs 1 -max_hsps 1 -num_threads 4" > $BASENAME.bac.taxids.txt
cat $BASENAME.euk.rRNA_transcripts.fasta | parallel -N 1 -j 6 --recstart ">" --pipe "blastn -db nt -query - -outfmt '6 qseqid staxids sseqid pident qlen length mismatch gapope evalue bitscore' -max_target_seqs 1 -max_hsps 1 -num_threads 4" > $BASENAME.euk.taxids.txt
cat $BASENAME.arc.rRNA_transcripts.fasta | parallel -N 1 -j 6 --recstart ">" --pipe "blastn -db nt -query - -outfmt '6 qseqid staxids sseqid pident qlen length mismatch gapope evalue bitscore' -max_target_seqs 1 -max_hsps 1 -num_threads 4" > $BASENAME.arc.taxids.txt
echo "...done"

#get lineage information for taxonids
echo "Get lineage information for taxonids"
retrieve_lineage_from_taxid.R $BASENAME.bac.taxids.txt > $BASENAME.bac.taxid_and_lineage_info.txt
retrieve_lineage_from_taxid.R $BASENAME.euk.taxids.txt > $BASENAME.euk.taxid_and_lineage_info.txt
retrieve_lineage_from_taxid.R $BASENAME.arc.taxids.txt > $BASENAME.arc.taxid_and_lineage_info.txt
echo "...done"
