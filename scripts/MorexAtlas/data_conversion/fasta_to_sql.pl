#!/usr/bin/perl


# Converts a fasta file into sql for adding to the transcript_sequences table

# Usage:
# ./fasta_to_sql.pl iso_rna_hvmx_merged.fasta > iso_rna_hvmx_merged_fasta.sql"

use strict;

print "use morexgeneatlas;\n";


$/ = "\n>";

open (FASTA, "$ARGV[0]");
my $count;

while(<FASTA>){

    $_ =~ s/\n>//;

    my($transcript_id, $sequence) = split(/\n/, $_, 2);

    $transcript_id =~ s/>//;

    $sequence =~ s/\s+//g;

    my $seq_length = length($sequence); 

    my($gene_id, $transcript_number) = split(/\.\d+$/, $transcript_id, 2);

    my($chr, $id) = split(/G/, $gene_id, 2);

    $chr =~ s/Hv_Mx_//;

    # extra id for each transcript for the database
    $count++;

    print "insert into transcript_sequences(transcript_id, dataset_name, gene_id, chr_id, transcript_sequence, seq_length, t_id) values 
    (\"$transcript_id\", \"HvMxRTD\", \"$gene_id\", \"$chr\", \"$sequence\", \"$seq_length\", \"$count\");\n";

}