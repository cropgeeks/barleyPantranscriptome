#!/usr/bin/perl

# Convert PANNZER output to SQL for database

# Usage:

# ./pannzer2go.pl Hv_Mx.transcriptAnno_sorted.txt > Hv_Mx.transcriptAnno_sorted.sql

use strict;

print "use morexgeneatlas;\n";

open (ANNOT, "$ARGV[0]");

while(<ANNOT>){

    $_ =~ s/\s+$//;

    my($transcript_id, $pannzer_annotation, $go_ids, $go_terms) = split(/\t/, $_, 11);

    if($pannzer_annotation eq ""){
        next;
    }

    my ($gene_id, $rest) = split(/\.\d+$/, $transcript_id, 2);
    

    $pannzer_annotation =~ s/\;$//;
    $go_terms =~ s/\;$//;
    $go_ids =~ s/\;$//;

    $go_terms =~ s/\'/\\\'/g;
    $pannzer_annotation =~ s/\'/\\\'/g;

    print "insert into gene_annotation(transcript_id, gene_id, dataset_name, pannzer_annotation, go_ids, go_terms) values (\"$transcript_id\", \"$gene_id\", \"HvMxRTD\", \"$pannzer_annotation\", \"$go_ids\", \"$go_terms\");\n";

}