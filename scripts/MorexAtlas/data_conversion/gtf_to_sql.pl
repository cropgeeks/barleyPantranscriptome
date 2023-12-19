#!/usr/bin/perl

# Converts a gtf file into sql for adding to the transcript_structure tables

# Usage: 

#  ./gtf_to_sql.pl iso_rna_hvmx_merged_sorted.gtf > iso_rna_hvmx_merged_sorted_gtf.sql


use strict;

print "use morexgeneatlas;\n";

my $exon_number;
my $prev_transcript_id;

open (GTF, "$ARGV[0]");

while(<GTF>){

    $_ =~ s/\s+$//;

    my($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);

    if ($feature eq "exon"){

      $attribute =~ s/\"//g;

      my ($transcript_id, $gene_id) = split(/\; /, $attribute, 2);
    
      $gene_id =~ s/gene_id //;
      $transcript_id =~ s/transcript_id //;

      if($transcript_id ne $prev_transcript_id){
        $exon_number = 0;
      }

      $exon_number++;

      print "insert into transcript_structures(transcript_id, dataset_name, f_start, f_end, gene_id, chr_id, strand, exon_number) values 
      (\"$transcript_id\", \"HvMxRTD\", \"$start\", \"$end\", \"$gene_id\", \"$chr\", \"$strand\", \"$exon_number\");\n";

      $prev_transcript_id = $transcript_id;


    }

}