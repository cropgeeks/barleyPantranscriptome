#!/usr/bin/perl

# Converts the merged HvMxRTD file from GTF to GFF for input to JBrowse

# Usage: 

# ./hvmxmerged_gtf2gff.pl iso_rna_hvmx_merged_edited.gtf > iso_rna_hvmx_merged.gff3


use strict;

print "##gff-version 3\n";

open(GTF, "$ARGV[0]") || die;

my %gene_list;

my %gtf_exons;

my $exon_count;

my @all;

while (<GTF>){
    
    $_ =~ s/\s+$//;

    if($_ =~ /^#/){
        next;
    }

    my($contig_id, $method, $type, $start, $stop, $dot1, $strand, $dot2, $comment) = split(/\t/, $_, 9);

    $comment =~ s/\"//g;

    my ($transcript_id, $gene_id, $rest) = split(/\;/, $comment, 3);

    $transcript_id =~ s/transcript_id //;
    
    $gene_id =~ s/ gene_id //;
    $gene_id =~ s/\;//;
    $gene_id =~ s/\s+//g;

    my $gff_id = $transcript_id . ".path1";
    my $mrna_id = $transcript_id . ".mrna1";

    if(!(defined $gene_list{$transcript_id})){

        $gene_list{$transcript_id} = 1;
        $exon_count = 0;
    }

    $exon_count++;

    my $exon_id = $transcript_id . ".exon." . $exon_count;

    $gtf_exons{$transcript_id}{$start} = "$contig_id\t$exon_count\t$stop\t$strand\t$exon_id";

}


foreach my $transcript_id(sort keys %gtf_exons){

    my $exon_text;
    my @exon_positions;
    my $strand;
    my $contig_id;
    my $start;
    my $stop;

    my ($gene_id, $mrna_number) = split(/\.\d+$/, $transcript_id, 2);

    my $mrna_id = $transcript_id . ".mrna1";

    foreach my $start(sort {$a <=> $b} keys %{$gtf_exons{$transcript_id}}){

        ($contig_id, $exon_count, $stop, $strand, my $exon_id) = split(/\t/, $gtf_exons{$transcript_id}{$start}, 6);

        push(@exon_positions, $start);
        push(@exon_positions, $stop);

        $exon_text .= "$contig_id\thvmx_merged\texon\t$start\t$stop\t.\t$strand\t.\tID=$exon_id;Parent=$transcript_id\n";

    }

    my @gene_locations = sort{$a <=> $b} (@exon_positions);

    my $gene_start = $gene_locations[0];

    my $gene_stop = $gene_locations[@gene_locations-1];

    print "$contig_id\thvmx_merged\ttranscript\t$gene_start\t$gene_stop\t.\t$strand\t.\tID=$transcript_id;Name=$transcript_id\n";

    print $exon_text;


}




