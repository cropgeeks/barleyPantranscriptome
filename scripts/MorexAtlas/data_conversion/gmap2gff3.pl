#!/usr/bin/perl

# Usage:
# ./gmap2gff3.pl mloc_transcriptome_gmap_morex.gff

# GFF input created by gmap with parameters set to -f 2, -n 1, --gff3-add-separators 0,--min-trimmed-coverage=0.8,--min-identity=0.9

use strict;


open(GMAP, "$ARGV[0]") || die;

my $outfile = $ARGV[0];
$outfile =~ s/\.gff/\.gff3/;

open(GFF3, ">$outfile");

print GFF3 "##gff-version 3\n";



while (<GMAP>){
    
    $_ =~ s/\s+$//;

    if($_ =~ /gff-version/){
        next;
    }

    $_ =~ s/Full=//g;
    $_ =~ s/Short=/\, /g;

    my($contig_id, $method, $type, $start, $stop, $dot1, $strand, $dot2, $comment) = split(/\t/, $_, 9);

    if($comment =~ /^\"/){
        $comment =~ s/\"ID=/ID=/;
        $comment =~ s/\"$//;
    }

    $comment=~ s/\>//g;

    ($comment, my $note) = split(/\;Dir=/, $comment, 2);


    print GFF3 "$contig_id\t$method\t$type\t$start\t$stop\t.\t$strand\t.\t$comment\n";

}






