#!/usr/bin/perl

use strict;

# Go through all TMAP files generated by gffcompare in the current directory and parse into tab-delimited text for loading into database

# Usage:

# ./convert_tmap_files_for_mysql.pl > convert_tmap_files_for_mysql.txt

my $files = `ls *cmp*tmap`;

my @tmap_files = split(/\n/, $files);

my %done;

foreach my $file(@tmap_files){

    my ($text, $rest) = split(/\./, $file, 2);

    my ($set, $ref) = split(/\_cmp\_/, $text);

    &reformat_tmap_file($set, $ref, $file);

}


sub reformat_tmap_file{
    my ($set, $ref, $file) = @_;

    my %gmap_matches;

    open(TMAP, "$file") || die ("Cannot find file '$file'\n");

    my $tran_dataset = $set;
    my $ref_dataset = $ref;

    while(<TMAP>){

        $_ =~ s/\s+$//;

        if($_ =~ /^ref_gene_id/){
            next;
        }

        $_ =~ s/\.path1//g;
        $_ =~ s/\.mrna1//g;

        my ($ref_gene_id, $ref_id, $class_code, $qry_gene_id, $qry_id, $num_exons, $FPKM, $TPM, $cov, $len, $major_iso_id, $ref_match_len) = split(/\t/, $_, 12);


        $qry_id =~ s/^\>//;

        my $qry_gene_id = $qry_id;
        $qry_gene_id =~ s/\.\d+$//;
        my $ref_gene_id = $ref_id;
        $ref_gene_id =~ s/\.\d+$//;

        print "$qry_id\t$qry_gene_id\t$tran_dataset\t$ref_id\t$ref_gene_id\t$ref_dataset\t$ref_match_len\t$num_exons\t$len\t$major_iso_id\t$class_code\n";

    }
}




