#!/usr/bin/perl

use strict;

# Pull all mRNA comment lines from the GFF output from GMAP using:
# grep mRNA *gff > all_gmap_morex_mrna_info.txt
# Convert to tab delimited text for the database

# Usage: 

# ./convert_gmap_info_for_sql.pl > convert_gmap_info_for_sql.txt

open(GMAP, "all_gmap_morex_mrna_info.txt");

my %transcript_datasets;

while(<GMAP>){

    my($contig_id, $method, $type, $start, $stop, $dot1, $strand, $dot2, $comment) = split(/\t/, $_, 9);

    my($source_gmap_file, $chr) = split(/:/, $contig_id, 2);

    my($trans_id, $transcript_name, $parent, $gmap_info) = split(/\;/, $_, 4);

    $transcript_name =~ s/Name=//;

    $transcript_name =~ s/\>//;

    # Skip the HvMxRTD ones - need to use the original gtf file for this
    if($transcript_name =~ /^Hv_Mx/){ 
        next;
    }

    my($gene_name, $rest) = split(/\./, $transcript_name, 2);

    $gmap_info =~ s/\=//g;
    $gmap_info =~ s/[a-zA-Z]+//g;

    my $query_dataset;

    if ($source_gmap_file eq "bartv1_transcriptome_gmap_morex.gff"){
        $query_dataset = "BARTv1";
    }
    elsif ($source_gmap_file eq "bartv2_transcriptome_gmap_morex.gff"){
        $query_dataset = "BARTv2";
    }
    elsif ($source_gmap_file eq "horvu_morexv1_gmap_morex.gff"){
        $query_dataset = "HORVU1";
    }
    elsif ($source_gmap_file eq "horvu_morexv3_gmap_morex.gff"){
        $query_dataset = "HORVU3";
    }
    elsif ($source_gmap_file eq "jloc_transcriptome_gmap_morex.gff"){
        $query_dataset = "JLOC";
    }
    elsif ($source_gmap_file eq "mloc_transcriptome_gmap_morex.gff"){
        $query_dataset = "MLOC";
    }
    elsif ($source_gmap_file eq "morex_gsrtd_transcriptome_gmap_morex.gff"){
        $query_dataset = "MxGSRTD";
    }
    elsif ($source_gmap_file eq "panbart20_transcriptome_gmap_morex.gff"){
        $query_dataset = "PanBart20";
    }


    $transcript_datasets{$source_gmap_file} = $query_dataset;

    my($dir, $coverage, $identity, $matches, $mismatches, $indels, $unk) = split(/\;/, $gmap_info, 7);

    print "$transcript_name\t$gene_name\t$query_dataset\t$source_gmap_file\tMorexGS\t$chr\t$start\t$stop\t$strand\t$coverage\t$identity\t$matches\t$mismatches\t$indels\n";

}

foreach my $dataset(sort keys %transcript_datasets){
    print STDERR "$dataset, $transcript_datasets{$dataset}\n";
}