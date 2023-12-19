#!/usr/bin/perl

use strict;
use CGI;
use GD;

use lib '/var/www/html/morexgeneatlas';
use morexgeneatlas;
use morexgeneatlas_wrapper;

my $cgi_query = CGI->new();
my $list      = $cgi_query->param("list");


print "Content-type: text/html\n\n";

morexgeneatlas::printHeader("Transcript Relationships");

my $db_query = morexgeneatlas_wrapper->new('morexgeneatlas');

my $temp_dir = '/var/www/html/morexgeneatlas/temp';


print "<div class='container-fluid'>
      <div class='row'>
        <nav class='col-md-2 d-none d-md-block grey sidebar'>
          <div class='sidebar-sticky'>
            <ul class='nav flex-column'>

                <li class='nav-item nav-link'  id='nav-title'><img class='img-fluid' src='images/morexgeneatlas-logo.png' width=150'></li>

                <li class='nav-item nav-link'  id='nav-header'><b>Utilities Menu</b></li>
                <div class='list-group'>
                  <li class='nav-item'><a class='nav-link' href='index.html'>Home</a></li>
                  <li class='nav-item'><a class='nav-link' href='blast.html'>Homology Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='keyword.html'>Annotation Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='download.html'>Bulk Data Download</a></li>
                  <li class='nav-item'><a class='nav-link active' href='relator.html'>Transcript Relationships</a></li>
                </div>
                <li class='nav-item nav-link'  id='nav-title'></li>
                <li class='nav-item nav-link'  id='nav-header'><b>Links</b></li>
                <li class='nav-item nav-link'  id='nav-title'></li>
                <li class='nav-item'><a class='nav-link ext-links' href='contact.html'>Citation / About Authors</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://ics.hutton.ac.uk/barleyrtd'>BarleyRTD Website</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://ics.hutton.ac.uk'>Information and Computing Sciences \@Hutton</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://www.hutton.ac.uk'>The James Hutton Institute</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='https://www.barleyhub.org/'>International Barley Hub</a></li>
                <li class='nav-item nav-link'  id='nav-title'></li>

            </ul>
          </div>
        </nav>
        <main role='main' class='col-md-9 ml-sm-auto col-lg-10 pt-3 px-4'>\n\n";

if (($list eq "") || ($list =~ /^\s+$/)){
  # if nothing is in the text box
  print "No ID's submited<BR>\n";
  morexgeneatlas::printFooter();
  die("No ID's given");

}

my @list = split(/\n/, $list);

print "<div class='data-header'>Transcriptome Relationships</div>\n<div class='data-body'> <p>These comparisons between the various barley transcriptomes were produced by first using <a href='https://pubmed.ncbi.nlm.nih.gov/27008021/'>GMAP</a> to map the transcriptomes against the <a download href='downloads/210316_Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.zip'>Morex Reference Genome</a> assembly and generate locations. The GMAP mapping was run with parameters restricting the results to the single best match with 90% sequence identity and 80% transcript coverage. Then <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7222033/'>gffcompare</a> was used to compare the GMAP locations of the transcriptomes to the locations of the reference transcript dataset (HvMxRTD annotation). </p>\n";

print "<div id=\"results\">\n";

print "<table class='table-bordered' width='100%'><tbody><tr><th>Query ID</th><th>Related Transcript ID</th><th>Dataset</th><th>Gffcompare Class Code</th></tr>\n";

foreach my $id(@list){

    $id =~ s/\s+//g;

    if($id eq ""){
        next;
    }

    my %related= %{$db_query->getRelatedGenes($id)};

    my @trackDetails = keys(%related);

    my %gene_relation;
    my %related_datasets;



    if(scalar(@trackDetails)<1){

        
        print "<tr><td>$id</td><td>None</td><td>-</td><td>-</td></tr>\n";
        

    } else {

        foreach my $related_transcript(sort keys %related){

            foreach my $transcript_dataset(sort {$b <=> $a}keys %{$related{$related_transcript}}){

                my ($query_gene_id, $reference_tran_id, $reference_gene_id, $reference_dataset, $reference_match_length, $query_tran_exon_number, $query_tran_length, $major_isoform, $gffcompare_class_code) = split(/\t/, $related{$related_transcript}{$transcript_dataset}, 9);

                my $transcript_number = $related_transcript;
                $transcript_number =~ s/$query_gene_id\.//;

                $gene_relation{$query_gene_id}{$transcript_dataset}{$transcript_number} = $gffcompare_class_code;

                $related_datasets{$query_gene_id} = $transcript_dataset;

            }
        }

        foreach my $query_transcript(sort keys %gene_relation) {

            my $trans_list;

            foreach my $transcript_dataset(sort {$a <=> $b} keys %{$gene_relation{$query_transcript}}){

                foreach my $transcript_number(sort {$a <=> $b} keys %{$gene_relation{$query_transcript}{$transcript_dataset}}){

                    my $code_link = &morexgeneatlas::get_code_link($query_transcript, $transcript_number, $gene_relation{$query_transcript}{$transcript_dataset}{$transcript_number});
                    
                    $trans_list .="$code_link ,";

                }

                $trans_list =~ s/,$//;

                print "<tr><td>$id</td><td><a href='search_ids.cgi?seq_name=$query_transcript'>$query_transcript</a></td><td>$related_datasets{$query_transcript}</td><td>$trans_list</td></tr>\n";

            }
            
        }

    }
   
}
print "</table><br><br>\n";
print "</div>\n";

morexgeneatlas::printFooter();


sub gffcompare_class_code {

    my ($code) = @_;

    my $code_def;

    if($code eq "="){
        $code_def = "complete, exact intron chain match";
    } elsif ($code eq "complete"){
        $code_def = "exact intron chain match";
    } elsif ($code eq "c"){
        $code_def = "contained in reference transcript (intron compatible)";
    } elsif ($code eq "k"){
        $code_def = "contains reference transcript (reverse containment)";
    } elsif ($code eq "m"){
        $code_def = "retained intron(s) compared to reference, full intron chain match everywhere else";
    } elsif ($code eq "n"){
        $code_def = "completely overlaps intron from reference transcript, partial or no intron chain match everywhere else";
    } elsif ($code eq "j"){
        $code_def = "multi-exon with at least one junction match";
    } elsif ($code eq "e"){
        $code_def = "single exon that partially covers an intron from reference";
    } elsif ($code eq "o"){
        $code_def = "other same strand overlap with reference exons";
    } elsif ($code eq "s"){
        $code_def = "intron match on the opposite strand (likely a mapping error)";
    } elsif ($code eq "x"){
        $code_def = "exonic overlap on the opposite strand";
    } elsif ($code eq "i"){
        $code_def = "fully contained within a reference intron";
    } elsif ($code eq "y"){
        $code_def = "contains a reference within its intron(s)";
    } elsif ($code eq "p"){
        $code_def = "possible polymerase run-on (close to reference but no overlap)";
    } elsif ($code eq "r"){
        $code_def = "repeat (at least 50% bases are soft-masked)";
    } elsif ($code eq "u"){
        $code_def = "none of the above (unknown, intergenic)";
    }

    return $code_def;
}


sub make_link {

    my ($tran_id) = @_;

    my $link = $tran_id;

    if ($tran_id =~ /^BaRT2v18/) {

    } elsif ($tran_id =~ /^BART1/) {

        $link = "<a href='https://ics.hutton.ac.uk/eorna/transcript.cgi?seq_name=$tran_id&dataset=150831_barley_pseudomolecules'>$tran_id</a>";
    
    
    } elsif ($tran_id =~ /^HORVU[0-9]H/) {

        $link = "<a href='https://ics.hutton.ac.uk/eorna/transcript.cgi?seq_name=HORVU1Hr1G000010.1&dataset=150831_barley_pseudomolecules'>$tran_id</a>";
    
    
    } elsif ($tran_id =~ /^HORVU\./) {
    
    
    } elsif ($tran_id =~ /^JLOC/) {
    
        $link = "<a href='https://ics.hutton.ac.uk/barleyGenes/view_cds_sequence2.cgi?seq_name=$tran_id&dataset=Cuffmerge_version4'>$tran_id</a>";
    
    } elsif ($tran_id =~ /^MLOC/) {

        $link = "<a href='https://ics.hutton.ac.uk/morexGenes/view_cds_sequence2.cgi?seq_name=$tran_id&dataset=assembly3_WGSMorex_rbca.fasta'>$tran_id</a>";
    
    
    } else {
        #not a recognised nomenclature
    }

    return $link;

}
