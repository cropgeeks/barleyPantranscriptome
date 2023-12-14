#!/usr/bin/perl

use CGI;
use GD;
use strict;

use lib '/var/www/html/morexgeneatlas';
use morexgeneatlas_wrapper;
use morexgeneatlas;

my $cgi_query = CGI->new();
my $keywords    = $cgi_query->param("keywords");
my $separator   = $cgi_query->param("sep");
my $neg         = $cgi_query->param("neg");
my $dataset     = $cgi_query->param("dataset");
my $output_file = $cgi_query->param("output");

print "Content-type: text/html\n\n";

morexgeneatlas::printHeader("Keyword Search Results");


my $query = morexgeneatlas_wrapper->new('morexgeneatlas');

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
                  <li class='nav-item'><a class='nav-link active' href='keyword.html'>Annotation Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='download.html'>Bulk Data Download</a></li>
                  <li class='nav-item'><a class='nav-link' href='relator.html'>Transcriptome Relationships</a></li>
                </div>
                <li class='nav-item nav-link'  id='nav-title'></li>
                <li class='nav-item nav-link'  id='nav-header'><b>Links</b></li>
                <li class='nav-item nav-link'  id='nav-title'></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://ics.hutton.ac.uk/barleyrtd'>BarleyRTD Website</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://ics.hutton.ac.uk'>Information and Computing Sciences \@Hutton</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://www.hutton.ac.uk'>The James Hutton Institute</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='https://www.barleyhub.org/'>International Barley Hub</a></li>
                <li class='nav-item nav-link'  id='nav-title'></li>

            </ul>
          </div>
        </nav>
        <main role='main' class='col-md-9 ml-sm-auto col-lg-10 pt-3 px-4'>\n\n";



if (($keywords eq "") || ($keywords =~ /^\s+$/)){
  # if nothing is in the text box
  print "No keywords submitted<BR>\n";
  morexgeneatlas::printFooter();
  die("No keywords given");

}

my $and_chosen = "";
my $or_chosen = "";
if ($separator eq "AND") {
  $and_chosen = "checked";
} elsif ($separator eq "OR"){
  $or_chosen = "checked";
}


my @keywords = split(/\s+/, $keywords);
my $n = 0;
my @clean_keywords;
foreach my $word(@keywords) {

  if ($word ne "") {
    $word =~ s/\W+//g;
    $clean_keywords[$n] = "\%$word\%";
    $n++;
  }
  
}
my $clean_keywords = join(" $separator ", @clean_keywords);


my @clean_negs;

if ($neg ne "") {

  my @neg_keywords = split(/\s+/, $neg);
  my $n = 0;

  foreach my $word(@neg_keywords) {
    if ($word ne "") {
      $word =~ s/\W+//g;
      $clean_negs[$n] = "\%$word\%";
      $n++;
    }
  }
}

my $go_phrase = "(pannzer_annotation like '$clean_keywords[0]'";

for(my $i = 1; $i < @clean_keywords; $i++) {

  $go_phrase .= " $separator pannzer_annotation like '$clean_keywords[$i]'";

}
$go_phrase = $go_phrase . ") ";



foreach my $neg(@clean_negs){

  $go_phrase .= " AND pannzer_annotation not like '$neg'";
  
}


my $blast_phrase = "(description like '$clean_keywords[0]'";

for(my $i = 1; $i < @clean_keywords; $i++) {

  $blast_phrase .= " $separator description like '$clean_keywords[$i]'";

}
$blast_phrase = $blast_phrase . ") ";



foreach my $neg(@clean_negs){

  $blast_phrase .= " AND description not like '$neg'";
  
}


  my $total;

  my %matchingSequences = %{$query->searchAllAnnotation($go_phrase, $blast_phrase)};


  my @trackDetails = keys(%matchingSequences);

  if(scalar(@trackDetails)<1){

    print "<div class='data-header' id='go_ann'>Annotation Search Results</div>\n<div class='data-body'>
      <h6><font color='#cc3300'><b>Sorry there are no transcripts with annotation matching those keywords</b></font>
      </h6></div><br>\n";

  } else {

     my $search_results = "<div id=\"keyword_results\">
  <table class='table-bordered' width='100%'><tbody><th>Transcript ID</th><th>Search Type</th><th>Match</th><th>Match Description</th></tr>\n";
  
    my $count_matches;
  
    foreach my $gene_id(sort (keys %matchingSequences)) {

      foreach my $transcript_number(sort {$a <=> $b} keys %{$matchingSequences{$gene_id}}){

        $count_matches++;

        my $add_colour;
        if(($count_matches % 2) == 0){
          $add_colour = "bgcolor='#dddddd'";
        }

        foreach my $type(sort {$a <=> $b} keys %{$matchingSequences{$gene_id}{$transcript_number}}){
      
          my ($dataset_name, $match_name, $description) = split(/\t/, $matchingSequences{$gene_id}{$transcript_number}{$type}, 3);

          my $transcript_id = $gene_id . "." . $transcript_number;

          my $link;
          if($type eq "TAIRPP10"){
            $link = "<a href='http://www.araport.org/locus/$match_name' target='_blank'>$match_name</a>";
          } elsif ($type eq "RICEPP7"){
            $link = "<a href='http://rice.uga.edu/cgi-bin/gbrowse/rice/?name=$match_name' target='_blank'>$match_name</a>";
          } else {
            $link = $match_name;
          }

          $search_results .= "<tr><td $add_colour><a href='gene.cgi?seq_name=$gene_id&dataset=$dataset_name'>$transcript_id</a></td><td $add_colour>$type</td><td $add_colour>$link</td><td $add_colour>$description</td></tr>\n";

          
        }
      }
    }

    print "<div class='data-header' id='go_ann'>Annotation Search Results</div>\n<div class='data-body'>
    <br><br><div id=\"return\">Searching annotation has returned $count_matches transcripts <br><br>
    $search_results
    </tbody></table></div><br>\n\n";
}


morexgeneatlas::printFooter();
 
 
sub make_data_file{
  
  my($phrase, $dataset) = @_;
  
  my %allmatchingSequences = %{$query->searchGOAnnotation($phrase, $dataset)};


  my @trackDetails = keys(%allmatchingSequences);
  
  my $output_file;
  my $total;

  if(scalar(@trackDetails)<1){
  
  } else {
  
  
    my $random_number = int(rand 1000);
  
    my $datestamp = `date '+%d%m%y_%H%M%S'`;
    chomp $datestamp;

    $output_file = $datestamp . "_" . $random_number . ".txt";
  
    if (-e "$temp_dir/$output_file"){
      
      print "<br><div id='return'><font color='#cc3300'>There has been an error in creating the output file</font></div><br><br>\n";
    
      eorna::printFooter();
      
      die("Error: file '$temp_dir/$output_file' already exists\n");
    }
    
    open (OUT, ">$temp_dir/$output_file");
    
    
    print OUT "gene_id\tchr_id\tgene_start\tgene_end\tdirection\tnumber_of_transcripts\tpannzer_annotation\n";
  
    
  
    foreach my $gene_id(sort (keys %allmatchingSequences)) {
    
      my ($dataset_name, $chr_id, $number_of_transcripts, $gene_start, $gene_end, $strand, $pannzer_annotation) = split(/\t/, $allmatchingSequences{$gene_id});
      
      print OUT "$gene_id\t$chr_id\t$gene_start\t$gene_end\t$strand\t$number_of_transcripts\t$pannzer_annotation\n";
      
      $total++;
    
    }
    
    close OUT;
    
    

  }

  return ($output_file, $total);


}


